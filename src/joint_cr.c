// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "joint_cr.h"

#include "body.h"
#include "core.h"
#include "stack_allocator.h"
#include "world.h"

#include <float.h>
#include <math.h>

#define S2_CR_MAX_MEASURE_ITERATIONS 1000

// Stop once the preconditioned residual energy falls this far below its start, about float precision on the
// residual. Past it roundoff drives both terms of the step length toward zero and the ratio can blow up.
#define S2_CR_TOLERANCE 1.0e-12

// Float matches a production solver. Double separates roundoff from truncation.
#define S2_CR_DOUBLE 0

#if S2_CR_DOUBLE
typedef double s2CRReal;
#else
typedef float s2CRReal;
#endif

typedef struct s2CRVector
{
	s2CRReal row[3];
} s2CRVector;

typedef struct s2JointBlock
{
	// Compliance and regularization, the part of each joint's diagonal block not carried by the bodies
	s2CRReal extra[3][3];

	// Inverse of the whole diagonal block, the block Jacobi preconditioner
	s2CRReal inverse[3][3];
} s2JointBlock;

typedef struct s2CRWork
{
	s2World* world;
	const s2JointFactor* factor;

	// Body velocity change
	s2CRReal* vx;
	s2CRReal* vy;
	s2CRReal* w;

	s2JointBlock* blocks;
	s2CRVector* vectors;
	s2CRVector* b;
	s2CRVector* x;
	s2CRVector* z;
	s2CRVector* p;
	s2CRVector* az;
	s2CRVector* ap;
	s2CRVector* q;
	s2CRVector* error;
	s2CRVector* product;
	s2CRVector* exact;
} s2CRWork;

// Tracks the energy norm error against an exact solution while iterating
typedef struct s2CRTracker
{
	s2CRReal exactNorm;
	int iterations1;
	int iterations2;
	float error8;
	float error32;
} s2CRTracker;

static void s2InvertBlock(const s2CRReal B[3][3], int n, s2CRReal inverse[3][3])
{
	for (int r = 0; r < 3; ++r)
	{
		for (int c = 0; c < 3; ++c)
		{
			inverse[r][c] = 0;
		}
	}

	s2CRReal det;
	if (n == 2)
	{
		det = B[0][0] * B[1][1] - B[0][1] * B[1][0];
	}
	else
	{
		det = B[0][0] * (B[1][1] * B[2][2] - B[1][2] * B[2][1]) - B[0][1] * (B[1][0] * B[2][2] - B[1][2] * B[2][0]) +
			  B[0][2] * (B[1][0] * B[2][1] - B[1][1] * B[2][0]);
	}

	s2CRReal scale = 1;
	for (int r = 0; r < n; ++r)
	{
		scale *= B[r][r];
	}

	// A body with fixed rotation leaves the block singular. Fall back to the diagonal.
	if (det <= FLT_EPSILON * scale)
	{
		for (int r = 0; r < n; ++r)
		{
			inverse[r][r] = B[r][r] > 0 ? 1 / B[r][r] : 0;
		}
		return;
	}

	s2CRReal invDet = 1 / det;
	if (n == 2)
	{
		inverse[0][0] = invDet * B[1][1];
		inverse[0][1] = -invDet * B[0][1];
		inverse[1][0] = -invDet * B[1][0];
		inverse[1][1] = invDet * B[0][0];
		return;
	}

	inverse[0][0] = invDet * (B[1][1] * B[2][2] - B[1][2] * B[2][1]);
	inverse[0][1] = invDet * (B[0][2] * B[2][1] - B[0][1] * B[2][2]);
	inverse[0][2] = invDet * (B[0][1] * B[1][2] - B[0][2] * B[1][1]);
	inverse[1][0] = invDet * (B[1][2] * B[2][0] - B[1][0] * B[2][2]);
	inverse[1][1] = invDet * (B[0][0] * B[2][2] - B[0][2] * B[2][0]);
	inverse[1][2] = invDet * (B[0][2] * B[1][0] - B[0][0] * B[1][2]);
	inverse[2][0] = invDet * (B[1][0] * B[2][1] - B[1][1] * B[2][0]);
	inverse[2][1] = invDet * (B[0][1] * B[2][0] - B[0][0] * B[2][1]);
	inverse[2][2] = invDet * (B[0][0] * B[1][1] - B[0][1] * B[1][0]);
}

static s2CRReal s2DotJoint(const s2CRVector* a, const s2CRVector* b, int count)
{
	s2CRReal sum = 0;
	for (int i = 0; i < count; ++i)
	{
		sum += a[i].row[0] * b[i].row[0] + a[i].row[1] * b[i].row[1] + a[i].row[2] * b[i].row[2];
	}
	return sum;
}

static void s2MulBlock(const s2CRReal M[3][3], const s2CRVector* x, s2CRVector* y)
{
	for (int r = 0; r < 3; ++r)
	{
		y->row[r] = M[r][0] * x->row[0] + M[r][1] * x->row[1] + M[r][2] * x->row[2];
	}
}

static void s2Precondition(const s2CRWork* w, const s2CRVector* x, s2CRVector* y)
{
	for (int i = 0; i < w->factor->jointCount; ++i)
	{
		s2MulBlock(w->blocks[i].inverse, x + i, y + i);
	}
}

// y = (D + G) x with D = J M^-1 J^T applied through body velocity changes, as a production solver would
static void s2ApplyJointOperator(const s2CRWork* w, const s2CRVector* x, s2CRVector* y)
{
	const s2Body* bodies = w->world->bodies;
	const s2JointFactor* factor = w->factor;
	s2CRReal* vx = w->vx;
	s2CRReal* vy = w->vy;
	s2CRReal* av = w->w;
	int jointCount = factor->jointCount;

	for (int i = 0; i < jointCount; ++i)
	{
		const s2JointRows* rows = factor->rows + i;
		vx[rows->indexA] = vy[rows->indexA] = av[rows->indexA] = 0;
		vx[rows->indexB] = vy[rows->indexB] = av[rows->indexB] = 0;
	}

	for (int i = 0; i < jointCount; ++i)
	{
		const s2JointRows* rows = factor->rows + i;
		const s2Body* bodyA = bodies + rows->indexA;
		const s2Body* bodyB = bodies + rows->indexB;

		s2CRReal px = 0, py = 0;
		for (int r = 0; r < rows->pointCount; ++r)
		{
			px += x[i].row[r] * rows->axes[r].x;
			py += x[i].row[r] * rows->axes[r].y;
		}
		s2CRReal angularImpulse = rows->rowCount > rows->pointCount ? x[i].row[rows->pointCount] : 0;

		s2CRReal mA = bodyA->invMass, iA = bodyA->invI;
		s2CRReal mB = bodyB->invMass, iB = bodyB->invI;
		s2CRReal crossA = (s2CRReal)rows->rA.x * py - (s2CRReal)rows->rA.y * px;
		s2CRReal crossB = (s2CRReal)rows->rB.x * py - (s2CRReal)rows->rB.y * px;

		vx[rows->indexA] -= mA * px;
		vy[rows->indexA] -= mA * py;
		av[rows->indexA] -= iA * (crossA + angularImpulse);
		vx[rows->indexB] += mB * px;
		vy[rows->indexB] += mB * py;
		av[rows->indexB] += iB * (crossB + angularImpulse);
	}

	for (int i = 0; i < jointCount; ++i)
	{
		const s2JointRows* rows = factor->rows + i;
		int a = rows->indexA, b = rows->indexB;

		// Relative velocity of the anchors, v + cross(w, r)
		s2CRReal dvx = (vx[b] - av[b] * rows->rB.y) - (vx[a] - av[a] * rows->rA.y);
		s2CRReal dvy = (vy[b] + av[b] * rows->rB.x) - (vy[a] + av[a] * rows->rA.x);
		s2CRReal dw = av[b] - av[a];

		s2MulBlock(w->blocks[i].extra, x + i, y + i);
		for (int r = 0; r < rows->rowCount; ++r)
		{
			y[i].row[r] += r < rows->pointCount ? rows->axes[r].x * dvx + rows->axes[r].y * dvy : dw;
		}
	}
}

static s2CRWork s2CreateCRWork(s2World* world, const s2JointFactor* factor)
{
	s2StackAllocator* stack = world->stackAllocator;
	int jointCount = factor->jointCount;
	int bodyCapacity = world->bodyPool.capacity;

	s2CRWork w = {0};
	w.world = world;
	w.factor = factor;
	w.vx = s2AllocateStackItem(stack, 3 * bodyCapacity * sizeof(s2CRReal), "cr body velocity");
	w.vy = w.vx + bodyCapacity;
	w.w = w.vy + bodyCapacity;
	w.blocks = s2AllocateStackItem(stack, jointCount * sizeof(s2JointBlock), "cr blocks");

	int vectorCount = 10;
	w.vectors = s2AllocateStackItem(stack, vectorCount * jointCount * sizeof(s2CRVector), "cr vectors");
	for (int i = 0; i < vectorCount * jointCount; ++i)
	{
		w.vectors[i] = (s2CRVector){0};
	}

	w.b = w.vectors;
	w.x = w.b + jointCount;
	w.z = w.x + jointCount;
	w.p = w.z + jointCount;
	w.az = w.p + jointCount;
	w.ap = w.az + jointCount;
	w.q = w.ap + jointCount;
	w.error = w.q + jointCount;
	w.product = w.error + jointCount;
	w.exact = w.product + jointCount;

	// Match the factored matrix: each joint's own block scaled by one plus its compliance, then the regularization
	// added to that scaled diagonal
	for (int i = 0; i < jointCount; ++i)
	{
		const s2JointRows* rows = factor->rows + i;
		float K[3][3];
		s2ComputeJointBlock(world, rows, rows, K);

		s2CRReal compliance = factor->compliance[i];
		s2CRReal B[3][3];
		for (int r = 0; r < 3; ++r)
		{
			for (int c = 0; c < 3; ++c)
			{
				bool active = r < rows->rowCount && c < rows->rowCount;
				s2CRReal extra = active ? compliance * K[r][c] : 0;
				if (active && r == c)
				{
					extra += factor->regularization * (1 + compliance) * K[r][r];
				}

				w.blocks[i].extra[r][c] = extra;
				B[r][c] = active ? K[r][c] + extra : 0;
			}
		}

		s2InvertBlock(B, rows->rowCount, w.blocks[i].inverse);
	}

	return w;
}

static void s2DestroyCRWork(s2CRWork* w)
{
	s2StackAllocator* stack = w->world->stackAllocator;
	s2FreeStackItem(stack, w->vectors);
	s2FreeStackItem(stack, w->blocks);
	s2FreeStackItem(stack, w->vx);
}

static void s2LoadVector(s2CRVector* y, const s2JointVector* x, int count)
{
	for (int i = 0; i < count; ++i)
	{
		for (int r = 0; r < 3; ++r)
		{
			y[i].row[r] = x[i].row[r];
		}
	}
}

// Preconditioned conjugate residual from Saad, Iterative Methods for Sparse Linear Systems, 2nd ed, p194.
// The residual is the preconditioned residual and the operator times the search direction follows by recurrence,
// so each iteration applies the operator once. Solves for the work's right hand side, leaving the solution in x.
static void s2RunCR(s2CRWork* w, int maxIterations, s2CRTracker* tracker)
{
	int jointCount = w->factor->jointCount;
	s2CRVector* x = w->x;
	s2CRVector* z = w->z;
	s2CRVector* p = w->p;
	s2CRVector* az = w->az;
	s2CRVector* ap = w->ap;
	s2CRVector* q = w->q;

	for (int i = 0; i < jointCount; ++i)
	{
		x[i] = (s2CRVector){0};
	}

	s2Precondition(w, w->b, z);
	s2ApplyJointOperator(w, z, az);
	for (int i = 0; i < jointCount; ++i)
	{
		p[i] = z[i];
		ap[i] = az[i];
	}

	s2CRReal rho = s2DotJoint(z, az, jointCount);
	s2CRReal rhoTolerance = (s2CRReal)S2_CR_TOLERANCE * rho;
	float relativeError = 1.0f;
	int k = 1;

	for (; k <= maxIterations; ++k)
	{
		s2Precondition(w, ap, q);

		s2CRReal denominator = s2DotJoint(ap, q, jointCount);
		if (rho <= rhoTolerance || denominator <= 0)
		{
			break;
		}

		s2CRReal alpha = rho / denominator;
		for (int i = 0; i < jointCount; ++i)
		{
			for (int r = 0; r < 3; ++r)
			{
				x[i].row[r] += alpha * p[i].row[r];
				z[i].row[r] -= alpha * q[i].row[r];
			}
		}

		if (tracker != NULL)
		{
			for (int i = 0; i < jointCount; ++i)
			{
				for (int r = 0; r < 3; ++r)
				{
					w->error[i].row[r] = x[i].row[r] - w->exact[i].row[r];
				}
			}

			s2ApplyJointOperator(w, w->error, w->product);
			s2CRReal norm2 = s2DotJoint(w->error, w->product, jointCount);
			relativeError = (float)(sqrt(norm2 > 0 ? norm2 : 0) / tracker->exactNorm);

			if (relativeError <= 0.1f && tracker->iterations1 > k)
			{
				tracker->iterations1 = k;
			}

			if (relativeError <= 0.01f && tracker->iterations2 > k)
			{
				tracker->iterations2 = k;
			}

			if (k == 8)
			{
				tracker->error8 = relativeError;
			}

			if (k == 32)
			{
				tracker->error32 = relativeError;
			}

			if (k >= 32 && tracker->iterations2 <= k)
			{
				break;
			}
		}

		if (k == maxIterations)
		{
			break;
		}

		s2ApplyJointOperator(w, z, az);
		s2CRReal rhoNew = s2DotJoint(z, az, jointCount);
		s2CRReal beta = rhoNew / rho;
		rho = rhoNew;

		for (int i = 0; i < jointCount; ++i)
		{
			for (int r = 0; r < 3; ++r)
			{
				p[i].row[r] = z[i].row[r] + beta * p[i].row[r];
				ap[i].row[r] = az[i].row[r] + beta * ap[i].row[r];
			}
		}
	}

	// Stopped before a checkpoint once the residual vanished, so it keeps the last error
	if (tracker != NULL)
	{
		if (k <= 8)
		{
			tracker->error8 = relativeError;
		}

		if (k <= 32)
		{
			tracker->error32 = relativeError;
		}
	}
}

void s2SolveJointCR(s2World* world, s2JointFactor* factor, int iterations)
{
	int jointCount = factor->jointCount;
	if (jointCount == 0)
	{
		return;
	}

	s2CRWork w = s2CreateCRWork(world, factor);
	s2LoadVector(w.b, factor->rhs, jointCount);
	s2RunCR(&w, iterations, NULL);

	for (int i = 0; i < jointCount; ++i)
	{
		for (int r = 0; r < 3; ++r)
		{
			factor->rhs[i].row[r] = (float)w.x[i].row[r];
		}
	}

	s2DestroyCRWork(&w);
}

void s2MeasureJointCR(s2World* world, const s2JointFactor* factor, const s2JointVector* rhs, s2JointCRStats* stats)
{
	int jointCount = factor->jointCount;
	if (jointCount == 0)
	{
		return;
	}

	s2CRWork w = s2CreateCRWork(world, factor);
	s2LoadVector(w.b, rhs, jointCount);
	s2LoadVector(w.exact, factor->rhs, jointCount);

	// Float floor: residual of the exact solution in the preconditioned norm
	s2ApplyJointOperator(&w, w.exact, w.product);
	s2CRReal exactNorm2 = s2DotJoint(w.exact, w.product, jointCount);
	for (int i = 0; i < jointCount; ++i)
	{
		for (int r = 0; r < 3; ++r)
		{
			w.error[i].row[r] = w.b[i].row[r] - w.product[i].row[r];
		}
	}

	s2Precondition(&w, w.error, w.q);
	s2Precondition(&w, w.b, w.z);
	s2CRReal rhsNorm2 = s2DotJoint(w.b, w.z, jointCount);
	s2CRReal residual2 = s2DotJoint(w.error, w.q, jointCount);
	float exactResidual = rhsNorm2 > 0 ? (float)sqrt((residual2 > 0 ? residual2 : 0) / rhsNorm2) : 0.0f;

	s2CRTracker tracker = {0};
	if (exactNorm2 > 0)
	{
		tracker.exactNorm = (s2CRReal)sqrt(exactNorm2);
		tracker.iterations1 = S2_CR_MAX_MEASURE_ITERATIONS + 1;
		tracker.iterations2 = S2_CR_MAX_MEASURE_ITERATIONS + 1;
		tracker.error8 = 1.0f;
		tracker.error32 = 1.0f;
		s2RunCR(&w, S2_CR_MAX_MEASURE_ITERATIONS, &tracker);
	}

	stats->solveCount += 1;
	stats->iterations1 = S2_MAX(stats->iterations1, tracker.iterations1);
	stats->iterations2 = S2_MAX(stats->iterations2, tracker.iterations2);
	stats->error8 = S2_MAX(stats->error8, tracker.error8);
	stats->error32 = S2_MAX(stats->error32, tracker.error32);
	stats->exactResidual = S2_MAX(stats->exactResidual, exactResidual);

	s2DestroyCRWork(&w);
}
