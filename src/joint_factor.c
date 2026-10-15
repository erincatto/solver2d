// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "joint_factor.h"

#include "body.h"
#include "core.h"
#include "joint.h"
#include "stack_allocator.h"
#include "world.h"

#include <float.h>
#include <math.h>

void s2CreateJointFactor(s2World* world, s2JointFactor* factor, int jointCount)
{
	s2StackAllocator* stack = world->stackAllocator;
	int count = S2_MAX(1, jointCount);

	*factor = (s2JointFactor){0};
	factor->jointCount = jointCount;
	factor->jointIndices = s2AllocateStackItem(stack, count * sizeof(int), "joint indices");
	factor->rows = s2AllocateStackItem(stack, count * sizeof(s2JointRows), "joint rows");
	factor->compliance = s2AllocateStackItem(stack, count * sizeof(float), "joint compliance");
	factor->rhs = s2AllocateStackItem(stack, count * sizeof(s2JointVector), "joint rhs");

	for (int i = 0; i < jointCount; ++i)
	{
		factor->compliance[i] = 0.0f;
	}
}

// Reverse Cuthill-McKee ordering of the joint graph for a small envelope.
// From George and Liu, Computer Solution of Large Sparse Positive Definite Systems, 1981.
static void s2OrderJoints(int count, const int* adjacencyStart, const int* adjacency, int* order, int* position)
{
	for (int i = 0; i < count; ++i)
	{
		position[i] = -1;
	}

	int orderedCount = 0;
	while (orderedCount < count)
	{
		// Start each connected component from a node of minimum degree
		int start = -1;
		int startDegree = INT32_MAX;
		for (int i = 0; i < count; ++i)
		{
			int degree = adjacencyStart[i + 1] - adjacencyStart[i];
			if (position[i] == -1 && degree < startDegree)
			{
				start = i;
				startDegree = degree;
			}
		}

		position[start] = orderedCount;
		order[orderedCount++] = start;

		int head = orderedCount - 1;
		while (head < orderedCount)
		{
			int node = order[head++];
			int levelStart = orderedCount;

			for (int k = adjacencyStart[node]; k < adjacencyStart[node + 1]; ++k)
			{
				int other = adjacency[k];
				if (position[other] == -1)
				{
					position[other] = orderedCount;
					order[orderedCount++] = other;
				}
			}

			// Visit new neighbors by increasing degree
			for (int k = levelStart + 1; k < orderedCount; ++k)
			{
				int node1 = order[k];
				int degree1 = adjacencyStart[node1 + 1] - adjacencyStart[node1];
				int m = k;
				while (m > levelStart && adjacencyStart[order[m - 1] + 1] - adjacencyStart[order[m - 1]] > degree1)
				{
					order[m] = order[m - 1];
					m -= 1;
				}
				order[m] = node1;
			}
		}
	}

	for (int i = 0; i < count / 2; ++i)
	{
		int temp = order[i];
		order[i] = order[count - 1 - i];
		order[count - 1 - i] = temp;
	}

	for (int i = 0; i < count; ++i)
	{
		position[order[i]] = i;
	}
}

void s2OrderJointFactor(s2World* world, s2JointFactor* factor)
{
	s2StackAllocator* stack = world->stackAllocator;
	s2Body* bodies = world->bodies;
	s2Joint* joints = world->joints;
	int jointCapacity = world->jointPool.capacity;
	int jointCount = factor->jointCount;
	int count = S2_MAX(1, jointCount);

	factor->order = s2AllocateStackItem(stack, count * sizeof(int), "joint order");
	factor->position = s2AllocateStackItem(stack, count * sizeof(int), "joint position");
	factor->rowStart = s2AllocateStackItem(stack, (jointCount + 1) * sizeof(int), "joint row start");

	factor->jointMap = s2AllocateStackItem(stack, jointCapacity * sizeof(int), "joint map");
	for (int i = 0; i < jointCapacity; ++i)
	{
		factor->jointMap[i] = S2_NULL_INDEX;
	}

	for (int i = 0; i < jointCount; ++i)
	{
		factor->jointMap[factor->jointIndices[i]] = i;
	}

	int* adjacencyStart = s2AllocateStackItem(stack, (jointCount + 1) * sizeof(int), "joint adjacency start");
	int adjacencyCapacity = 0;
	for (int i = 0; i < jointCount; ++i)
	{
		adjacencyStart[i] = adjacencyCapacity;
		s2Body* bodyA = bodies + factor->rows[i].indexA;
		s2Body* bodyB = bodies + factor->rows[i].indexB;
		adjacencyCapacity += bodyA->type == s2_dynamicBody ? bodyA->jointCount : 0;
		adjacencyCapacity += bodyB->type == s2_dynamicBody ? bodyB->jointCount : 0;
	}
	adjacencyStart[jointCount] = adjacencyCapacity;

	int* adjacency = s2AllocateStackItem(stack, S2_MAX(1, adjacencyCapacity) * sizeof(int), "joint adjacency");
	int adjacencyCount = 0;
	for (int i = 0; i < jointCount; ++i)
	{
		int capacity = adjacencyStart[i + 1];
		adjacencyStart[i] = adjacencyCount;

		int indices[2] = {factor->rows[i].indexA, factor->rows[i].indexB};
		for (int side = 0; side < 2; ++side)
		{
			s2Body* body = bodies + indices[side];
			if (body->type != s2_dynamicBody)
			{
				continue;
			}

			int jointKey = body->jointList;
			while (jointKey != S2_NULL_INDEX)
			{
				int jointIndex = jointKey >> 1;
				int other = factor->jointMap[jointIndex];
				if (other != S2_NULL_INDEX && other != i)
				{
					S2_ASSERT(adjacencyCount < capacity);
					adjacency[adjacencyCount++] = other;
				}
				jointKey = joints[jointIndex].edges[jointKey & 1].nextKey;
			}
		}
	}
	adjacencyStart[jointCount] = adjacencyCount;

	factor->adjacencyStart = adjacencyStart;
	factor->adjacency = adjacency;

	if (jointCount > 0)
	{
		s2OrderJoints(jointCount, adjacencyStart, adjacency, factor->order, factor->position);
	}

	int rowCount = 0;
	for (int p = 0; p < jointCount; ++p)
	{
		factor->rowStart[p] = rowCount;
		rowCount += factor->rows[factor->order[p]].rowCount;
	}
	factor->rowStart[jointCount] = rowCount;
	factor->rowCount = rowCount;

	// Fill stays inside the envelope, so it is fixed by the ordering
	factor->first = s2AllocateStackItem(stack, S2_MAX(1, rowCount) * sizeof(int), "joint first");
	factor->offset = s2AllocateStackItem(stack, (rowCount + 1) * sizeof(int), "joint offset");

	int size = 0;
	for (int p = 0; p < jointCount; ++p)
	{
		int node = factor->order[p];
		int firstBlock = p;
		for (int k = adjacencyStart[node]; k < adjacencyStart[node + 1]; ++k)
		{
			firstBlock = S2_MIN(firstBlock, factor->position[adjacency[k]]);
		}

		int firstRow = factor->rowStart[firstBlock];
		for (int row = factor->rowStart[p]; row < factor->rowStart[p + 1]; ++row)
		{
			factor->first[row] = firstRow;
			factor->offset[row] = size;
			size += row - firstRow + 1;
		}
	}
	factor->offset[rowCount] = size;

	factor->values = s2AllocateStackItem(stack, S2_MAX(1, size) * sizeof(float), "joint envelope");
	factor->scale = s2AllocateStackItem(stack, S2_MAX(1, rowCount) * sizeof(float), "joint scale");
	factor->work = s2AllocateStackItem(stack, S2_MAX(1, rowCount) * sizeof(float), "joint work");
}

void s2DestroyJointFactor(s2World* world, s2JointFactor* factor)
{
	s2StackAllocator* stack = world->stackAllocator;

	if (factor->order != NULL)
	{
		s2FreeStackItem(stack, factor->work);
		s2FreeStackItem(stack, factor->scale);
		s2FreeStackItem(stack, factor->values);
		s2FreeStackItem(stack, factor->offset);
		s2FreeStackItem(stack, factor->first);
		s2FreeStackItem(stack, factor->adjacency);
		s2FreeStackItem(stack, factor->adjacencyStart);
		s2FreeStackItem(stack, factor->jointMap);
		s2FreeStackItem(stack, factor->rowStart);
		s2FreeStackItem(stack, factor->position);
		s2FreeStackItem(stack, factor->order);
	}

	s2FreeStackItem(stack, factor->rhs);
	s2FreeStackItem(stack, factor->compliance);
	s2FreeStackItem(stack, factor->rows);
	s2FreeStackItem(stack, factor->jointIndices);
	*factor = (s2JointFactor){0};
}

// Jacobian of one row on one body. The sign is -1 for body A and 1 for body B.
static void s2GetRowJacobian(s2Vec2 r, float sign, int row, s2Vec2* linear, float* angular)
{
	if (row == 2)
	{
		*linear = s2Vec2_zero;
		*angular = sign;
		return;
	}

	// Angular Jacobian of the x and y rows is the perp of the anchor
	s2Vec2 perp = s2LeftPerp(r);
	*linear = row == 0 ? (s2Vec2){sign, 0.0f} : (s2Vec2){0.0f, sign};
	*angular = sign * (row == 0 ? perp.x : perp.y);
}

void s2ComputeJointBlock(s2World* world, const s2JointRows* p, const s2JointRows* q, float K[3][3])
{
	for (int r = 0; r < 3; ++r)
	{
		for (int c = 0; c < 3; ++c)
		{
			K[r][c] = 0.0f;
		}
	}

	int indicesP[2] = {p->indexA, p->indexB};
	int indicesQ[2] = {q->indexA, q->indexB};
	s2Vec2 anchorsP[2] = {p->rA, p->rB};
	s2Vec2 anchorsQ[2] = {q->rA, q->rB};
	float signs[2] = {-1.0f, 1.0f};

	for (int a = 0; a < 2; ++a)
	{
		for (int b = 0; b < 2; ++b)
		{
			if (indicesP[a] != indicesQ[b])
			{
				continue;
			}

			s2Body* body = world->bodies + indicesP[a];
			float m = body->invMass;
			float i = body->invI;

			for (int r = 0; r < p->rowCount; ++r)
			{
				s2Vec2 linearP;
				float angularP;
				s2GetRowJacobian(anchorsP[a], signs[a], r, &linearP, &angularP);

				for (int c = 0; c < q->rowCount; ++c)
				{
					s2Vec2 linearQ;
					float angularQ;
					s2GetRowJacobian(anchorsQ[b], signs[b], c, &linearQ, &angularQ);

					K[r][c] += m * s2Dot(linearP, linearQ) + i * angularP * angularQ;
				}
			}
		}
	}
}

static inline float* s2EnvelopeEntry(s2JointFactor* factor, int row, int column)
{
	S2_ASSERT(factor->first[row] <= column && column <= row);
	return factor->values + factor->offset[row] + column - factor->first[row];
}

void s2FactorJoints(s2World* world, s2JointFactor* factor, float regularization)
{
	int jointCount = factor->jointCount;
	int rowCount = factor->rowCount;
	const int* order = factor->order;
	const int* position = factor->position;
	const int* adjacencyStart = factor->adjacencyStart;
	const int* adjacency = factor->adjacency;

	int size = factor->offset[rowCount];
	for (int i = 0; i < size; ++i)
	{
		factor->values[i] = 0.0f;
	}

	for (int p = 0; p < jointCount; ++p)
	{
		int node = order[p];
		const s2JointRows* rowsP = factor->rows + node;
		int startP = factor->rowStart[p];

		float K[3][3];
		s2ComputeJointBlock(world, rowsP, rowsP, K);
		float diagonalScale = 1.0f + factor->compliance[node];
		for (int r = 0; r < rowsP->rowCount; ++r)
		{
			for (int c = 0; c <= r; ++c)
			{
				*s2EnvelopeEntry(factor, startP + r, startP + c) = diagonalScale * K[r][c];
			}
		}

		for (int k = adjacencyStart[node]; k < adjacencyStart[node + 1]; ++k)
		{
			int q = position[adjacency[k]];
			if (q >= p)
			{
				continue;
			}

			const s2JointRows* rowsQ = factor->rows + adjacency[k];
			int startQ = factor->rowStart[q];
			s2ComputeJointBlock(world, rowsP, rowsQ, K);
			for (int r = 0; r < rowsP->rowCount; ++r)
			{
				for (int c = 0; c < rowsQ->rowCount; ++c)
				{
					*s2EnvelopeEntry(factor, startP + r, startQ + c) = K[r][c];
				}
			}
		}
	}

	// Symmetric scaling to a unit diagonal preserves float precision across mass ratios
	for (int i = 0; i < rowCount; ++i)
	{
		factor->scale[i] = 1.0f / sqrtf(S2_ABS(*s2EnvelopeEntry(factor, i, i)) + FLT_EPSILON);
	}

	for (int i = 0; i < rowCount; ++i)
	{
		for (int j = factor->first[i]; j <= i; ++j)
		{
			*s2EnvelopeEntry(factor, i, j) *= factor->scale[i] * factor->scale[j];
		}
		*s2EnvelopeEntry(factor, i, i) += regularization;
	}

	// Row oriented Cholesky. Entries left of the envelope are zero in both the matrix and the factor.
	float pivotFloor = S2_MAX(regularization, FLT_EPSILON);
	for (int i = 0; i < rowCount; ++i)
	{
		int firstI = factor->first[i];
		float* rowI = factor->values + factor->offset[i] - firstI;

		for (int j = firstI; j <= i; ++j)
		{
			int firstJ = factor->first[j];
			float* rowJ = factor->values + factor->offset[j] - firstJ;

			float sum = rowI[j];
			for (int k = S2_MAX(firstI, firstJ); k < j; ++k)
			{
				sum -= rowI[k] * rowJ[k];
			}

			if (j < i)
			{
				rowI[j] = sum / rowJ[j];
			}
			else
			{
				// Roundoff on redundant joints can leave a nonpositive pivot
				rowI[i] = sqrtf(S2_MAX(sum, pivotFloor));
			}
		}
	}
}

void s2SolveJointFactor(s2JointFactor* factor)
{
	int jointCount = factor->jointCount;
	int rowCount = factor->rowCount;
	float* x = factor->work;
	const float* scale = factor->scale;

	for (int p = 0; p < jointCount; ++p)
	{
		int node = factor->order[p];
		int start = factor->rowStart[p];
		s2JointVector b = factor->rhs[node];
		x[start] = scale[start] * b.linear.x;
		x[start + 1] = scale[start + 1] * b.linear.y;
		if (factor->rows[node].rowCount == 3)
		{
			x[start + 2] = scale[start + 2] * b.angular;
		}
	}

	for (int i = 0; i < rowCount; ++i)
	{
		int firstI = factor->first[i];
		const float* rowI = factor->values + factor->offset[i] - firstI;
		float sum = x[i];
		for (int k = firstI; k < i; ++k)
		{
			sum -= rowI[k] * x[k];
		}
		x[i] = sum / rowI[i];
	}

	for (int i = rowCount - 1; i >= 0; --i)
	{
		int firstI = factor->first[i];
		const float* rowI = factor->values + factor->offset[i] - firstI;
		x[i] /= rowI[i];
		for (int k = firstI; k < i; ++k)
		{
			x[k] -= rowI[k] * x[i];
		}
	}

	for (int p = 0; p < jointCount; ++p)
	{
		int node = factor->order[p];
		int start = factor->rowStart[p];
		s2JointVector* solution = factor->rhs + node;
		solution->linear = (s2Vec2){scale[start] * x[start], scale[start + 1] * x[start + 1]};
		solution->angular = factor->rows[node].rowCount == 3 ? scale[start + 2] * x[start + 2] : 0.0f;
	}
}
