// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver2d/math.h"

typedef struct s2Joint s2Joint;
typedef struct s2World s2World;

// Equality rows of a revolute, weld, or prismatic joint. Point rows come first, then an optional angle row.
// Revolute: the x and y point rows. Weld: the point rows and the angle row. Prismatic: one point row along the
// perpendicular of the axis and the angle row.
// Cdot = dot(axis, vB + cross(wB, rB) - vA - cross(wA, rA))
// Cdot = wB - wA
typedef struct s2JointRows
{
	int indexA;
	int indexB;
	s2Vec2 rA;
	s2Vec2 rB;
	s2Vec2 axes[2];
	int pointCount;
	int rowCount;
} s2JointRows;

typedef struct s2JointVector
{
	float row[3];
} s2JointVector;

// Rows at the current body state. The prismatic point row acts on body A at the point coincident with the anchor
// on B.
void s2GetJointRows(s2World* world, s2Joint* joint, s2JointRows* rows);

// Where the joint accumulates the impulse of each row
void s2GetJointRowImpulses(s2Joint* joint, float* impulses[3]);

static inline s2Vec2 s2GetPointImpulse(const s2JointRows* rows, const float* lambda)
{
	s2Vec2 P = s2Vec2_zero;
	for (int r = 0; r < rows->pointCount; ++r)
	{
		P = s2MulAdd(P, lambda[r], rows->axes[r]);
	}
	return P;
}

static inline float s2GetAngleImpulse(const s2JointRows* rows, const float* lambda)
{
	return rows->rowCount > rows->pointCount ? lambda[rows->pointCount] : 0.0f;
}

// Direct solver for the coupled equality rows of many joints, J M^-1 J^T.
// Envelope Cholesky in reverse Cuthill-McKee order. Joints are adjacent when they share a dynamic body.
// Per joint arrays are indexed in the order the caller lists the joints.
typedef struct s2JointFactor
{
	int jointCount;
	int rowCount;

	// Filled by the caller. The joint indices and rows before ordering, the anchors and axes and compliance before
	// factoring, the right hand side before solving. The solve overwrites the right hand side with the solution.
	int* jointIndices;
	s2JointRows* rows;
	float* compliance;
	s2JointVector* rhs;

	// Set by the factorization
	float regularization;

	// Joint order[p] owns the rows starting at rowStart[p]
	int* order;
	int* position;
	int* rowStart;

	int* jointMap;
	int* adjacencyStart;
	int* adjacency;

	// Cholesky factor stored by rows from the first nonzero column to the diagonal
	int* first;
	int* offset;
	float* values;
	float* scale;
	float* work;
} s2JointFactor;

void s2CreateJointFactor(s2World* world, s2JointFactor* factor, int jointCount);

// Builds the joint graph and its ordering. Needs the joint indices and the row counts.
void s2OrderJointFactor(s2World* world, s2JointFactor* factor);

// Each joint's own block is scaled by one plus its compliance. Regularization is added to the unit diagonal
// after symmetric scaling and is also the pivot floor.
void s2FactorJoints(s2World* world, s2JointFactor* factor, float regularization);

void s2SolveJointFactor(s2JointFactor* factor);

void s2DestroyJointFactor(s2World* world, s2JointFactor* factor);

// Delassus block coupling the rows of joints p and q through their shared bodies
void s2ComputeJointBlock(s2World* world, const s2JointRows* p, const s2JointRows* q, float K[3][3]);
