// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver2d/math.h"

typedef struct s2World s2World;

// Point rows of a revolute joint
// Cdot = vB + cross(wB, rB) - vA - cross(wA, rA)
typedef struct s2PointRows
{
	int indexA;
	int indexB;
	s2Vec2 rA;
	s2Vec2 rB;
} s2PointRows;

// Direct solver for the coupled point rows of many joints, J M^-1 J^T with two rows per joint.
// Envelope Cholesky in reverse Cuthill-McKee order. Joints are adjacent when they share a dynamic body.
// Per joint arrays are indexed in the order the caller lists the joints.
typedef struct s2JointFactor
{
	int jointCount;

	// Filled by the caller. The joint indices before ordering, the rows and compliance before factoring,
	// the right hand side before solving. The solve overwrites the right hand side with the solution.
	int* jointIndices;
	s2PointRows* rows;
	float* compliance;
	s2Vec2* rhs;

	// Row pair p belongs to joint order[p]
	int* order;
	int* position;

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

// Builds the joint graph and its ordering. Needs the joint indices.
void s2OrderJointFactor(s2World* world, s2JointFactor* factor);

// Each joint's own block is scaled by one plus its compliance. Regularization is added to the unit diagonal
// after symmetric scaling and is also the pivot floor.
void s2FactorJoints(s2World* world, s2JointFactor* factor, float regularization);

void s2SolveJointFactor(s2JointFactor* factor);

void s2DestroyJointFactor(s2World* world, s2JointFactor* factor);

// Delassus block coupling the point rows of joints p and q through their shared bodies
void s2ComputeJointBlock(s2World* world, const s2PointRows* p, const s2PointRows* q, float K[2][2]);
