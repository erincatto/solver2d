// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver2d/math.h"

typedef struct s2World s2World;

// Equality rows of a revolute or weld joint. Rows 0 and 1 are the point constraint, row 2 is the weld angle.
// Cdot = vB + cross(wB, rB) - vA - cross(wA, rA)
// Cdot = wB - wA
typedef struct s2JointRows
{
	int indexA;
	int indexB;
	s2Vec2 rA;
	s2Vec2 rB;
	int rowCount;
} s2JointRows;

// One value per row, the angular value only for welds
typedef struct s2JointVector
{
	s2Vec2 linear;
	float angular;
} s2JointVector;

// Direct solver for the coupled equality rows of many joints, J M^-1 J^T.
// Envelope Cholesky in reverse Cuthill-McKee order. Joints are adjacent when they share a dynamic body.
// Per joint arrays are indexed in the order the caller lists the joints.
typedef struct s2JointFactor
{
	int jointCount;
	int rowCount;

	// Filled by the caller. The joint indices and rows before ordering, the anchors and compliance before
	// factoring, the right hand side before solving. The solve overwrites the right hand side with the solution.
	int* jointIndices;
	s2JointRows* rows;
	float* compliance;
	s2JointVector* rhs;

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
