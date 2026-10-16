// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "joint_factor.h"

// Matrix free conjugate residual on the system the joint factor solves, (D + G) x = b with the factor's rows,
// compliance, and regularization. Preconditioned by each joint's own block. Starts from zero.

// Solves the right hand side in the factor in place with a fixed iteration count. Needs no factorization.
void s2SolveJointCR(s2World* world, s2JointFactor* factor, int iterations);

// How conjugate residual does against the exact solve. Worst case over the solves since the last reset.
typedef struct s2JointCRStats
{
	int solveCount;

	// Iterations until the energy norm error relative to the exact solution falls below 0.1 and 0.01.
	// One past the iteration limit when it never does.
	int iterations1;
	int iterations2;

	// Relative energy norm error after 8 and 32 iterations
	float error8;
	float error32;

	// Relative residual of the exact solution, the float floor under the error
	float exactResidual;
} s2JointCRStats;

// Runs conjugate residual on the system the factor just solved. The right hand side is the one the solve consumed.
// The exact solution is what the solve left in the factor. Bodies are not touched.
void s2MeasureJointCR(s2World* world, const s2JointFactor* factor, const s2JointVector* rhs, s2JointCRStats* stats);
