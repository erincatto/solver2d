// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "joint.h"
#include "joint_factor.h"

#include <stdbool.h>

// Soft Step's per joint update dlambda = -m K^-1 (Cdot + b) - i lambda becomes, over all joints,
// (D + G) dlambda = -(J v + b) - G lambda with D = J M^-1 J^T and compliance G = K / c per joint.
// These pieces are shared by the solvers that solve the equality rows of revolute, weld, and prismatic joints
// together with the joint factor.

// What a direct solve needs from a revolute, weld, or prismatic joint
typedef struct s2DirectJoint
{
	float biasCoefficient;
	float massCoefficient;
	float impulseCoefficient;

	// Accumulated impulse of each equality row
	float* impulses[3];

	// Revolute motor and limits or prismatic spring and limits, warm started with the direct rows
	float axialImpulse;
} s2DirectJoint;

s2DirectJoint s2GetDirectJoint(s2Joint* base);

// Position error of each equality row at the current body state, even when the Jacobian anchors are fixed.
// Measured as each joint's own soft solve does.
void s2GetDirectError(s2World* world, s2Joint* base, float C[3]);

// Creates the factor over the revolute, weld, and prismatic joints that touch a dynamic body and orders it.
// Without all, only joints that opt in.
void s2CreateDirectFactor(s2World* world, s2JointFactor* factor, bool all);

void s2UpdateDirectRows(s2World* world, s2JointFactor* factor);

// Applies the stored impulses of the direct joints, equality rows with motor, spring, and limits, times the scale
void s2WarmStartDirect(s2World* world, s2JointFactor* factor, float scale);

// Right hand side of the equality rows from the current body velocities. With bias it adds the position bias and
// Soft Step's impulse leak on the stored impulses times the scale. The bias rate comes from the softness, or from
// each joint's own when it is null. The leak uses the factor's compliance.
void s2BuildDirectRhs(s2World* world, s2JointFactor* factor, bool useBias, const s2Softness* softness, float scale);

// Adds the solution left in the factor to the stored impulses, divided by the scale. Optionally applies it to the
// bodies.
void s2ApplyDirectSolution(s2World* world, s2JointFactor* factor, float scale, bool applyToBodies);
