// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver_data.h"
#include "stack_allocator.h"

#include "solver2d/callbacks.h"

typedef struct s2World s2World;
typedef struct s2StepContext s2StepContext;

typedef struct s2ContactSolver
{
	s2World* world;
	s2StepContext* context;
	struct s2ContactVelocityConstraint* velocityConstraints;
	struct s2ContactPositionConstraint* positionConstraints;
	int constraintCount;
} s2ContactSolver;

s2ContactSolver s2CreateContactSolver(s2World* world, s2StepContext* context);
void s2DestroyContactSolver(s2ContactSolver* solver, s2StackAllocator* alloc);

void s2ContactSolver_SolveVelocityConstraints(s2ContactSolver* solver);
void s2ContactSolver_ApplyRestitution(s2ContactSolver* solver);
void s2ContactSolver_StoreImpulses(s2ContactSolver* solver);
bool s2ContactSolver_SolvePositionConstraintsBlock(s2ContactSolver* solver);
