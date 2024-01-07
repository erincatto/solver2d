// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver_data.h"
#include "stack_allocator.h"

#include "solver2d/callbacks.h"

typedef struct s2ContactSolverDef
{
	const s2StepContext* context;
	struct s2World* world;
	int32_t contactList;
	int32_t contactCount;
} s2ContactSolverDef;

typedef struct s2ContactSolver
{
	const s2StepContext* context;
	struct s2World* world;
	struct s2ContactPositionConstraint* positionConstraints;
	struct s2ContactVelocityConstraint* velocityConstraints;
	int32_t contactList;
	int32_t contactCount;
	int32_t constraintCount;
} s2ContactSolver;

s2ContactSolver* s2CreateContactSolver(s2ContactSolverDef* def);

static inline void s2DestroyContactSolver(s2ContactSolver* solver, s2StackAllocator* alloc)
{
	s2FreeStackItem(alloc, solver->velocityConstraints);
	s2FreeStackItem(alloc, solver->positionConstraints);
	s2FreeStackItem(alloc, solver);
}

void s2ContactSolver_Initialize(s2ContactSolver* solver);
void s2ContactSolver_SolveVelocityConstraints(s2ContactSolver* solver);
void s2ContactSolver_ApplyRestitution(s2ContactSolver* solver);
void s2ContactSolver_StoreImpulses(s2ContactSolver* solver);
bool s2ContactSolver_SolvePositionConstraintsBlock(s2ContactSolver* solver);
