// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "broad_phase.h"
#include "pool.h"

#include "solver2d/callbacks.h"
#include "solver2d/timer.h"

#define B2_GRAPH_COLOR 1

typedef struct s2Contact s2Contact;

typedef struct s2World
{
	int16_t index;

	s2SolverType solverType;

	struct s2BlockAllocator* blockAllocator;
	struct s2StackAllocator* stackAllocator;

	s2BroadPhase broadPhase;

	s2Pool bodyPool;
	s2Pool contactPool;
	s2Pool jointPool;
	s2Pool shapePool;
	s2Pool islandPool;

	// These are sparse arrays that point into the pools above
	struct s2Body* bodies;
	struct s2Contact* contacts;
	struct s2Joint* joints;
	struct s2Shape* shapes;

	// Id that is incremented every time step
	uint64_t stepId;

	s2Vec2 gravity;
	float restitutionThreshold;
	uint16_t revision;
} s2World;

s2World* s2GetWorldFromId(s2WorldId id);
s2World* s2GetWorldFromIndex(int16_t index);

bool s2IsBodyIdValid(s2World* world, s2BodyId id);
