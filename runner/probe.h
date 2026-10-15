// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver2d/id.h"
#include "solver2d/types.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct ProbeStep
{
	int bodyCount;
	int contactCount;
	int pointCount;
	float maxPenetration;
	float meanPenetration;
	float normalImpulse;
	int jointCount;
	float maxJointError;
	float meanJointError;
	float maxLimitError;
	float kineticEnergy;
	float potentialEnergy;
	float maxSpeed;
	float maxAngularSpeed;
} ProbeStep;

typedef struct ProbeBody
{
	int index;
	s2Vec2 position;
	float angle;
	s2Vec2 linearVelocity;
	float angularVelocity;
} ProbeBody;

typedef struct Probe Probe;

Probe* CreateProbe(void);
void DestroyProbe(Probe* probe);

/// Record body poses. Call right before stepping so contact separation can be measured after the step.
void Probe_Snapshot(Probe* probe, s2WorldId worldId);

/// Measure the world state left by the step that followed the last snapshot.
ProbeStep Probe_Measure(const Probe* probe, s2WorldId worldId);

/// Gather dynamic bodies.
/// @return the dynamic body count, which may exceed capacity
int Probe_GetBodies(s2WorldId worldId, ProbeBody* bodies, int capacity);

#ifdef __cplusplus
}
#endif
