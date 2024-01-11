// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include <stdint.h>

typedef struct s2Body s2Body;
typedef struct s2Contact s2Contact;
typedef struct s2StepContext s2StepContext;
typedef struct s2World s2World;

typedef struct s2StepContext
{
	float dt;
	float inv_dt;
	int32_t velocityIterations;
	int32_t positionIterations;
	float restitutionThreshold;
	s2Body* bodies;
	int32_t bodyCapacity;
} s2StepContext;

typedef struct s2ContactConstraintPoint
{
	s2Vec2 rA, rB;
	s2Vec2 rAf, rBf;
	s2Vec2 localAnchorA, localAnchorB;
	float tangentSeparation;
	float separation;
	float normalImpulse;
	float tangentImpulse;
	float normalMass;
	float tangentMass;
	float gamma;
	float massCoefficient;
	float biasCoefficient;
	float impulseCoefficient;
	float baumgarte;
	bool frictionValid;
} s2ContactConstraintPoint;

typedef struct s2ContactConstraint
{
	s2Contact* contact;
	int indexA;
	int indexB;
	s2ContactConstraintPoint points[2];
	s2Vec2 normal;
	float friction;
	float kinematicFriction;
	float restitution;
	int pointCount;
} s2ContactConstraint;

// common
void s2IntegrateVelocities(s2World* world, float h);
void s2IntegratePositions(s2World* world, float h);
void s2WarmStartContacts(s2World* world, s2ContactConstraint* constraints, int constraintCount);
void s2StoreContactImpulses(s2ContactConstraint* constraints, int constraintCount);

void s2Solve_PGS_NGS_Block(s2World* world, s2StepContext* stepContext);
void s2Solve_PGS_NGS(s2World* world, s2StepContext* context);
void s2Solve_PGS_Soft(s2World* world, s2StepContext* stepContext);
void s2Solve_XPDB(s2World* world, s2StepContext* stepContext);
void s2Solve_TGS_Soft(s2World* world, s2StepContext* stepContext);
void s2Solve_TGS_Sticky(s2World* world, s2StepContext* stepContext);
