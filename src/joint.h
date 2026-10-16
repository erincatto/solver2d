// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver2d/id.h"
#include "solver2d/types.h"

#include "pool.h"

#include <stdint.h>

typedef struct s2DebugDraw s2DebugDraw;
typedef struct s2StepContext s2StepContext;
typedef struct s2World s2World;

typedef enum s2JointType
{
	s2_revoluteJoint,
	s2_mouseJoint,
	s2_weldJoint,
	s2_prismaticJoint,
} s2JointType;

typedef struct s2JointEdge
{
	int32_t bodyIndex;
	int32_t prevKey;
	int32_t nextKey;
} s2JointEdge;

typedef struct s2MouseJoint
{
	float hertz;
	float dampingRatio;

	s2Vec2 targetA;
	
	// Solver shared
	s2Vec2 impulse;
	float motorImpulse;

	float biasCoefficient;
	float massCoefficient;
	float impulseCoefficient;

	s2Vec2 localAnchorB;
	float invMassB;
	float invIB;
	s2Mat22 pivotMass;
	s2Vec2 centerDiff0;
} s2MouseJoint;

typedef struct s2RevoluteJoint
{
	// Solver shared
	s2Vec2 impulse;
	float motorImpulse;
	float lowerImpulse;
	float upperImpulse;
	bool enableMotor;
	float maxMotorTorque;
	float motorSpeed;
	bool enableLimit;
	bool enableDirect;
	float referenceAngle;
	float lowerAngle;
	float upperAngle;

	// Solver temp

	// Local anchors relative to center of mass
	s2Vec2 localAnchorA;
	s2Vec2 localAnchorB;
	s2Vec2 centerDiff0;
	float invMassA;
	float invMassB;
	float invIA;
	float invIB;
	s2Mat22 pivotMass;
	float biasCoefficient;
	float massCoefficient;
	float impulseCoefficient;
	float axialMass;

	float inertiaScale;
} s2RevoluteJoint;

typedef struct s2WeldJoint
{
	// Solver shared
	s2Vec2 linearImpulse;
	float angularImpulse;
	float referenceAngle;
	bool enableDirect;

	// Solver temp
	s2Vec2 localAnchorA;
	s2Vec2 localAnchorB;
	s2Vec2 centerDiff0;
	float invMassA;
	float invMassB;
	float invIA;
	float invIB;
	float axialMass;
	float biasCoefficient;
	float massCoefficient;
	float impulseCoefficient;
} s2WeldJoint;

typedef struct s2Softness
{
	float biasRate;
	float massScale;
	float impulseScale;
} s2Softness;

typedef struct s2PrismaticJoint
{
	// Solver shared
	// Perpendicular and angle impulses
	s2Vec2 impulse;
	float springImpulse;
	float lowerImpulse;
	float upperImpulse;
	s2Vec2 localAxisA;
	float referenceAngle;
	bool enableSpring;
	float hertz;
	float dampingRatio;
	float targetTranslation;
	bool enableLimit;
	bool enableDirect;
	float lowerTranslation;
	float upperTranslation;

	// Solver temp
	s2Vec2 localAnchorA;
	s2Vec2 localAnchorB;
	s2Vec2 centerDiff0;
	float invMassA;
	float invMassB;
	float invIA;
	float invIB;
	s2Softness softness;
	s2Softness limitSoftness;
	s2Softness springSoftness;
} s2PrismaticJoint;

/// The base joint class. Joints are used to constraint two bodies together in
/// various fashions. Some joints also feature limits and motors.
typedef struct s2Joint
{
	s2Object object;
	s2JointType type;
	s2JointEdge edges[2];
	s2Vec2 localOriginAnchorA;
	s2Vec2 localOriginAnchorB;
	float drawSize;

	union
	{
		s2MouseJoint mouseJoint;
		s2RevoluteJoint revoluteJoint;
		s2WeldJoint weldJoint;
		s2PrismaticJoint prismaticJoint;
	};

	bool collideConnected;
} s2Joint;

// From Box2D. Zero hertz gives no force.
s2Softness s2MakeSoft(float hertz, float zeta, float h);

// shared
void s2WarmStartJoint(s2Joint* joint, s2StepContext* context);
void s2SolveJointPosition(s2Joint* joint, s2StepContext* context);

void s2PrepareJoint(s2Joint* joint, s2StepContext* context, bool warmStart);
void s2SolveJoint(s2Joint* joint, s2StepContext* context, float h);

void s2PrepareJoint_Soft(s2Joint* joint, s2StepContext* context, float h, float hertz, bool warmStart);
void s2SolveJoint_Soft(s2Joint* joint, s2StepContext* context, float h, float inv_h, bool useBias);
void s2SolveRevoluteAxial_Soft(s2Joint* base, s2StepContext* context, float h, float inv_h, bool useBias);
void s2SolvePrismaticAxial_Soft(s2Joint* base, s2StepContext* context, float inv_h, bool useBias);
void s2PreparePrismaticLimit_Soft(s2Joint* base, float h, float hertz);

void s2SolveJoint_Baumgarte(s2Joint* joint, s2StepContext* context, float h, float inv_h, bool useBias);

void s2PrepareJoint_XPBD(s2Joint* joint, s2StepContext* context);
void s2SolveJoint_XPBD(s2Joint* joint, s2StepContext* context, float inv_h);
void s2SolveJointVelocity_XPBD(s2Joint* joint, s2StepContext* context, float h);

void s2DrawJoint(s2DebugDraw* draw, s2World* world, s2Joint* joint);
