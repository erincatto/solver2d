// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver2d/types.h"


typedef struct s2MouseJointDef
{
	/// The first attached body.
	s2BodyId bodyIdA;

	/// The second attached body.
	s2BodyId bodyIdB;

	/// The initial target point in world space
	s2Vec2 target;

	/// Stiffness in hertz
	float hertz;

	/// Damping ratio, non-dimensional
	float dampingRatio;
} s2MouseJointDef;

static inline struct s2MouseJointDef s2DefaultMouseJointDef(void)
{
	s2MouseJointDef def = {0};
	def.bodyIdA = s2_nullBodyId;
	def.bodyIdB = s2_nullBodyId;
	def.target = S2_LITERAL(s2Vec2){0.0f, 0.0f};
	def.hertz = 15.0f;
	def.dampingRatio = 1.0f;
	return def;
}

typedef struct s2RevoluteJointDef
{
	/// The first attached body.
	s2BodyId bodyIdA;

	/// The second attached body.
	s2BodyId bodyIdB;

	/// The local anchor point relative to bodyA's origin.
	s2Vec2 localAnchorA;

	/// The local anchor point relative to bodyB's origin.
	s2Vec2 localAnchorB;

	/// The bodyB angle minus bodyA angle in the reference state (radians).
	/// This defines the zero angle for the joint limit.
	float referenceAngle;

	/// A flag to enable joint limits.
	bool enableLimit;

	/// The lower angle for the joint limit (radians).
	float lowerAngle;

	/// The upper angle for the joint limit (radians).
	float upperAngle;

	/// A flag to enable the joint motor.
	bool enableMotor;

	/// The desired motor speed. Usually in radians per second.
	float motorSpeed;

	/// The maximum motor torque used to achieve the desired motor speed.
	/// Usually in N-m.
	float maxMotorTorque;

	float drawSize;

	/// Set this flag to true if the attached bodies should collide.
	bool collideConnected;
} s2RevoluteJointDef;

static inline struct s2RevoluteJointDef s2DefaultRevoluteJointDef(void)
{
	s2RevoluteJointDef def = {0};
	def.bodyIdA = s2_nullBodyId;
	def.bodyIdB = s2_nullBodyId;
	def.localAnchorA = S2_LITERAL(s2Vec2){0.0f, 0.0f};
	def.localAnchorB = S2_LITERAL(s2Vec2){0.0f, 0.0f};
	def.referenceAngle = 0.0f;
	def.lowerAngle = 0.0f;
	def.upperAngle = 0.0f;
	def.maxMotorTorque = 0.0f;
	def.motorSpeed = 0.0f;
	def.enableLimit = false;
	def.enableMotor = false;
	def.drawSize = 1.0f;
	def.collideConnected = false;
	return def;
}
