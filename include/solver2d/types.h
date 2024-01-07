// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver2d/color.h"
#include "solver2d/constants.h"
#include "solver2d/id.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
#define S2_LITERAL(T) T
#else
#define S2_LITERAL(T) (T)
#endif

#define S2_ARRAY_COUNT(A) (int)(sizeof(A) / sizeof(A[0]))
#define S2_MAYBE_UNUSED(x) ((void)(x))
#define S2_NULL_INDEX (-1)

/// 2D vector
typedef struct s2Vec2
{
	float x, y;
} s2Vec2;

/// 2D rotation
typedef struct s2Rot
{
	/// Sine and cosine
	float s, c;
} s2Rot;

/// A 2D rigid transform
typedef struct s2Transform
{
	s2Vec2 p;
	s2Rot q;
} s2Transform;

/// A 2-by-2 Matrix
typedef struct s2Mat22
{
	/// columns
	s2Vec2 cx, cy;
} s2Mat22;

/// Axis-aligned bounding box
typedef struct s2Box
{
	s2Vec2 lowerBound;
	s2Vec2 upperBound;
} s2Box;

/// Ray-cast input data. The ray extends from p1 to p1 + maxFraction * (p2 - p1).
typedef struct s2RayCastInput
{
	s2Vec2 p1, p2;
	float maxFraction;
} s2RayCastInput;

/// Ray-cast output data. The ray hits at p1 + fraction * (p2 - p1), where p1 and p2 come from s2RayCastInput.
typedef struct s2RayCastOutput
{
	s2Vec2 normal;
	s2Vec2 point;
	float fraction;
	int32_t iterations;
	bool hit;
} s2RayCastOutput;

typedef enum s2SolverType
{
	s2_solverPGS_NGS,
	s2_solverPGS_Soft,
	s2_solverTGS_Soft,
	s2_solverTGS_Soft_Sticky,
} s2SolverType;

typedef struct s2WorldDef
{
	s2SolverType solverType;
	s2Vec2 gravity;
	float restitutionThreshold;
	bool enableSleep;
	int32_t bodyCapacity;
	int32_t shapeCapacity;
	int32_t contactCapacity;
	int32_t jointCapacity;
	int32_t stackAllocatorCapacity;
} s2WorldDef;

/// The body type.
/// static: zero mass, zero velocity, may be manually moved
/// kinematic: zero mass, non-zero velocity set by user, moved by solver
/// dynamic: positive mass, non-zero velocity determined by forces, moved by solver
typedef enum s2BodyType
{
	s2_staticBody = 0,
	s2_kinematicBody = 1,
	s2_dynamicBody = 2,
	s2_bodyTypeCount
} s2BodyType;

/// A body definition holds all the data needed to construct a rigid body.
/// You can safely re-use body definitions. Shapes are added to a body after construction.
typedef struct s2BodyDef
{
	/// The body type: static, kinematic, or dynamic.
	/// Note: if a dynamic body would have zero mass, the mass is set to one.
	enum s2BodyType type;

	/// The world position of the body. Avoid creating bodies at the origin
	/// since this can lead to many overlapping shapes.
	s2Vec2 position;

	/// The world angle of the body in radians.
	float angle;

	/// The linear velocity of the body's origin in world co-ordinates.
	s2Vec2 linearVelocity;

	/// The angular velocity of the body.
	float angularVelocity;

	/// Linear damping is use to reduce the linear velocity. The damping parameter
	/// can be larger than 1.0f but the damping effect becomes sensitive to the
	/// time step when the damping parameter is large.
	float linearDamping;

	/// Angular damping is use to reduce the angular velocity. The damping parameter
	/// can be larger than 1.0f but the damping effect becomes sensitive to the
	/// time step when the damping parameter is large.
	float angularDamping;

	/// Scale the gravity applied to this body.
	float gravityScale;

	/// Use this to store application specific body data.
	void* userData;

} s2BodyDef;

/// Used to create a shape
typedef struct s2ShapeDef
{
	/// Use this to store application specific shape data.
	void* userData;

	/// The friction coefficient, usually in the range [0,1].
	float friction;

	/// The restitution (elasticity) usually in the range [0,1].
	float restitution;

	/// The density, usually in kg/m^2.
	float density;
} s2ShapeDef;

/// Make a world definition with default values.
static inline s2WorldDef s2DefaultWorldDef(void)
{
	s2WorldDef def = {0};
	def.solverType = s2_solverPGS_NGS;
	def.gravity = S2_LITERAL(s2Vec2){0.0f, -10.0f};
	def.restitutionThreshold = 1.0f;
	def.bodyCapacity = 8;
	def.shapeCapacity = 8;
	def.contactCapacity = 8;
	def.jointCapacity = 8;
	def.stackAllocatorCapacity = 1024 * 1024;
	return def;
}

/// Make a body definition with default values.
static inline s2BodyDef s2DefaultBodyDef(void)
{
	s2BodyDef def = {0};
	def.type = s2_staticBody;
	def.position = S2_LITERAL(s2Vec2){0.0f, 0.0f};
	def.angle = 0.0f;
	def.linearVelocity = S2_LITERAL(s2Vec2){0.0f, 0.0f};
	def.angularVelocity = 0.0f;
	def.linearDamping = 0.0f;
	def.angularDamping = 0.0f;
	def.gravityScale = 1.0f;
	def.userData = NULL;
	return def;
}

static inline struct s2ShapeDef s2DefaultShapeDef(void)
{
	s2ShapeDef def = {0};
	def.friction = 0.6f;
	def.restitution = 0.0f;
	def.density = 0.0f;
	return def;
}
