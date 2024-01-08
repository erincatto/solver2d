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
#define S2_ZERO_INIT {}
#else
#define S2_LITERAL(T) (T)
#define S2_ZERO_INIT {0}
#endif

#define S2_ARRAY_COUNT(A) (int)(sizeof(A) / sizeof(A[0]))
#define S2_MAYBE_UNUSED(x) ((void)(x))
#define S2_NULL_INDEX (-1)

typedef struct s2Vec2
{
	float x, y;
} s2Vec2;

typedef struct s2Rot
{
	/// Sine and cosine
	float s, c;
} s2Rot;

typedef struct s2Transform
{
	s2Vec2 p;
	s2Rot q;
} s2Transform;

typedef struct s2Mat22
{
	/// columns
	s2Vec2 cx, cy;
} s2Mat22;

typedef struct s2Box
{
	s2Vec2 lowerBound;
	s2Vec2 upperBound;
} s2Box;

typedef struct s2RayCastInput
{
	s2Vec2 p1, p2;
	float maxFraction;
} s2RayCastInput;

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
	s2_solverPGS_NGS_Block = 0,
	s2_solverPGS_NGS,
	s2_solverPGS_Soft,
	s2_solverXPBD,
	//s2_solverTGS_Soft,
	//s2_solverTGS_Soft_Sticky,
	s2_solverTypeCount,
} s2SolverType;

typedef struct s2WorldDef
{
	enum s2SolverType solverType;
	s2Vec2 gravity;
	float restitutionThreshold;
	bool enableSleep;
	int32_t bodyCapacity;
	int32_t shapeCapacity;
	int32_t contactCapacity;
	int32_t jointCapacity;
	int32_t stackAllocatorCapacity;
} s2WorldDef;

typedef enum s2BodyType
{
	s2_staticBody = 0,
	s2_kinematicBody = 1,
	s2_dynamicBody = 2,
	s2_bodyTypeCount
} s2BodyType;

typedef struct s2BodyDef
{
	enum s2BodyType type;
	s2Vec2 position;
	float angle;
	s2Vec2 linearVelocity;
	float angularVelocity;
	void* userData;

} s2BodyDef;

typedef struct s2ShapeDef
{
	void* userData;
	float friction;
	float restitution;
	float density;
} s2ShapeDef;

static inline s2WorldDef s2DefaultWorldDef(void)
{
	s2WorldDef def = S2_ZERO_INIT;
	def.solverType = s2_solverPGS_NGS_Block;
	def.gravity = S2_LITERAL(s2Vec2){0.0f, -10.0f};
	def.restitutionThreshold = 1.0f;
	def.bodyCapacity = 8;
	def.shapeCapacity = 8;
	def.contactCapacity = 8;
	def.jointCapacity = 8;
	def.stackAllocatorCapacity = 1024 * 1024;
	return def;
}

static inline s2BodyDef s2DefaultBodyDef(void)
{
	s2BodyDef def = S2_ZERO_INIT;
	def.type = s2_staticBody;
	def.position = S2_LITERAL(s2Vec2){0.0f, 0.0f};
	def.angle = 0.0f;
	def.linearVelocity = S2_LITERAL(s2Vec2){0.0f, 0.0f};
	def.angularVelocity = 0.0f;
	def.userData = NULL;
	return def;
}

static inline struct s2ShapeDef s2DefaultShapeDef(void)
{
	s2ShapeDef def = S2_ZERO_INIT;
	def.friction = 0.6f;
	def.restitution = 0.0f;
	def.density = 1.0f;
	return def;
}
