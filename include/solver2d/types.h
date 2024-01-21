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
	#define S2_ZERO_INIT                                                                                                         \
		{                                                                                                                        \
		}
#else
	#define S2_LITERAL(T) (T)
	#define S2_ZERO_INIT                                                                                                         \
		{                                                                                                                        \
			0                                                                                                                    \
		}
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
	// sine and cosine
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
	s2_solverTGS_Soft,
	s2_solverTGS_Sticky,
	s2_solverTGS_NGS,
	s2_solverPGS,
	s2_solverTypeCount,
} s2SolverType;

typedef struct s2WorldDef
{
	enum s2SolverType solverType;
} s2WorldDef;

static const s2WorldDef s2_defaultWorldDef = {
	s2_solverPGS_NGS_Block,
};

typedef enum s2BodyType
{
	s2_staticBody = 0,
	s2_kinematicBody = 1,
	s2_dynamicBody = 2,
	s2_bodyTypeCount
} s2BodyType;

typedef struct s2BodyDef
{
	s2BodyType type;
	s2Vec2 position;
	float angle;
	s2Vec2 linearVelocity;
	float angularVelocity;
	float linearDamping;
	float angularDamping;
	float gravityScale;
	void* userData;
} s2BodyDef;

static const s2BodyDef s2_defaultBodyDef = {
	s2_staticBody, // type
	{0.0f, 0.0f},  // position
	0.0f,		   // angle
	{0.0f, 0.0f},  // linearVelocity
	0.0f,		   // angularVelocity
	0.0f,		   // linearDamping
	0.0f,		   // angularDamping
	1.0f,		   // gravity scale
	NULL,		   // userData
};

typedef struct s2Filter
{
	uint32_t categoryBits;
	uint32_t maskBits;
	int32_t groupIndex;
} s2Filter;

static const s2Filter s2_defaultFilter = {0x00000001, 0xFFFFFFFF, 0};

typedef struct s2ShapeDef
{
	void* userData;
	float friction;
	float restitution;
	float density;
	s2Filter filter;
} s2ShapeDef;

static const s2ShapeDef s2_defaultShapeDef = {
	NULL,						 // userData
	0.6f,						 // friction
	0.0f,						 // restitution
	1.0f,						 // density
	{0x00000001, 0xFFFFFFFF, 0}, // filter
};

static inline s2WorldDef s2DefaultWorldDef(void)
{
	s2WorldDef def = S2_ZERO_INIT;
	def.solverType = s2_solverPGS_NGS_Block;
	return def;
}
