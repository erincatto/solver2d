// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include <stdint.h>

/// These ids serve as handles to internal Box2D objects. These should be considered opaque data and passed by value.
/// Include this header if you need the id definitions and not the whole Box2D API.

/// References a world instance
typedef struct s2WorldId
{
	int16_t index;
	uint16_t revision;
} s2WorldId;

/// References a rigid body instance
typedef struct s2BodyId
{
	int32_t index;
	int16_t world;
	uint16_t revision;
} s2BodyId;

/// References a shape instance
typedef struct s2ShapeId
{
	int32_t index;
	int16_t world;
	uint16_t revision;
} s2ShapeId;

/// References a joint instance
typedef struct s2JointId
{
	int32_t index;
	int16_t world;
	uint16_t revision;
} s2JointId;

static const s2WorldId s2_nullWorldId = {-1, 0};
static const s2BodyId s2_nullBodyId = {-1, -1, 0};
static const s2ShapeId s2_nullShapeId = {-1, -1, 0};
static const s2JointId s2_nullJointId = {-1, -1, 0};

#define S2_IS_NULL(ID) (ID.index == -1)
#define S2_NON_NULL(ID) (ID.index != -1)
