// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "pool.h"

#include "solver2d/distance.h"
#include "solver2d/geometry.h"
#include "solver2d/types.h"

typedef struct s2BroadPhase s2BroadPhase;

typedef enum s2ShapeType
{
	s2_capsuleShape,
	s2_circleShape,
	s2_polygonShape,
	s2_segmentShape,
	s2_shapeTypeCount
} s2ShapeType;

typedef struct s2Shape
{
	s2Object object;
	int32_t bodyIndex;
	int32_t nextShapeIndex;
	enum s2ShapeType type;
	float density;
	float friction;
	float restitution;
	s2Filter filter;

	s2Box aabb;
	s2Box fatAABB;
	int32_t proxyKey;

	void* userData;
	bool enlargedAABB;

	union
	{
		s2Capsule capsule;
		s2Circle circle;
		s2Polygon polygon;
		s2Segment segment;
	};
} s2Shape;

s2MassData s2Shape_ComputeMass(const s2Shape* shape);
s2Box s2Shape_ComputeAABB(const s2Shape* shape, s2Transform xf);

void s2Shape_CreateProxy(s2Shape* shape, s2BroadPhase* bp, s2BodyType type, s2Transform xf);
void s2Shape_DestroyProxy(s2Shape* shape, s2BroadPhase* bp);

s2DistanceProxy s2Shape_MakeDistanceProxy(const s2Shape* shape);
