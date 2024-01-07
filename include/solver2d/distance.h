// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "constants.h"
#include "types.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct s2SegmentDistanceResult
{
	s2Vec2 closest1;
	s2Vec2 closest2;
	float fraction1;
	float fraction2;
	float distanceSquared;
} s2SegmentDistanceResult;

/// Compute the distance between two line segments, clamping at the end points if needed.
s2SegmentDistanceResult s2SegmentDistance(s2Vec2 p1, s2Vec2 q1, s2Vec2 p2, s2Vec2 q2);

/// A distance proxy is used by the GJK algorithm.
/// It encapsulates any shape.
typedef struct s2DistanceProxy
{
	s2Vec2 vertices[s2_maxPolygonVertices];
	int32_t count;
	float radius;
} s2DistanceProxy;

/// Used to warm start s2Distance.
/// Set count to zero on first call.
typedef struct s2DistanceCache
{
	float metric; ///< length or area
	uint16_t count;
	uint8_t indexA[3]; ///< vertices on shape A
	uint8_t indexB[3]; ///< vertices on shape B
} s2DistanceCache;

static const s2DistanceCache s2_emptyDistanceCache = {0};

/// Input for s2Distance.
/// You have to option to use the shape radii
/// in the computation. Even
typedef struct s2DistanceInput
{
	s2DistanceProxy proxyA;
	s2DistanceProxy proxyB;
	s2Transform transformA;
	s2Transform transformB;
	bool useRadii;
} s2DistanceInput;

/// Output for s2Distance.
typedef struct s2DistanceOutput
{
	s2Vec2 pointA; ///< closest point on shapeA
	s2Vec2 pointB; ///< closest point on shapeB
	float distance;
	int32_t iterations; ///< number of GJK iterations used
} s2DistanceOutput;

/// Compute the closest points between two shapes. Supports any combination of:
/// s2Circle, s2Polygon, s2EdgeShape. The simplex cache is input/output.
/// On the first call set s2SimplexCache.count to zero.
s2DistanceOutput s2ShapeDistance(s2DistanceCache* cache, const s2DistanceInput* input);

s2DistanceProxy s2MakeProxy(const s2Vec2* vertices, int32_t count, float radius);

#ifdef __cplusplus
}
#endif
