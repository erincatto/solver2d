// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver2d/types.h"

#define s2_nullFeature UCHAR_MAX

typedef struct s2Circle s2Circle;
typedef struct s2Capsule s2Capsule;
typedef struct s2DistanceCache s2DistanceCache;
typedef struct s2Polygon s2Polygon;
typedef struct s2Segment s2Segment;
typedef struct s2SmoothSegment s2SmoothSegment;

#define S2_MAKE_ID(A, B) ((uint8_t)(A) << 8 | (uint8_t)(B))

typedef struct s2ManifoldPoint
{
	// world coordinates of contact point
	// todo accuracy problem
	//s2Vec2 point;

	// todo more accurate
	s2Vec2 localAnchorA;
	s2Vec2 localAnchorB;

	// Friction anchors
	s2Vec2 frictionAnchorA;
	s2Vec2 frictionAnchorB;
	s2Vec2 frictionNormalA;
	s2Vec2 frictionNormalB;

	float separation;
	float normalImpulse;
	float tangentImpulse;
	uint16_t id;

	bool persisted;
} s2ManifoldPoint;

typedef struct s2Manifold
{
	s2ManifoldPoint points[2];
	s2Vec2 normal;
	int32_t pointCount;
	int32_t constraintIndex;
	bool frictionPersisted;
} s2Manifold;

static const s2Manifold s2_emptyManifold = {0};

#ifdef __cplusplus
extern "C"
{
#endif

/// Compute the collision manifold between two circles.
s2Manifold s2CollideCircles(const s2Circle* circleA, s2Transform xfA, const s2Circle* circleB, s2Transform xfB);

/// Compute the collision manifold between a capsule and circle
s2Manifold s2CollideCapsuleAndCircle(const s2Capsule* capsuleA, s2Transform xfA, const s2Circle* circleB, s2Transform xfB);

/// Compute the collision manifold between an segment and a circle.
s2Manifold s2CollideSegmentAndCircle(const s2Segment* segmentA, s2Transform xfA, const s2Circle* circleB, s2Transform xfB);

/// Compute the collision manifold between a polygon and a circle.
s2Manifold s2CollidePolygonAndCircle(const s2Polygon* polygonA, s2Transform xfA, const s2Circle* circleB, s2Transform xfB);

/// Compute the collision manifold between a capsule and circle
s2Manifold s2CollideCapsules(const s2Capsule* capsuleA, s2Transform xfA, const s2Capsule* capsuleB, s2Transform xfB,
							 s2DistanceCache* cache);

/// Compute the collision manifold between an segment and a capsule.
s2Manifold s2CollideSegmentAndCapsule(const s2Segment* segmentA, s2Transform xfA, const s2Capsule* capsuleB, s2Transform xfB,
									  s2DistanceCache* cache);

/// Compute the collision manifold between a polygon and capsule
s2Manifold s2CollidePolygonAndCapsule(const s2Polygon* polygonA, s2Transform xfA, const s2Capsule* capsuleB, s2Transform xfB,
									  s2DistanceCache* cache);

/// Compute the collision manifold between two polygons.
s2Manifold s2CollidePolygons(const s2Polygon* polyA, s2Transform xfA, const s2Polygon* polyB, s2Transform xfB,
							 s2DistanceCache* cache);

/// Compute the collision manifold between an segment and a polygon.
s2Manifold s2CollideSegmentAndPolygon(const s2Segment* segmentA, s2Transform xfA, const s2Polygon* polygonB, s2Transform xfB,
									  s2DistanceCache* cache);

#ifdef __cplusplus
}
#endif
