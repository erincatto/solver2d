// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver2d/constants.h"
#include "solver2d/types.h"

typedef struct s2Hull s2Hull;
typedef struct s2RayCastOutput s2RayCastOutput;
typedef struct s2RayCastInput s2RayCastInput;

/// This holds the mass data computed for a shape.
typedef struct s2MassData
{
	/// The mass of the shape, usually in kilograms.
	float mass;

	/// The position of the shape's centroid relative to the shape's origin.
	s2Vec2 center;

	/// The rotational inertia of the shape about the local origin.
	float I;
} s2MassData;

/// A solid circle
typedef struct s2Circle
{
	s2Vec2 point;
	float radius;
} s2Circle;

/// A solid capsule
typedef struct s2Capsule
{
	s2Vec2 point1, point2;
	float radius;
} s2Capsule;

/// A solid convex polygon. It is assumed that the interior of the polygon is to
/// the left of each edge.
/// Polygons have a maximum number of vertices equal to s2_maxPolygonVertices.
/// In most cases you should not need many vertices for a convex polygon.
typedef struct s2Polygon
{
	s2Vec2 vertices[s2_maxPolygonVertices];
	s2Vec2 normals[s2_maxPolygonVertices];
	float radius;
	int32_t count;
} s2Polygon;

/// A line segment with two-sided collision.
typedef struct s2Segment
{
	s2Vec2 point1, point2;
} s2Segment;

/// A smooth line segment with one-sided collision. Only collides on the right side.
/// Normally these are generated from a chain shape.
/// ghost1 -> point1 -> point2 -> ghost2
/// This is only relevant for contact manifolds, otherwise use a regular segment.
typedef struct s2SmoothSegment
{
	/// The tail ghost vertex
	s2Vec2 ghost1;

	/// The line segment
	s2Vec2 point1, point2;

	/// The head ghost vertex
	s2Vec2 ghost2;
} s2SmoothSegment;

#ifdef __cplusplus
extern "C"
{
#endif

bool s2IsValidRay(const s2RayCastInput* input);

/// Helper functions to make convex polygons
s2Polygon s2MakePolygon(const s2Hull* hull);
s2Polygon s2MakeSquare(float h);
s2Polygon s2MakeBox(float hx, float hy);
s2Polygon s2MakeOffsetBox(float hx, float hy, s2Vec2 center, float angle);
s2Polygon s2MakeCapsule(s2Vec2 p1, s2Vec2 p2, float radius);

/// Compute mass properties
s2MassData s2ComputeCircleMass(const s2Circle* shape, float density);
s2MassData s2ComputeCapsuleMass(const s2Capsule* shape, float density);
s2MassData s2ComputePolygonMass(const s2Polygon* shape, float density);

/// These compute the bounding box in world space
s2Box s2ComputeCircleAABB(const s2Circle* shape, s2Transform xf);
s2Box s2ComputeCapsuleAABB(const s2Capsule* shape, s2Transform xf);
s2Box s2ComputePolygonAABB(const s2Polygon* shape, s2Transform xf);
s2Box s2ComputeSegmentAABB(const s2Segment* shape, s2Transform xf);

/// Test a point in local space
bool s2PointInCircle(s2Vec2 point, const s2Circle* shape);
bool s2PointInCapsule(s2Vec2 point, const s2Capsule* shape);
bool s2PointInPolygon(s2Vec2 point, const s2Polygon* shape);

// Ray cast versus shape in shape local space. Initial overlap is treated as a miss.
s2RayCastOutput s2RayCastCircle(const s2RayCastInput* input, const s2Circle* shape);
s2RayCastOutput s2RayCastCapsule(const s2RayCastInput* input, const s2Capsule* shape);
s2RayCastOutput s2RayCastSegment(const s2RayCastInput* input, const s2Segment* shape);
s2RayCastOutput s2RayCastPolygon(const s2RayCastInput* input, const s2Polygon* shape);

#ifdef __cplusplus
}
#endif
