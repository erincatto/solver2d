// SPDX-FileCopyrightText: 2023 Erin Catto
// SPDX-License-Identifier: MIT

#include "solver2d/manifold.h"

#include "solver2d/distance.h"
#include "solver2d/geometry.h"
#include "solver2d/math.h"
#include "core.h"

#include <float.h>
#include <string.h>

#define S2_MAKE_ID(A, B) ((uint8_t)(A) << 8 | (uint8_t)(B))

s2Manifold s2CollideCircles(const s2Circle* circleA, s2Transform xfA, const s2Circle* circleB, s2Transform xfB)
{
	s2Manifold manifold = {0};

	s2Transform xf = s2InvMulTransforms(xfA, xfB);

	s2Vec2 pointA = circleA->point;
	s2Vec2 pointB = s2TransformPoint(xf, circleB->point);

	float distance;
	s2Vec2 normal = s2GetLengthAndNormalize(&distance, s2Sub(pointB, pointA));

	float radiusA = circleA->radius;
	float radiusB = circleB->radius;

	float separation = distance - radiusA - radiusB;
	if (separation > s2_speculativeDistance)
	{
		return manifold;
	}

	s2Vec2 cA = s2MulAdd(pointA, radiusA, normal);
	s2Vec2 cB = s2MulAdd(pointB, -radiusB, normal);
	s2Vec2 contactPointA = s2Lerp(cA, cB, 0.5f);

	manifold.normal = s2RotateVector(xfA.q, normal);
	manifold.points[0].localAnchorA = contactPointA;
	manifold.points[0].localAnchorB = s2InvTransformPoint(xf, contactPointA);
	manifold.points[0].separation = separation;
	manifold.points[0].id = 0;
	manifold.pointCount = 1;
	return manifold;
}

/// Compute the collision manifold between a capsule and circle
s2Manifold s2CollideCapsuleAndCircle(const s2Capsule* capsuleA, s2Transform xfA, const s2Circle* circleB, s2Transform xfB)
{
	s2Manifold manifold = {0};

	s2Transform xf = s2InvMulTransforms(xfA, xfB);

	// Compute circle position in the frame of the capsule.
	s2Vec2 pB = s2TransformPoint(xf, circleB->point);

	// Compute closest point
	s2Vec2 p1 = capsuleA->point1;
	s2Vec2 p2 = capsuleA->point2;

	s2Vec2 e = s2Sub(p2, p1);

	// dot(p - pA, e) = 0
	// dot(p - (p1 + s1 * e), e) = 0
	// s1 = dot(p - p1, e)
	s2Vec2 pA;
	float s1 = s2Dot(s2Sub(pB, p1), e);
	float s2 = s2Dot(s2Sub(p2, pB), e);
	if (s1 < 0.0f)
	{
		// p1 region
		pA = p1;
	}
	else if (s2 < 0.0f)
	{
		// p2 region
		pA = p2;
	}
	else
	{
		// circle colliding with segment interior
		float s = s1 / s2Dot(e, e);
		pA = s2MulAdd(p1, s, e);
	}

	float distance;
	s2Vec2 normal = s2GetLengthAndNormalize(&distance, s2Sub(pB, pA));

	float radiusA = capsuleA->radius;
	float radiusB = circleB->radius;
	float separation = distance - radiusA - radiusB;
	if (separation > s2_speculativeDistance)
	{
		return manifold;
	}

	s2Vec2 cA = s2MulAdd(pA, radiusA, normal);
	s2Vec2 cB = s2MulAdd(pB, -radiusB, normal);
	s2Vec2 contactPointA = s2Lerp(cA, cB, 0.5f);

	manifold.normal = s2RotateVector(xfA.q, normal);
	manifold.points[0].localAnchorA = contactPointA;
	manifold.points[0].localAnchorB = s2InvTransformPoint(xf, contactPointA);
	manifold.points[0].separation = separation;
	manifold.points[0].id = 0;
	manifold.pointCount = 1;
	return manifold;
}

s2Manifold s2CollidePolygonAndCircle(const s2Polygon* polygonA, s2Transform xfA, const s2Circle* circleB, s2Transform xfB)
{
	s2Manifold manifold = {0};

	s2Transform xf = s2InvMulTransforms(xfA, xfB);

	// Compute circle position in the frame of the polygon.
	s2Vec2 c = s2TransformPoint(xf, circleB->point);
	float radiusA = polygonA->radius;
	float radiusB = circleB->radius;
	float radius = radiusA + radiusB;

	// Find the min separating edge.
	int normalIndex = 0;
	float separation = -FLT_MAX;
	int vertexCount = polygonA->count;
	const s2Vec2* vertices = polygonA->vertices;
	const s2Vec2* normals = polygonA->normals;

	for (int i = 0; i < vertexCount; ++i)
	{
		float s = s2Dot(normals[i], s2Sub(c, vertices[i]));
		if (s > separation)
		{
			separation = s;
			normalIndex = i;
		}
	}

	if (separation > radius + s2_speculativeDistance)
	{
		return manifold;
	}

	// Vertices of the reference edge.
	int vertIndex1 = normalIndex;
	int vertIndex2 = vertIndex1 + 1 < vertexCount ? vertIndex1 + 1 : 0;
	s2Vec2 v1 = vertices[vertIndex1];
	s2Vec2 v2 = vertices[vertIndex2];

	// Compute barycentric coordinates
	float u1 = s2Dot(s2Sub(c, v1), s2Sub(v2, v1));
	float u2 = s2Dot(s2Sub(c, v2), s2Sub(v1, v2));

	if (u1 < 0.0f && separation > FLT_EPSILON)
	{
		// Circle center is closest to v1 and safely outside the polygon
		s2Vec2 normal = s2Normalize(s2Sub(c, v1));
		separation = s2Dot(s2Sub(c, v1), normal);
		if (separation > radius + s2_speculativeDistance)
		{
			return manifold;
		}

		s2Vec2 cA = s2MulAdd(v1, radiusA, normal);
		s2Vec2 cB = s2MulSub(c, radiusB, normal);
		s2Vec2 contactPointA = s2Lerp(cA, cB, 0.5f);

		manifold.normal = s2RotateVector(xfA.q, normal);
		manifold.points[0].localAnchorA = contactPointA;
		manifold.points[0].localAnchorB = s2InvTransformPoint(xf, contactPointA);
		manifold.points[0].separation = s2Dot(s2Sub(cB, cA), normal);
		manifold.points[0].id = 0;
		manifold.pointCount = 1;
	}
	else if (u2 < 0.0f && separation > FLT_EPSILON)
	{
		// Circle center is closest to v2 and safely outside the polygon
		s2Vec2 normal = s2Normalize(s2Sub(c, v2));
		separation = s2Dot(s2Sub(c, v2), normal);
		if (separation > radius + s2_speculativeDistance)
		{
			return manifold;
		}

		s2Vec2 cA = s2MulAdd(v2, radiusA, normal);
		s2Vec2 cB = s2MulSub(c, radiusB, normal);
		s2Vec2 contactPointA = s2Lerp(cA, cB, 0.5f);

		manifold.normal = s2RotateVector(xfA.q, normal);
		manifold.points[0].localAnchorA = contactPointA;
		manifold.points[0].localAnchorB = s2InvTransformPoint(xf, contactPointA);
		manifold.points[0].separation = s2Dot(s2Sub(cB, cA), normal);
		manifold.points[0].id = 0;
		manifold.pointCount = 1;
	}
	else
	{
		// Circle center is between v1 and v2. Center may be inside polygon
		s2Vec2 normal = normals[normalIndex];
		manifold.normal = s2RotateVector(xfA.q, normal);

		// cA is the projection of the circle center onto to the reference edge
		s2Vec2 cA = s2MulAdd(c, radiusA - s2Dot(s2Sub(c, v1), normal), normal);

		// cB is the deepest point on the circle with respect to the reference edge
		s2Vec2 cB = s2MulSub(c, radiusB, normal);

		s2Vec2 contactPointA = s2Lerp(cA, cB, 0.5f);

		// The contact point is the midpoint in world space
		manifold.points[0].localAnchorA = contactPointA;
		manifold.points[0].localAnchorB = s2InvTransformPoint(xf, contactPointA);
		manifold.points[0].separation = separation - radius;
		manifold.points[0].id = 0;
		manifold.pointCount = 1;
	}

	return manifold;
}

s2Manifold s2CollideCapsules(const s2Capsule* capsuleA, s2Transform xfA, const s2Capsule* capsuleB, s2Transform xfB,
							 s2DistanceCache* cache)
{
	s2Polygon polyA = s2MakeCapsule(capsuleA->point1, capsuleA->point2, capsuleA->radius);
	s2Polygon polyB = s2MakeCapsule(capsuleB->point1, capsuleB->point2, capsuleB->radius);
	return s2CollidePolygons(&polyA, xfA, &polyB, xfB, cache);
}

s2Manifold s2CollideSegmentAndCapsule(const s2Segment* segmentA, s2Transform xfA, const s2Capsule* capsuleB, s2Transform xfB,
									  s2DistanceCache* cache)
{
	s2Polygon polyA = s2MakeCapsule(segmentA->point1, segmentA->point2, 0.0f);
	s2Polygon polyB = s2MakeCapsule(capsuleB->point1, capsuleB->point2, capsuleB->radius);
	return s2CollidePolygons(&polyA, xfA, &polyB, xfB, cache);
}

s2Manifold s2CollidePolygonAndCapsule(const s2Polygon* polygonA, s2Transform xfA, const s2Capsule* capsuleB, s2Transform xfB,
									  s2DistanceCache* cache)
{
	s2Polygon polyB = s2MakeCapsule(capsuleB->point1, capsuleB->point2, capsuleB->radius);
	return s2CollidePolygons(polygonA, xfA, &polyB, xfB, cache);
}

// Polygon clipper used by GJK and SAT to compute contact points when there are potentially two contact points.
static s2Manifold s2ClipPolygons(const s2Polygon* polyA, const s2Polygon* polyB, int edgeA, int edgeB, bool flip)
{
	s2Manifold manifold = {0};

	// reference polygon
	const s2Polygon* poly1;
	int i11, i12;

	// incident polygon
	const s2Polygon* poly2;
	int i21, i22;

	if (flip)
	{
		poly1 = polyB;
		poly2 = polyA;
		i11 = edgeB;
		i12 = edgeB + 1 < polyB->count ? edgeB + 1 : 0;
		i21 = edgeA;
		i22 = edgeA + 1 < polyA->count ? edgeA + 1 : 0;
	}
	else
	{
		poly1 = polyA;
		poly2 = polyB;
		i11 = edgeA;
		i12 = edgeA + 1 < polyA->count ? edgeA + 1 : 0;
		i21 = edgeB;
		i22 = edgeB + 1 < polyB->count ? edgeB + 1 : 0;
	}

	s2Vec2 normal = poly1->normals[i11];

	// Reference edge vertices
	s2Vec2 v11 = poly1->vertices[i11];
	s2Vec2 v12 = poly1->vertices[i12];

	// Incident edge vertices
	s2Vec2 v21 = poly2->vertices[i21];
	s2Vec2 v22 = poly2->vertices[i22];

	s2Vec2 tangent = s2CrossSV(1.0f, normal);

	float lower1 = 0.0f;
	float upper1 = s2Dot(s2Sub(v12, v11), tangent);

	// Incident edge points opposite of tangent due to CCW winding
	float upper2 = s2Dot(s2Sub(v21, v11), tangent);
	float lower2 = s2Dot(s2Sub(v22, v11), tangent);

	// This check can fail slightly due to mismatch with GJK code.
	// Perhaps fallback to a single point here? Otherwise we get two coincident points.
	// if (upper2 < lower1 || upper1 < lower2)
	//{
	//	// numeric failure
	//	S2_ASSERT(false);
	//	return manifold;
	//}

	s2Vec2 vLower;
	if (lower2 < lower1 && upper2 - lower2 > FLT_EPSILON)
	{
		vLower = s2Lerp(v22, v21, (lower1 - lower2) / (upper2 - lower2));
	}
	else
	{
		vLower = v22;
	}

	s2Vec2 vUpper;
	if (upper2 > upper1 && upper2 - lower2 > FLT_EPSILON)
	{
		vUpper = s2Lerp(v22, v21, (upper1 - lower2) / (upper2 - lower2));
	}
	else
	{
		vUpper = v21;
	}

	// TODO_ERIN vLower can be very close to vUpper, reduce to one point?

	float separationLower = s2Dot(s2Sub(vLower, v11), normal);
	float separationUpper = s2Dot(s2Sub(vUpper, v11), normal);

	float r1 = poly1->radius;
	float r2 = poly2->radius;

	// Put contact points at midpoint, accounting for radii
	vLower = s2MulAdd(vLower, 0.5f * (r1 - r2 - separationLower), normal);
	vUpper = s2MulAdd(vUpper, 0.5f * (r1 - r2 - separationUpper), normal);

	float radius = r1 + r2;

	if (flip == false)
	{
		manifold.normal = normal;
		s2ManifoldPoint* cp = manifold.points + 0;

		{
			cp->localAnchorA = vLower;
			cp->separation = separationLower - radius;
			if (cp->separation < -0.5f)
			{
				cp->separation += 0.0f;
			}
			cp->id = S2_MAKE_ID(i11, i22);
			manifold.pointCount += 1;
			cp += 1;
		}

		{
			cp->localAnchorA = vUpper;
			cp->separation = separationUpper - radius;
			if (cp->separation < -0.5f)
			{
				cp->separation += 0.0f;
			}
			cp->id = S2_MAKE_ID(i12, i21);
			manifold.pointCount += 1;
		}
	}
	else
	{
		manifold.normal = s2Neg(normal);
		s2ManifoldPoint* cp = manifold.points + 0;

		{
			cp->localAnchorA = vUpper;
			cp->separation = separationUpper - radius;
			if (cp->separation < -0.5f)
			{
				cp->separation += 0.0f;
			}
			cp->id = S2_MAKE_ID(i21, i12);
			manifold.pointCount += 1;
			cp += 1;
		}

		{
			cp->localAnchorA = vLower;
			cp->separation = separationLower - radius;
			if (cp->separation < -0.5f)
			{
				cp->separation += 0.0f;
			}
			cp->id = S2_MAKE_ID(i22, i11);
			manifold.pointCount += 1;
		}
	}

	return manifold;
}

// Find the max separation between poly1 and poly2 using edge normals from poly1.
static float s2FindMaxSeparation(int* edgeIndex, const s2Polygon* poly1, const s2Polygon* poly2)
{
	int count1 = poly1->count;
	int count2 = poly2->count;
	const s2Vec2* n1s = poly1->normals;
	const s2Vec2* v1s = poly1->vertices;
	const s2Vec2* v2s = poly2->vertices;

	int bestIndex = 0;
	float maxSeparation = -FLT_MAX;
	for (int i = 0; i < count1; ++i)
	{
		// Get poly1 normal in frame2.
		s2Vec2 n = n1s[i];
		s2Vec2 v1 = v1s[i];

		// Find deepest point for normal i.
		float si = FLT_MAX;
		for (int j = 0; j < count2; ++j)
		{
			float sij = s2Dot(n, s2Sub(v2s[j], v1));
			if (sij < si)
			{
				si = sij;
			}
		}

		if (si > maxSeparation)
		{
			maxSeparation = si;
			bestIndex = i;
		}
	}

	*edgeIndex = bestIndex;
	return maxSeparation;
}

// This function assumes there is overlap
static s2Manifold s2PolygonSAT(const s2Polygon* polyA, const s2Polygon* polyB)
{
	int edgeA = 0;
	float separationA = s2FindMaxSeparation(&edgeA, polyA, polyB);

	int edgeB = 0;
	float separationB = s2FindMaxSeparation(&edgeB, polyB, polyA);

	bool flip;

	if (separationB > separationA)
	{
		flip = true;
		s2Vec2 searchDirection = polyB->normals[edgeB];

		// Find the incident edge on polyA
		int count = polyA->count;
		const s2Vec2* normals = polyA->normals;
		edgeA = 0;
		float minDot = FLT_MAX;
		for (int i = 0; i < count; ++i)
		{
			float dot = s2Dot(searchDirection, normals[i]);
			if (dot < minDot)
			{
				minDot = dot;
				edgeA = i;
			}
		}
	}
	else
	{
		flip = false;
		s2Vec2 searchDirection = polyA->normals[edgeA];

		// Find the incident edge on polyB
		int count = polyB->count;
		const s2Vec2* normals = polyB->normals;
		edgeB = 0;
		float minDot = FLT_MAX;
		for (int i = 0; i < count; ++i)
		{
			float dot = s2Dot(searchDirection, normals[i]);
			if (dot < minDot)
			{
				minDot = dot;
				edgeB = i;
			}
		}
	}

	return s2ClipPolygons(polyA, polyB, edgeA, edgeB, flip);
}

// Due to speculation, every polygon is rounded
// Algorithm:
// compute distance
// if distance <= 0.1f * s2_linearSlop
//   SAT
// else
//   find closest features from GJK
//   expect 2-1 or 1-1 or 1-2 features
//   if 2-1 or 1-2
//     clip
//   else
//     vertex-vertex
//   end
// end
s2Manifold s2CollidePolygons(const s2Polygon* polyA, s2Transform xfA, const s2Polygon* polyB, s2Transform xfB,
							 s2DistanceCache* cache)
{
	s2Manifold manifold = {0};
	float radius = polyA->radius + polyB->radius;

	s2Transform xf = s2InvMulTransforms(xfA, xfB);

	// Put polyB in polyA's frame to reduce round-off error
	s2Polygon localPolyB;
	localPolyB.count = polyB->count;
	localPolyB.radius = polyB->radius;
	for (int i = 0; i < localPolyB.count; ++i)
	{
		localPolyB.vertices[i] = s2TransformPoint(xf, polyB->vertices[i]);
		localPolyB.normals[i] = s2RotateVector(xf.q, polyB->normals[i]);
	}

	s2DistanceInput input;
	input.proxyA = s2MakeProxy(polyA->vertices, polyA->count, 0.0f);
	input.proxyB = s2MakeProxy(localPolyB.vertices, localPolyB.count, 0.0f);
	input.transformA = s2Transform_identity;
	input.transformB = s2Transform_identity;
	input.useRadii = false;

	s2DistanceOutput output = s2ShapeDistance(cache, &input);

	if (output.distance > radius + s2_speculativeDistance)
	{
		// no contact
		return manifold;
	}

	if (output.distance < 0.1f * s2_linearSlop)
	{
		// distance is small or zero, fallback to SAT
		manifold = s2PolygonSAT(polyA, &localPolyB);

		if (manifold.pointCount > 0)
		{
			manifold.normal = s2RotateVector(xfA.q, manifold.normal);
			for (int i = 0; i < manifold.pointCount; ++i)
			{
				manifold.points[i].localAnchorB = s2InvTransformPoint(xf, manifold.points[i].localAnchorA);
			}
		}

		return manifold;
	}

	if (cache->count == 1)
	{
		// vertex-vertex collision
		s2Vec2 pA = output.pointA;
		s2Vec2 pB = output.pointB;

		float distance = output.distance;
		s2Vec2 normal = s2Normalize(s2Sub(pB, pA));
		s2Vec2 contactPointA = s2MulAdd(pB, 0.5f * (polyA->radius - localPolyB.radius - distance), normal);

		manifold.normal = s2RotateVector(xfA.q, normal);
		s2ManifoldPoint* cp = manifold.points + 0;
		cp->localAnchorA = contactPointA;
		cp->localAnchorB = s2InvTransformPoint(xf, contactPointA);
		cp->separation = distance - radius;
		cp->id = S2_MAKE_ID(cache->indexA[0], cache->indexB[0]);
		manifold.pointCount = 1;
		return manifold;
	}

	// vertex-edge collision
	S2_ASSERT(cache->count == 2);
	bool flip;
	int countA = polyA->count;
	int countB = localPolyB.count;
	int edgeA, edgeB;

	int a1 = cache->indexA[0];
	int a2 = cache->indexA[1];
	int b1 = cache->indexB[0];
	int s2x = cache->indexB[1];

	if (a1 == a2)
	{
		// 1 point on A, expect 2 points on B
		S2_ASSERT(b1 != s2x);

		// Find reference edge that most aligns with vector between closest points.
		// This works for capsules and polygons
		s2Vec2 axis = s2Sub(output.pointA, output.pointB);
		float dot1 = s2Dot(axis, localPolyB.normals[b1]);
		float dot2 = s2Dot(axis, localPolyB.normals[s2x]);
		edgeB = dot1 > dot2 ? b1 : s2x;

		flip = true;

		// Get the normal of the reference edge in polyA's frame.
		axis = localPolyB.normals[edgeB];

		// Find the incident edge on polyA
		// Limit search to edges adjacent to closest vertex on A
		int edgeA1 = a1;
		int edgeA2 = edgeA1 == 0 ? countA - 1 : edgeA1 - 1;
		dot1 = s2Dot(axis, polyA->normals[edgeA1]);
		dot2 = s2Dot(axis, polyA->normals[edgeA2]);
		edgeA = dot1 < dot2 ? edgeA1 : edgeA2;
	}
	else
	{
		// Find reference edge that most aligns with vector between closest points.
		// This works for capsules and polygons
		s2Vec2 axis = s2Sub(output.pointB, output.pointA);
		float dot1 = s2Dot(axis, polyA->normals[a1]);
		float dot2 = s2Dot(axis, polyA->normals[a2]);
		edgeA = dot1 > dot2 ? a1 : a2;

		flip = false;

		// Get the normal of the reference edge in polyB's frame.
		axis = polyA->normals[edgeA];

		// Find the incident edge on polyB
		// Limit search to edges adjacent to closest vertex
		int edgeB1 = b1;
		int edgeB2 = edgeB1 == 0 ? countB - 1 : edgeB1 - 1;
		dot1 = s2Dot(axis, localPolyB.normals[edgeB1]);
		dot2 = s2Dot(axis, localPolyB.normals[edgeB2]);
		edgeB = dot1 < dot2 ? edgeB1 : edgeB2;
	}

	manifold = s2ClipPolygons(polyA, &localPolyB, edgeA, edgeB, flip);
	if (manifold.pointCount > 0)
	{
		manifold.normal = s2RotateVector(xfA.q, manifold.normal);
		for (int i = 0; i < manifold.pointCount; ++i)
		{
			manifold.points[i].localAnchorB = s2InvTransformPoint(xf, manifold.points[i].localAnchorA);
		}
	}

	return manifold;
}

s2Manifold s2CollideSegmentAndCircle(const s2Segment* segmentA, s2Transform xfA, const s2Circle* circleB, s2Transform xfB)
{
	s2Capsule capsuleA = {segmentA->point1, segmentA->point2, 0.0f};
	return s2CollideCapsuleAndCircle(&capsuleA, xfA, circleB, xfB);
}

s2Manifold s2CollideSegmentAndPolygon(const s2Segment* segmentA, s2Transform xfA, const s2Polygon* polygonB, s2Transform xfB,
									  s2DistanceCache* cache)
{
	s2Polygon polygonA = s2MakeCapsule(segmentA->point1, segmentA->point2, 0.0f);
	return s2CollidePolygons(&polygonA, xfA, polygonB, xfB, cache);
}
