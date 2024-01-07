
// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "solver2d/distance.h"

#include "core.h"

#include "solver2d/constants.h"
#include "solver2d/math.h"
#include "solver2d/timer.h"

#include <float.h>

/// Follows Ericson 5.1.9 Closest Points of Two Line Segments
s2SegmentDistanceResult s2SegmentDistance(s2Vec2 p1, s2Vec2 q1, s2Vec2 p2, s2Vec2 q2)
{
	s2SegmentDistanceResult result = {0};

	s2Vec2 d1 = s2Sub(q1, p1);
	s2Vec2 d2 = s2Sub(q2, p2);
	s2Vec2 r = s2Sub(p1, p2);
	float dd1 = s2Dot(d1, d1);
	float dd2 = s2Dot(d2, d2);
	float rd2 = s2Dot(r, d2);
	float rd1 = s2Dot(r, d1);

	const float epsSqr = FLT_EPSILON * FLT_EPSILON;

	if (dd1 < epsSqr || dd2 < epsSqr)
	{
		// Handle all degeneracies
		if (dd1 >= epsSqr)
		{
			// Segment 2 is degenerate
			result.fraction1 = S2_CLAMP(-rd1 / dd1, 0.0f, 1.0f);
			result.fraction2 = 0.0f;
		}
		else if (dd2 >= epsSqr)
		{
			// Segment 1 is degenerate
			result.fraction1 = 0.0f;
			result.fraction2 = S2_CLAMP(rd2 / dd2, 0.0f, 1.0f);
		}
		else
		{
			result.fraction1 = 0.0f;
			result.fraction2 = 0.0f;
		}
	}
	else
	{
		// Non-degenerate segments
		float d12 = s2Dot(d1, d2);

		float denom = dd1 * dd2 - d12 * d12;

		// Fraction on segment 1
		float f1 = 0.0f;
		if (denom != 0.0f)
		{
			// not parallel
			f1 = S2_CLAMP((d12 * rd2 - rd1 * dd2) / denom, 0.0f, 1.0f);
		}

		// Compute point on segment 2 closest to p1 + f1 * d1
		float f2 = (d12 * f1 + rd2) / dd2;

		// Clamping of segment 2 requires a do over on segment 1
		if (f2 < 0.0f)
		{
			f2 = 0.0f;
			f1 = S2_CLAMP(-rd1 / dd1, 0.0f, 1.0f);
		}
		else if (f2 > 1.0f)
		{
			f2 = 1.0f;
			f1 = S2_CLAMP((d12 - rd1) / dd1, 0.0f, 1.0f);
		}

		result.fraction1 = f1;
		result.fraction2 = f2;
	}

	result.closest1 = s2MulAdd(p1, result.fraction1, d1);
	result.closest2 = s2MulAdd(p2, result.fraction2, d2);
	result.distanceSquared = s2DistanceSquared(result.closest1, result.closest2);
	return result;
}

// GJK using Voronoi regions (Christer Ericson) and Barycentric coordinates.

s2DistanceProxy s2MakeProxy(const s2Vec2* vertices, int32_t count, float radius)
{
	count = S2_MIN(count, s2_maxPolygonVertices);
	s2DistanceProxy proxy;
	for (int32_t i = 0; i < count; ++i)
	{
		proxy.vertices[i] = vertices[i];
	}
	proxy.count = count;
	proxy.radius = radius;
	return proxy;
}

static s2Vec2 s2Weight2(float a1, s2Vec2 w1, float a2, s2Vec2 w2)
{
	return (s2Vec2){a1 * w1.x + a2 * w2.x, a1 * w1.y + a2 * w2.y};
}

static s2Vec2 s2Weight3(float a1, s2Vec2 w1, float a2, s2Vec2 w2, float a3, s2Vec2 w3)
{
	return (s2Vec2){a1 * w1.x + a2 * w2.x + a3 * w3.x, a1 * w1.y + a2 * w2.y + a3 * w3.y};
}

static int32_t s2FindSupport(const s2DistanceProxy* proxy, s2Vec2 direction)
{
	int32_t bestIndex = 0;
	float bestValue = s2Dot(proxy->vertices[0], direction);
	for (int32_t i = 1; i < proxy->count; ++i)
	{
		float value = s2Dot(proxy->vertices[i], direction);
		if (value > bestValue)
		{
			bestIndex = i;
			bestValue = value;
		}
	}

	return bestIndex;
}

typedef struct s2SimplexVertex
{
	s2Vec2 wA;		// support point in proxyA
	s2Vec2 wB;		// support point in proxyB
	s2Vec2 w;		// wB - wA
	float a;		// barycentric coordinate for closest point
	int32_t indexA; // wA index
	int32_t indexB; // wB index
} s2SimplexVertex;

typedef struct s2Simplex
{
	s2SimplexVertex v1, v2, v3;
	int32_t count;
} s2Simplex;

static float s2Simplex_Metric(const s2Simplex* s)
{
	switch (s->count)
	{
		case 0:
			S2_ASSERT(false);
			return 0.0f;

		case 1:
			return 0.0f;

		case 2:
			return s2Distance(s->v1.w, s->v2.w);

		case 3:
			return s2Cross(s2Sub(s->v2.w, s->v1.w), s2Sub(s->v3.w, s->v1.w));

		default:
			S2_ASSERT(false);
			return 0.0f;
	}
}

static s2Simplex s2MakeSimplexFromCache(const s2DistanceCache* cache, const s2DistanceProxy* proxyA, s2Transform transformA,
										const s2DistanceProxy* proxyB, s2Transform transformB)
{
	S2_ASSERT(cache->count <= 3);
	s2Simplex s;

	// Copy data from cache.
	s.count = cache->count;

	s2SimplexVertex* vertices[] = {&s.v1, &s.v2, &s.v3};
	for (int32_t i = 0; i < s.count; ++i)
	{
		s2SimplexVertex* v = vertices[i];
		v->indexA = cache->indexA[i];
		v->indexB = cache->indexB[i];
		s2Vec2 wALocal = proxyA->vertices[v->indexA];
		s2Vec2 wBLocal = proxyB->vertices[v->indexB];
		v->wA = s2TransformPoint(transformA, wALocal);
		v->wB = s2TransformPoint(transformB, wBLocal);
		v->w = s2Sub(v->wB, v->wA);

		// invalid
		v->a = -1.0f;
	}

	// If the cache is empty or invalid ...
	if (s.count == 0)
	{
		s2SimplexVertex* v = vertices[0];
		v->indexA = 0;
		v->indexB = 0;
		s2Vec2 wALocal = proxyA->vertices[0];
		s2Vec2 wBLocal = proxyB->vertices[0];
		v->wA = s2TransformPoint(transformA, wALocal);
		v->wB = s2TransformPoint(transformB, wBLocal);
		v->w = s2Sub(v->wB, v->wA);
		v->a = 1.0f;
		s.count = 1;
	}

	return s;
}

static void s2MakeSimplexCache(s2DistanceCache* cache, const s2Simplex* simplex)
{
	cache->metric = s2Simplex_Metric(simplex);
	cache->count = (uint16_t)simplex->count;
	const s2SimplexVertex* vertices[] = {&simplex->v1, &simplex->v2, &simplex->v3};
	for (int32_t i = 0; i < simplex->count; ++i)
	{
		cache->indexA[i] = (uint8_t)vertices[i]->indexA;
		cache->indexB[i] = (uint8_t)vertices[i]->indexB;
	}
}

s2Vec2 s2ComputeSimplexSearchDirection(const s2Simplex* simplex)
{
	switch (simplex->count)
	{
		case 1:
			return s2Neg(simplex->v1.w);

		case 2:
		{
			s2Vec2 e12 = s2Sub(simplex->v2.w, simplex->v1.w);
			float sgn = s2Cross(e12, s2Neg(simplex->v1.w));
			if (sgn > 0.0f)
			{
				// Origin is left of e12.
				return s2CrossSV(1.0f, e12);
			}
			else
			{
				// Origin is right of e12.
				return s2CrossVS(e12, 1.0f);
			}
		}

		default:
			S2_ASSERT(false);
			return s2Vec2_zero;
	}
}

s2Vec2 s2ComputeSimplexClosestPoint(const s2Simplex* s)
{
	switch (s->count)
	{
		case 0:
			S2_ASSERT(false);
			return s2Vec2_zero;

		case 1:
			return s->v1.w;

		case 2:
			return s2Weight2(s->v1.a, s->v1.w, s->v2.a, s->v2.w);

		case 3:
			return s2Vec2_zero;

		default:
			S2_ASSERT(false);
			return s2Vec2_zero;
	}
}

void s2ComputeSimplexWitnessPoints(s2Vec2* a, s2Vec2* b, const s2Simplex* s)
{
	switch (s->count)
	{
		case 0:
			S2_ASSERT(false);
			break;

		case 1:
			*a = s->v1.wA;
			*b = s->v1.wB;
			break;

		case 2:
			*a = s2Weight2(s->v1.a, s->v1.wA, s->v2.a, s->v2.wA);
			*b = s2Weight2(s->v1.a, s->v1.wB, s->v2.a, s->v2.wB);
			break;

		case 3:
			*a = s2Weight3(s->v1.a, s->v1.wA, s->v2.a, s->v2.wA, s->v3.a, s->v3.wA);
			// TODO_ERIN why are these not equal?
			//*b = s2Weight3(s->v1.a, s->v1.wB, s->v2.a, s->v2.wB, s->v3.a, s->v3.wB);
			*b = *a;
			break;

		default:
			S2_ASSERT(false);
			break;
	}
}

// Solve a line segment using barycentric coordinates.
//
// p = a1 * w1 + a2 * w2
// a1 + a2 = 1
//
// The vector from the origin to the closest point on the line is
// perpendicular to the line.
// e12 = w2 - w1
// dot(p, e) = 0
// a1 * dot(w1, e) + a2 * dot(w2, e) = 0
//
// 2-by-2 linear system
// [1      1     ][a1] = [1]
// [w1.e12 w2.e12][a2] = [0]
//
// Define
// d12_1 =  dot(w2, e12)
// d12_2 = -dot(w1, e12)
// d12 = d12_1 + d12_2
//
// Solution
// a1 = d12_1 / d12
// a2 = d12_2 / d12
void s2SolveSimplex2(s2Simplex* s)
{
	s2Vec2 w1 = s->v1.w;
	s2Vec2 w2 = s->v2.w;
	s2Vec2 e12 = s2Sub(w2, w1);

	// w1 region
	float d12_2 = -s2Dot(w1, e12);
	if (d12_2 <= 0.0f)
	{
		// a2 <= 0, so we clamp it to 0
		s->v1.a = 1.0f;
		s->count = 1;
		return;
	}

	// w2 region
	float d12_1 = s2Dot(w2, e12);
	if (d12_1 <= 0.0f)
	{
		// a1 <= 0, so we clamp it to 0
		s->v2.a = 1.0f;
		s->count = 1;
		s->v1 = s->v2;
		return;
	}

	// Must be in e12 region.
	float inv_d12 = 1.0f / (d12_1 + d12_2);
	s->v1.a = d12_1 * inv_d12;
	s->v2.a = d12_2 * inv_d12;
	s->count = 2;
}

void s2SolveSimplex3(s2Simplex* s)
{
	s2Vec2 w1 = s->v1.w;
	s2Vec2 w2 = s->v2.w;
	s2Vec2 w3 = s->v3.w;

	// Edge12
	// [1      1     ][a1] = [1]
	// [w1.e12 w2.e12][a2] = [0]
	// a3 = 0
	s2Vec2 e12 = s2Sub(w2, w1);
	float w1e12 = s2Dot(w1, e12);
	float w2e12 = s2Dot(w2, e12);
	float d12_1 = w2e12;
	float d12_2 = -w1e12;

	// Edge13
	// [1      1     ][a1] = [1]
	// [w1.e13 w3.e13][a3] = [0]
	// a2 = 0
	s2Vec2 e13 = s2Sub(w3, w1);
	float w1e13 = s2Dot(w1, e13);
	float w3e13 = s2Dot(w3, e13);
	float d13_1 = w3e13;
	float d13_2 = -w1e13;

	// Edge23
	// [1      1     ][a2] = [1]
	// [w2.e23 w3.e23][a3] = [0]
	// a1 = 0
	s2Vec2 e23 = s2Sub(w3, w2);
	float w2e23 = s2Dot(w2, e23);
	float w3e23 = s2Dot(w3, e23);
	float d23_1 = w3e23;
	float d23_2 = -w2e23;

	// Triangle123
	float n123 = s2Cross(e12, e13);

	float d123_1 = n123 * s2Cross(w2, w3);
	float d123_2 = n123 * s2Cross(w3, w1);
	float d123_3 = n123 * s2Cross(w1, w2);

	// w1 region
	if (d12_2 <= 0.0f && d13_2 <= 0.0f)
	{
		s->v1.a = 1.0f;
		s->count = 1;
		return;
	}

	// e12
	if (d12_1 > 0.0f && d12_2 > 0.0f && d123_3 <= 0.0f)
	{
		float inv_d12 = 1.0f / (d12_1 + d12_2);
		s->v1.a = d12_1 * inv_d12;
		s->v2.a = d12_2 * inv_d12;
		s->count = 2;
		return;
	}

	// e13
	if (d13_1 > 0.0f && d13_2 > 0.0f && d123_2 <= 0.0f)
	{
		float inv_d13 = 1.0f / (d13_1 + d13_2);
		s->v1.a = d13_1 * inv_d13;
		s->v3.a = d13_2 * inv_d13;
		s->count = 2;
		s->v2 = s->v3;
		return;
	}

	// w2 region
	if (d12_1 <= 0.0f && d23_2 <= 0.0f)
	{
		s->v2.a = 1.0f;
		s->count = 1;
		s->v1 = s->v2;
		return;
	}

	// w3 region
	if (d13_1 <= 0.0f && d23_1 <= 0.0f)
	{
		s->v3.a = 1.0f;
		s->count = 1;
		s->v1 = s->v3;
		return;
	}

	// e23
	if (d23_1 > 0.0f && d23_2 > 0.0f && d123_1 <= 0.0f)
	{
		float inv_d23 = 1.0f / (d23_1 + d23_2);
		s->v2.a = d23_1 * inv_d23;
		s->v3.a = d23_2 * inv_d23;
		s->count = 2;
		s->v1 = s->v3;
		return;
	}

	// Must be in triangle123
	float inv_d123 = 1.0f / (d123_1 + d123_2 + d123_3);
	s->v1.a = d123_1 * inv_d123;
	s->v2.a = d123_2 * inv_d123;
	s->v3.a = d123_3 * inv_d123;
	s->count = 3;
}

#define S2_GJK_DEBUG 0

// Warning: writing to these globals significantly slows multi-threading performance
#if S2_GJK_DEBUG
int32_t s2_gjkCalls;
int32_t s2_gjkIters;
int32_t s2_gjkMaxIters;
#endif

s2DistanceOutput s2ShapeDistance(s2DistanceCache* cache, const s2DistanceInput* input)
{
#if S2_GJK_DEBUG
	++s2_gjkCalls;
#endif

	s2DistanceOutput output = {0};

	const s2DistanceProxy* proxyA = &input->proxyA;
	const s2DistanceProxy* proxyB = &input->proxyB;

	s2Transform transformA = input->transformA;
	s2Transform transformB = input->transformB;

	// Initialize the simplex.
	s2Simplex simplex = s2MakeSimplexFromCache(cache, proxyA, transformA, proxyB, transformB);

	// Get simplex vertices as an array.
	s2SimplexVertex* vertices[] = {&simplex.v1, &simplex.v2, &simplex.v3};
	const int32_t k_maxIters = 20;

	// These store the vertices of the last simplex so that we
	// can check for duplicates and prevent cycling.
	int32_t saveA[3], saveB[3];
	int32_t saveCount = 0;

	// Main iteration loop.
	int32_t iter = 0;
	while (iter < k_maxIters)
	{
		// Copy simplex so we can identify duplicates.
		saveCount = simplex.count;
		for (int32_t i = 0; i < saveCount; ++i)
		{
			saveA[i] = vertices[i]->indexA;
			saveB[i] = vertices[i]->indexB;
		}

		switch (simplex.count)
		{
			case 1:
				break;

			case 2:
				s2SolveSimplex2(&simplex);
				break;

			case 3:
				s2SolveSimplex3(&simplex);
				break;

			default:
				S2_ASSERT(false);
		}

		// If we have 3 points, then the origin is in the corresponding triangle.
		if (simplex.count == 3)
		{
			break;
		}

		// Get search direction.
		s2Vec2 d = s2ComputeSimplexSearchDirection(&simplex);

		// Ensure the search direction is numerically fit.
		if (s2Dot(d, d) < FLT_EPSILON * FLT_EPSILON)
		{
			// The origin is probably contained by a line segment
			// or triangle. Thus the shapes are overlapped.

			// We can't return zero here even though there may be overlap.
			// In case the simplex is a point, segment, or triangle it is difficult
			// to determine if the origin is contained in the CSO or very close to it.
			break;
		}

		// Compute a tentative new simplex vertex using support points.
		s2SimplexVertex* vertex = vertices[simplex.count];
		vertex->indexA = s2FindSupport(proxyA, s2InvRotateVector(transformA.q, s2Neg(d)));
		vertex->wA = s2TransformPoint(transformA, proxyA->vertices[vertex->indexA]);
		vertex->indexB = s2FindSupport(proxyB, s2InvRotateVector(transformB.q, d));
		vertex->wB = s2TransformPoint(transformB, proxyB->vertices[vertex->indexB]);
		vertex->w = s2Sub(vertex->wB, vertex->wA);

		// Iteration count is equated to the number of support point calls.
		++iter;

#if S2_GJK_DEBUG
		++s2_gjkIters;
#endif

		// Check for duplicate support points. This is the main termination criteria.
		bool duplicate = false;
		for (int32_t i = 0; i < saveCount; ++i)
		{
			if (vertex->indexA == saveA[i] && vertex->indexB == saveB[i])
			{
				duplicate = true;
				break;
			}
		}

		// If we found a duplicate support point we must exit to avoid cycling.
		if (duplicate)
		{
			break;
		}

		// New vertex is ok and needed.
		++simplex.count;
	}

#if S2_GJK_DEBUG
	s2_gjkMaxIters = S2_MAX(s2_gjkMaxIters, iter);
#endif

	// Prepare output
	s2ComputeSimplexWitnessPoints(&output.pointA, &output.pointB, &simplex);
	output.distance = s2Distance(output.pointA, output.pointB);
	output.iterations = iter;

	// Cache the simplex
	s2MakeSimplexCache(cache, &simplex);

	// Apply radii if requested
	if (input->useRadii)
	{
		if (output.distance < FLT_EPSILON)
		{
			// Shapes are too close to safely compute normal
			s2Vec2 p = (s2Vec2){0.5f * (output.pointA.x + output.pointB.x), 0.5f * (output.pointA.y + output.pointB.y)};
			output.pointA = p;
			output.pointB = p;
			output.distance = 0.0f;
		}
		else
		{
			// Keep closest points on perimeter even if overlapped, this way
			// the points move smoothly.
			float rA = proxyA->radius;
			float rB = proxyB->radius;
			output.distance = S2_MAX(0.0f, output.distance - rA - rB);
			s2Vec2 normal = s2Normalize(s2Sub(output.pointB, output.pointA));
			s2Vec2 offsetA = (s2Vec2){rA * normal.x, rA * normal.y};
			s2Vec2 offsetB = (s2Vec2){rB * normal.x, rB * normal.y};
			output.pointA = s2Add(output.pointA, offsetA);
			output.pointB = s2Sub(output.pointB, offsetB);
		}
	}

	return output;
}
