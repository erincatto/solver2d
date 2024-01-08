// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "solver2d/hull.h"

#include "core.h"

#include "solver2d/aabb.h"
#include "solver2d/math.h"

#include <float.h>

// quickhull recursion
static s2Hull s2RecurseHull(s2Vec2 p1, s2Vec2 p2, s2Vec2* ps, int32_t count)
{
	s2Hull hull;
	hull.count = 0;

	if (count == 0)
	{
		return hull;
	}

	// create an edge vector pointing from p1 to p2
	s2Vec2 e = s2Normalize(s2Sub(p2, p1));

	// discard points left of e and find point furthest to the right of e
	s2Vec2 rightPoints[s2_maxPolygonVertices];
	int32_t rightCount = 0;

	int32_t bestIndex = 0;
	float bestDistance = s2Cross(s2Sub(ps[bestIndex], p1), e);
	if (bestDistance > 0.0f)
	{
		rightPoints[rightCount++] = ps[bestIndex];
	}

	for (int32_t i = 1; i < count; ++i)
	{
		float distance = s2Cross(s2Sub(ps[i], p1), e);
		if (distance > bestDistance)
		{
			bestIndex = i;
			bestDistance = distance;
		}

		if (distance > 0.0f)
		{
			rightPoints[rightCount++] = ps[i];
		}
	}

	if (bestDistance < 2.0f * s2_linearSlop)
	{
		return hull;
	}

	s2Vec2 bestPoint = ps[bestIndex];

	// compute hull to the right of p1-bestPoint
	s2Hull hull1 = s2RecurseHull(p1, bestPoint, rightPoints, rightCount);

	// compute hull to the right of bestPoint-p2
	s2Hull hull2 = s2RecurseHull(bestPoint, p2, rightPoints, rightCount);

	// stitch together hulls
	for (int32_t i = 0; i < hull1.count; ++i)
	{
		hull.points[hull.count++] = hull1.points[i];
	}

	hull.points[hull.count++] = bestPoint;

	for (int32_t i = 0; i < hull2.count; ++i)
	{
		hull.points[hull.count++] = hull2.points[i];
	}

	S2_ASSERT(hull.count < s2_maxPolygonVertices);

	return hull;
}

// quickhull algorithm
// - merges vertices based on s2_linearSlop
// - removes collinear points using s2_linearSlop
// - returns an empty hull if it fails
s2Hull s2ComputeHull(const s2Vec2* points, int32_t count)
{
	s2Hull hull;
	hull.count = 0;

	if (count < 3 || count > s2_maxPolygonVertices)
	{
		// check your data
		return hull;
	}

	count = S2_MIN(count, s2_maxPolygonVertices);

	s2Box aabb = {{FLT_MAX, FLT_MAX}, {-FLT_MAX, -FLT_MAX}};

	// Perform aggressive point welding. First point always remains.
	// Also compute the bounding box for later.
	s2Vec2 ps[s2_maxPolygonVertices];
	int32_t n = 0;
	const float tolSqr = 16.0f * s2_linearSlop * s2_linearSlop;
	for (int32_t i = 0; i < count; ++i)
	{
		aabb.lowerBound = s2Min(aabb.lowerBound, points[i]);
		aabb.upperBound = s2Max(aabb.upperBound, points[i]);

		s2Vec2 vi = points[i];

		bool unique = true;
		for (int32_t j = 0; j < i; ++j)
		{
			s2Vec2 vj = points[j];

			float distSqr = s2DistanceSquared(vi, vj);
			if (distSqr < tolSqr)
			{
				unique = false;
				break;
			}
		}

		if (unique)
		{
			ps[n++] = vi;
		}
	}

	if (n < 3)
	{
		// all points very close together, check your data and check your scale
		return hull;
	}

	// Find an extreme point as the first point on the hull
	s2Vec2 c = s2AABB_Center(aabb);
	int32_t f1 = 0;
	float dsq1 = s2DistanceSquared(c, ps[f1]);
	for (int32_t i = 1; i < n; ++i)
	{
		float dsq = s2DistanceSquared(c, ps[i]);
		if (dsq > dsq1)
		{
			f1 = i;
			dsq1 = dsq;
		}
	}

	// remove p1 from working set
	s2Vec2 p1 = ps[f1];
	ps[f1] = ps[n - 1];
	n = n - 1;

	int32_t f2 = 0;
	float dsq2 = s2DistanceSquared(p1, ps[f2]);
	for (int32_t i = 1; i < n; ++i)
	{
		float dsq = s2DistanceSquared(p1, ps[i]);
		if (dsq > dsq2)
		{
			f2 = i;
			dsq2 = dsq;
		}
	}

	// remove p2 from working set
	s2Vec2 p2 = ps[f2];
	ps[f2] = ps[n - 1];
	n = n - 1;

	// split the points into points that are left and right of the line p1-p2.
	s2Vec2 rightPoints[s2_maxPolygonVertices - 2];
	int32_t rightCount = 0;

	s2Vec2 leftPoints[s2_maxPolygonVertices - 2];
	int32_t leftCount = 0;

	s2Vec2 e = s2Normalize(s2Sub(p2, p1));

	for (int32_t i = 0; i < n; ++i)
	{
		float d = s2Cross(s2Sub(ps[i], p1), e);

		// slop used here to skip points that are very close to the line p1-p2
		if (d >= 2.0f * s2_linearSlop)
		{
			rightPoints[rightCount++] = ps[i];
		}
		else if (d <= -2.0f * s2_linearSlop)
		{
			leftPoints[leftCount++] = ps[i];
		}
	}

	// compute hulls on right and left
	s2Hull hull1 = s2RecurseHull(p1, p2, rightPoints, rightCount);
	s2Hull hull2 = s2RecurseHull(p2, p1, leftPoints, leftCount);

	if (hull1.count == 0 && hull2.count == 0)
	{
		// all points collinear
		return hull;
	}

	// stitch hulls together, preserving CCW winding order
	hull.points[hull.count++] = p1;

	for (int32_t i = 0; i < hull1.count; ++i)
	{
		hull.points[hull.count++] = hull1.points[i];
	}

	hull.points[hull.count++] = p2;

	for (int32_t i = 0; i < hull2.count; ++i)
	{
		hull.points[hull.count++] = hull2.points[i];
	}

	S2_ASSERT(hull.count <= s2_maxPolygonVertices);

	// merge collinear
	bool searching = true;
	while (searching && hull.count > 2)
	{
		searching = false;

		for (int32_t i = 0; i < hull.count; ++i)
		{
			int32_t i1 = i;
			int32_t i2 = (i + 1) % hull.count;
			int32_t i3 = (i + 2) % hull.count;

			s2Vec2 s1 = hull.points[i1];
			s2Vec2 s2 = hull.points[i2];
			s2Vec2 s3 = hull.points[i3];

			// unit edge vector for s1-s3
			s2Vec2 r = s2Normalize(s2Sub(s3, s1));

			float distance = s2Cross(s2Sub(s2, s1), r);
			if (distance <= 2.0f * s2_linearSlop)
			{
				// remove midpoint from hull
				for (int32_t j = i2; j < hull.count - 1; ++j)
				{
					hull.points[j] = hull.points[j + 1];
				}
				hull.count -= 1;

				// continue searching for collinear points
				searching = true;

				break;
			}
		}
	}

	if (hull.count < 3)
	{
		// all points collinear, shouldn't be reached since this was validated above
		hull.count = 0;
	}

	return hull;
}

bool s2ValidateHull(const s2Hull* hull)
{
	if (hull->count < 3 || s2_maxPolygonVertices < hull->count)
	{
		return false;
	}

	// test that every point is behind every edge
	for (int32_t i = 0; i < hull->count; ++i)
	{
		// create an edge vector
		int32_t i1 = i;
		int32_t i2 = i < hull->count - 1 ? i1 + 1 : 0;
		s2Vec2 p = hull->points[i1];
		s2Vec2 e = s2Normalize(s2Sub(hull->points[i2], p));

		for (int32_t j = 0; j < hull->count; ++j)
		{
			// skip points that subtend the current edge
			if (j == i1 || j == i2)
			{
				continue;
			}

			float distance = s2Cross(s2Sub(hull->points[j], p), e);
			if (distance >= 0.0f)
			{
				return false;
			}
		}
	}

	// test for collinear points
	for (int32_t i = 0; i < hull->count; ++i)
	{
		int32_t i1 = i;
		int32_t i2 = (i + 1) % hull->count;
		int32_t i3 = (i + 2) % hull->count;

		s2Vec2 p1 = hull->points[i1];
		s2Vec2 p2 = hull->points[i2];
		s2Vec2 p3 = hull->points[i3];

		s2Vec2 e = s2Normalize(s2Sub(p3, p1));

		float distance = s2Cross(s2Sub(p2, p1), e);
		if (distance <= s2_linearSlop)
		{
			// p1-p2-p3 are collinear
			return false;
		}
	}

	return true;
}
