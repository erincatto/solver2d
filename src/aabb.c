// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "solver2d/aabb.h"
#include "solver2d/constants.h"
#include "solver2d/math.h"

#include <float.h>

bool s2AABB_IsValid(s2Box a)
{
	s2Vec2 d = s2Sub(a.upperBound, a.lowerBound);
	bool valid = d.x >= 0.0f && d.y >= 0.0f;
	valid = valid && s2IsValidVec2(a.lowerBound) && s2IsValidVec2(a.upperBound);
	return valid;
}

// From Real-time Collision Detection, p179.
s2RayCastOutput s2AABB_RayCast(s2Box a, s2Vec2 p1, s2Vec2 p2)
{
	// Radius not handled
	s2RayCastOutput output = {0};

	float tmin = -FLT_MAX;
	float tmax = FLT_MAX;

	s2Vec2 p = p1;
	s2Vec2 d = s2Sub(p2, p1);
	s2Vec2 absD = s2Abs(d);

	s2Vec2 normal = s2Vec2_zero;

	// x-coordinate
	if (absD.x < FLT_EPSILON)
	{
		// parallel
		if (p.x < a.lowerBound.x || a.upperBound.x < p.x)
		{
			return output;
		}
	}
	else
	{
		float inv_d = 1.0f / d.x;
		float t1 = (a.lowerBound.x - p.x) * inv_d;
		float t2 = (a.upperBound.x - p.x) * inv_d;

		// Sign of the normal vector.
		float s = -1.0f;

		if (t1 > t2)
		{
			float tmp = t1;
			t1 = t2;
			t2 = tmp;
			s = 1.0f;
		}

		// Push the min up
		if (t1 > tmin)
		{
			normal.y = 0.0f;
			normal.x = s;
			tmin = t1;
		}

		// Pull the max down
		tmax = S2_MIN(tmax, t2);

		if (tmin > tmax)
		{
			return output;
		}
	}

	// y-coordinate
	if (absD.y < FLT_EPSILON)
	{
		// parallel
		if (p.y < a.lowerBound.y || a.upperBound.y < p.y)
		{
			return output;
		}
	}
	else
	{
		float inv_d = 1.0f / d.y;
		float t1 = (a.lowerBound.y - p.y) * inv_d;
		float t2 = (a.upperBound.y - p.y) * inv_d;

		// Sign of the normal vector.
		float s = -1.0f;

		if (t1 > t2)
		{
			float tmp = t1;
			t1 = t2;
			t2 = tmp;
			s = 1.0f;
		}

		// Push the min up
		if (t1 > tmin)
		{
			normal.x = 0.0f;
			normal.y = s;
			tmin = t1;
		}

		// Pull the max down
		tmax = S2_MIN(tmax, t2);

		if (tmin > tmax)
		{
			return output;
		}
	}
	
	// Does the ray start inside the box?
	// Does the ray intersect beyond the max fraction?
	if (tmin < 0.0f || 1.0f < tmin)
	{
		return output;
	}

	// Intersection.
	output.fraction = tmin;
	output.normal = normal;
	output.point = s2Lerp(p1, p2, tmin);
	output.hit = true;
	return output;
}
