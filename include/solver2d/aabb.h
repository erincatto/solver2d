// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "constants.h"
#include "math.h"
#include "types.h"

// TODO_ERIN a lot of this could be internal

#ifdef __cplusplus
extern "C"
{
#endif

/// Verify that the bounds are sorted.
bool s2AABB_IsValid(s2Box a);

/// Ray cast an AABB
s2RayCastOutput s2AABB_RayCast(s2Box a, s2Vec2 p1, s2Vec2 p2);

#ifdef __cplusplus
}
#endif

/// Get the center of the AABB.
static inline s2Vec2 s2AABB_Center(s2Box a)
{
	s2Vec2 b = {0.5f * (a.lowerBound.x + a.upperBound.x), 0.5f * (a.lowerBound.y + a.upperBound.y)};
	return b;
}

/// Get the extents of the AABB (half-widths).
static inline s2Vec2 s2AABB_Extents(s2Box a)
{
	s2Vec2 b = {0.5f * (a.upperBound.x - a.lowerBound.x), 0.5f * (a.upperBound.y - a.lowerBound.y)};
	return b;
}

/// Get the perimeter length
static inline float s2AABB_Perimeter(s2Box a)
{
	float wx = a.upperBound.x - a.lowerBound.x;
	float wy = a.upperBound.y - a.lowerBound.y;
	return 2.0f * (wx + wy);
}

/// Union of two AABBs
static inline s2Box s2AABB_Union(s2Box a, s2Box b)
{
	s2Box c;
	c.lowerBound.x = S2_MIN(a.lowerBound.x, b.lowerBound.x);
	c.lowerBound.y = S2_MIN(a.lowerBound.y, b.lowerBound.y);
	c.upperBound.x = S2_MAX(a.upperBound.x, b.upperBound.x);
	c.upperBound.y = S2_MAX(a.upperBound.y, b.upperBound.y);
	return c;
}

/// Enlarge a to contain b
/// @return true if the AABB grew
static inline bool s2AABB_Enlarge(s2Box* a, s2Box b)
{
	bool changed = false;
	if (b.lowerBound.x < a->lowerBound.x)
	{
		a->lowerBound.x = b.lowerBound.x;
		changed = true;
	}

	if (b.lowerBound.y < a->lowerBound.y)
	{
		a->lowerBound.y = b.lowerBound.y;
		changed = true;
	}

	if (a->upperBound.x < b.upperBound.x)
	{
		a->upperBound.x = b.upperBound.x;
		changed = true;
	}

	if (a->upperBound.y < b.upperBound.y)
	{
		a->upperBound.y = b.upperBound.y;
		changed = true;
	}

	return changed;
}

/// Does a fully contain b
static inline bool s2AABB_Contains(s2Box a, s2Box b)
{
	bool s = true;
	s = s && a.lowerBound.x <= b.lowerBound.x;
	s = s && a.lowerBound.y <= b.lowerBound.y;
	s = s && b.upperBound.x <= a.upperBound.x;
	s = s && b.upperBound.y <= a.upperBound.y;
	return s;
}

static inline bool s2AABB_ContainsWithMargin(s2Box a, s2Box b, float margin)
{
	bool s = (a.lowerBound.x <= b.lowerBound.x - margin) & (a.lowerBound.y <= b.lowerBound.y - margin) &
			 (b.upperBound.x + margin <= a.upperBound.x) & (b.upperBound.y + margin <= a.upperBound.y);
	return s;
}

/// Do a and b overlap
static inline bool s2AABB_Overlaps(s2Box a, s2Box b)
{
	s2Vec2 d1 = {b.lowerBound.x - a.upperBound.x, b.lowerBound.y - a.upperBound.y};
	s2Vec2 d2 = {a.lowerBound.x - b.upperBound.x, a.lowerBound.y - b.upperBound.y};

	if (d1.x > 0.0f || d1.y > 0.0f)
		return false;

	if (d2.x > 0.0f || d2.y > 0.0f)
		return false;

	return true;
}
