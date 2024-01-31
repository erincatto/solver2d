// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "core.h"

#include "solver2d/constants.h"
#include "solver2d/math.h"

#include <float.h>

bool s2IsValid(float a)
{
	if (isnan(a))
	{
		return false;
	}

	if (isinf(a))
	{
		return false;
	}

	return true;
}

bool s2IsValidVec2(s2Vec2 v)
{
	if (isnan(v.x) || isnan(v.y))
	{
		return false;
	}

	if (isinf(v.x) || isinf(v.y))
	{
		return false;
	}

	return true;
}

s2Vec2 s2Normalize(s2Vec2 v)
{
	float length = s2Length(v);
	if (length < 0.001f * FLT_EPSILON)
	{
		return s2Vec2_zero;
	}

	float invLength = 1.0f / length;
	s2Vec2 n = {invLength * v.x, invLength * v.y};
	return n;
}

s2Vec2 s2NormalizeChecked(s2Vec2 v)
{
	float length = s2Length(v);
	if (length < FLT_EPSILON)
	{
		S2_ASSERT(false);
		return s2Vec2_zero;
	}

	float invLength = 1.0f / length;
	s2Vec2 n = {invLength * v.x, invLength * v.y};
	return n;
}

s2Vec2 s2GetLengthAndNormalize(float* length, s2Vec2 v)
{
	*length = s2Length(v);
	if (*length < FLT_EPSILON)
	{
		return s2Vec2_zero;
	}

	float invLength = 1.0f / *length;
	s2Vec2 n = {invLength * v.x, invLength * v.y};
	return n;
}
