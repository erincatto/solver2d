// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "constants.h"
#include "types.h"

typedef struct s2Hull
{
	s2Vec2 points[s2_maxPolygonVertices];
	int32_t count;
} s2Hull;

#ifdef __cplusplus
extern "C"
{
#endif

s2Hull s2ComputeHull(const s2Vec2* points, int32_t count);
bool s2ValidateHull(const s2Hull* hull);

#ifdef __cplusplus
}
#endif
