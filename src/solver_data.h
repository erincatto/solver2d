// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver2d/types.h"

typedef struct s2StepContext
{
	float dt;
	float inv_dt;
	int32_t velocityIterations;
	int32_t positionIterations;
	float restitutionThreshold;
	struct s2Body* bodies;
	int32_t bodyCapacity;
} s2StepContext;
