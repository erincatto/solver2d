// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "types.h"

typedef struct s2Statistics
{
	int32_t bodyCount;
	int32_t contactCount;
	int32_t jointCount;
	int32_t proxyCount;
	int32_t treeHeight;
	int32_t stackCapacity;
	int32_t stackUsed;
} s2Statistics;
