// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include <stdbool.h>
#include <stdint.h>

// todo add feature to look up elements by half-key, replace s2JointEdge and s2ContactEdge

#define S2_SHAPE_PAIR_KEY(K1, K2) K1 < K2 ? (uint64_t)K1 << 32 | (uint64_t)K2 : (uint64_t)K2 << 32 | (uint64_t)K1

typedef struct s2SetItem
{
	uint64_t key;
	uint32_t hash;
} s2SetItem;

typedef struct s2Set
{
	s2SetItem* items;
	uint32_t capacity;
	uint32_t count;
} s2Set;

s2Set s2CreateSet(int32_t capacity);
void s2DestroySet(s2Set* set);

void s2ClearSet(s2Set* set);

	// Returns true if key was already in set
bool s2AddKey(s2Set* set, uint64_t key);

// Returns true if the key was found
bool s2RemoveKey(s2Set* set, uint64_t key);

bool s2ContainsKey(const s2Set* set, uint64_t key);
