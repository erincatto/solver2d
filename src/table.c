// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "table.h"

#include "allocate.h"
#include "core.h"

#include "solver2d/types.h"

#include <stdbool.h>
#include <string.h>

#if _DEBUG
int32_t g_probeCount;
#endif

static inline bool s2IsPowerOf2(uint32_t x)
{
	return (x & (x - 1)) == 0;
}

static inline uint32_t s2RoundUpPowerOf2(uint32_t x)
{
	if (s2IsPowerOf2(x))
	{
		return x;
	}

	// Hacker's Delight p48
	// Can also use ctz, but perf is not needed in this use case.
	x -= 1;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;

	return x + 1;
}

s2Set s2CreateSet(int32_t capacity)
{
	s2Set set = {0};

	// Capacity must be a power of 2
	if (capacity > 16)
	{
		set.capacity = s2RoundUpPowerOf2(capacity);
	}
	else
	{
		set.capacity = 16;
	}

	set.count = 0;
	set.items = s2Alloc(capacity * sizeof(s2SetItem));
	memset(set.items, 0, capacity * sizeof(s2SetItem));

	return set;
}

void s2DestroySet(s2Set* set)
{
	s2Free(set->items, set->capacity * sizeof(s2SetItem));
	set->items = NULL;
	set->count = 0;
	set->capacity = 0;
}

void s2ClearSet(s2Set* set)
{
	set->count = 0;
	memset(set->items, 0, set->capacity * sizeof(s2SetItem));
}

// I need a good hash because the keys are built from pairs of increasing integers.
// A simple hash like hash = (integer1 XOR integer2) has many collisions.
// https://lemire.me/blog/2018/08/15/fast-strongly-universal-64-bit-hashing-everywhere/
// https://preshing.com/20130107/this-hash-set-is-faster-than-a-judy-array/
// TODO_ERIN try: https://www.jandrewrogers.com/2019/02/12/fast-perfect-hashing/
static inline uint32_t s2KeyHash(uint64_t key)
{
	uint64_t h = key;
	h ^= h >> 33;
	h *= 0xff51afd7ed558ccdL;
	h ^= h >> 33;
	h *= 0xc4ceb9fe1a85ec53L;
	h ^= h >> 33;

	return (uint32_t)h;
}

#if _DEBUG
int32_t g_probeCount;
#endif

int32_t s2FindSlot(const s2Set* set, uint64_t key, uint32_t hash)
{
	uint32_t capacity = set->capacity;
	int32_t index = hash & (capacity - 1);
	const s2SetItem* items = set->items;
	while (items[index].hash != 0 && items[index].key != key)
	{
#if _DEBUG
		g_probeCount += 1;
#endif
		index = (index + 1) & (capacity - 1);
	}

	return index;
}

static void s2AddKeyHaveCapacity(s2Set* set, uint64_t key, uint32_t hash)
{
	int32_t index = s2FindSlot(set, key, hash);
	s2SetItem* items = set->items;
	S2_ASSERT(items[index].hash == 0);

	items[index].key = key;
	items[index].hash = hash;
	set->count += 1;
}

static void s2GrowTable(s2Set* set)
{
	uint32_t oldCount = set->count;
	S2_MAYBE_UNUSED(oldCount);

	uint32_t oldCapacity = set->capacity;
	s2SetItem* oldItems = set->items;

	set->count = 0;
	// Capacity must be a power of 2
	set->capacity = 2 * oldCapacity;
	set->items = s2Alloc(set->capacity * sizeof(s2SetItem));
	memset(set->items, 0, set->capacity * sizeof(s2SetItem));

	// Transfer items into new array
	for (uint32_t i = 0; i < oldCapacity; ++i)
	{
		s2SetItem* item = oldItems + i;
		if (item->hash == 0)
		{
			// this item was empty
			continue;
		}

		s2AddKeyHaveCapacity(set, item->key, item->hash);
	}

	S2_ASSERT(set->count == oldCount);

	s2Free(oldItems, oldCapacity * sizeof(s2SetItem));
}

bool s2ContainsKey(const s2Set* set, uint64_t key)
{
	uint32_t hash = s2KeyHash(key);
	int32_t index = s2FindSlot(set, key, hash);
	return set->items[index].key == key;
}

bool s2AddKey(s2Set* set, uint64_t key)
{
	uint32_t hash = s2KeyHash(key);
	int32_t index = s2FindSlot(set, key, hash);
	if (set->items[index].hash != 0)
	{
		// Already in set
		S2_ASSERT(set->items[index].hash == hash && set->items[index].key == key);
		return true;
	}

	if (2 * set->count >= set->capacity)
	{
		s2GrowTable(set);
	}

	s2AddKeyHaveCapacity(set, key, hash);
	return false;
}

// See https://en.wikipedia.org/wiki/Open_addressing
bool s2RemoveKey(s2Set* set, uint64_t key)
{
	uint32_t hash = s2KeyHash(key);
	int32_t i = s2FindSlot(set, key, hash);
	s2SetItem* items = set->items;
	if (items[i].hash == 0)
	{
		// Not in set
		return false;
	}

	// Mark item i as unoccupied
	items[i].key = 0;
	items[i].hash = 0;

	S2_ASSERT(set->count > 0);
	set->count -= 1;

	// Attempt to fill item i
	int32_t j = i;
	uint32_t capacity = set->capacity;
	for (;;)
	{
		j = (j + 1) & (capacity - 1);
		if (items[j].hash == 0)
		{
			break;
		}

		// k is the first item for the hash of j
		int32_t k = items[j].hash & (capacity - 1);

		// determine if k lies cyclically in (i,j]
		// i <= j: | i..k..j |
		// i > j: |.k..j  i....| or |....j     i..k.|
		if (i <= j)
		{
			if (i < k && k <= j)
			{
				continue;
			}
		}
		else
		{
			if (i < k || k <= j)
			{
				continue;
			}
		}

		// Move j into i
		items[i] = items[j];

		// Mark item j as unoccupied
		items[j].key = 0;
		items[j].hash = 0;

		i = j;
	}

	return true;
}
