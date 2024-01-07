// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "array.h"

#include "allocate.h"
#include "core.h"

#include <string.h>

void* s2CreateArray(int32_t elementSize, int32_t capacity)
{
	void* result = (s2ArrayHeader*)s2Alloc(sizeof(s2ArrayHeader) + elementSize * capacity) + 1;
	s2Array(result).count = 0;
	s2Array(result).capacity = capacity;
	return result;
}

void s2DestroyArray(void* a, int32_t elementSize)
{
	int32_t capacity = s2Array(a).capacity;
	int32_t size = sizeof(s2ArrayHeader) + elementSize * capacity;
	s2Free(((s2ArrayHeader*)a) - 1, size);
}

void s2Array_Grow(void** a, int32_t elementSize)
{
	int32_t capacity = s2Array(*a).capacity;
	S2_ASSERT(capacity == s2Array(*a).count);

	// grow by 50%
	int32_t newCapacity = capacity + (capacity >> 1);
	newCapacity = newCapacity >= 2 ? newCapacity : 2;
	void* tmp = *a;
	*a = (s2ArrayHeader*)s2Alloc(sizeof(s2ArrayHeader) + elementSize * newCapacity) + 1;
	s2Array(*a).capacity = newCapacity;
	s2Array(*a).count = capacity;
	memcpy(*a, tmp, capacity * elementSize);
	s2DestroyArray(tmp, elementSize);
}
