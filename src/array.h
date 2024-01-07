// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "core.h"

#include <stdint.h>

typedef struct s2ArrayHeader
{
	int32_t count;
	int32_t capacity;
} s2ArrayHeader;

#define s2Array(a) ((s2ArrayHeader*)(a))[-1]

void* s2CreateArray(int32_t elementSize, int32_t capacity);
void s2DestroyArray(void* a, int32_t elementSize);
void s2Array_Grow(void** a, int32_t elementSize);

#define s2Array_Check(a, index) S2_ASSERT(0 <= index && index < s2Array(a).count)

#define s2Array_Clear(a) s2Array(a).count = 0

#define s2Array_Push(a, element) \
	if (s2Array(a).count == s2Array(a).capacity) s2Array_Grow((void**)&a, sizeof(element)); \
	S2_ASSERT(s2Array(a).count < s2Array(a).capacity); \
	a[s2Array(a).count++] = element

#define s2Array_RemoveSwap(a, index) \
	S2_ASSERT(0 <= index && index < s2Array(a).count); \
	a[index] = a[s2Array(a).count - 1]; \
	s2Array(a).count--

#define s2Array_Last(a) (a)[s2Array(a).count - 1];

#define s2Array_Pop(a) \
	S2_ASSERT(0 < s2Array(a).count); \
	s2Array(a).count--
