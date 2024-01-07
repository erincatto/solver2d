// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "array.h"

#include <stdbool.h>
#include <stdint.h>

// A pooled object has a stable index, but not a stable pointer.
 
// Any pooled struct must have this as the first member.
typedef struct s2Object
{
	int32_t index;
	int32_t next;	
	uint16_t revision;
} s2Object;

typedef struct s2Pool
{
	char* memory;
	int32_t objectSize;
	int32_t capacity;
	int32_t count;
	int32_t freeList;
} s2Pool;

s2Pool s2CreatePool(int32_t objectSize, int32_t capacity);
void s2DestroyPool(s2Pool* pool);

s2Object* s2AllocObject(s2Pool* pool);
void s2FreeObject(s2Pool* pool, s2Object* object);

void s2GrowPool(s2Pool* pool, int32_t capacity);

static inline bool s2ObjectValid(const s2Object* object)
{
	// this means the object is not on the free list
	return object->index == object->next;
}

static inline bool s2IsFree(const s2Object* object)
{
	return object->index != object->next;
}
