// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "pool.h"

#include "allocate.h"
#include "core.h"
#include "solver2d/types.h"

#include <string.h>

void s2ValidatePool(const s2Pool* pool)
{
#if S2_VALIDATE
	int32_t freeCount = 0;
	int32_t freeIndex = pool->freeList;
	int32_t objectSize = pool->objectSize;
	int32_t capacity = pool->capacity;

	while (freeIndex != S2_NULL_INDEX)
	{
		S2_ASSERT(0 <= freeIndex && freeIndex < pool->capacity);
		s2Object* object = (s2Object*)(pool->memory + freeIndex * objectSize);
		S2_ASSERT(object->next != object->index);
		freeIndex = object->next;
		freeCount += 1;
	}

	S2_ASSERT(freeCount + pool->count == capacity);
#else
	S2_MAYBE_UNUSED(pool);
#endif
}

s2Pool s2CreatePool(int32_t objectSize, int32_t capacity)
{
	S2_ASSERT(objectSize >= (int32_t)sizeof(s2Object));

	s2Pool pool;
	pool.objectSize = objectSize;
	pool.capacity = capacity > 1 ? capacity : 1;
	pool.count = 0;
	pool.memory = (char*)s2Alloc(pool.capacity * objectSize);

	pool.freeList = 0;
	for (int32_t i = 0; i < pool.capacity - 1; ++i)
	{
		s2Object* object = (s2Object*)(pool.memory + i * objectSize);
		object->index = i;
		object->next = i + 1;
		object->revision = 0;
	}

	s2Object* object = (s2Object*)(pool.memory + (pool.capacity - 1) * objectSize);
	object->index = pool.capacity - 1;
	object->next = S2_NULL_INDEX;
	object->revision = 0;

	return pool;
}

void s2DestroyPool(s2Pool* pool)
{
	s2Free(pool->memory, pool->capacity * pool->objectSize);
	pool->memory = NULL;
	pool->capacity = 0;
	pool->count = 0;
	pool->freeList = S2_NULL_INDEX;
	pool->objectSize = 0;
}

void s2GrowPool(s2Pool* pool, int32_t capacity)
{
	int32_t oldCapacity = pool->capacity;
	if (oldCapacity >= capacity)
	{
		return;
	}

	int32_t newCapacity = capacity > 2 ? capacity : 2;
	pool->capacity = newCapacity;
	char* newMemory = (char*)s2Alloc(pool->capacity * pool->objectSize);
	memcpy(newMemory, pool->memory, oldCapacity * pool->objectSize);
	s2Free(pool->memory, oldCapacity * pool->objectSize);
	pool->memory = newMemory;

	int32_t oldFreeList = pool->freeList;
	pool->freeList = oldCapacity;
	for (int32_t i = oldCapacity; i < newCapacity - 1; ++i)
	{
		s2Object* object = (s2Object*)(pool->memory + i * pool->objectSize);
		object->index = i;
		object->next = i + 1;
		object->revision = 0;
	}

	// Tail of free list
	s2Object* object = (s2Object*)(pool->memory + (newCapacity - 1) * pool->objectSize);
	object->index = newCapacity - 1;
	object->next = oldFreeList;
	object->revision = 0;

#if S2_VALIDATE
	s2ValidatePool(pool);
#endif
}

s2Object* s2AllocObject(s2Pool* pool)
{
	s2Object* newObject = NULL;
	if (pool->freeList != S2_NULL_INDEX)
	{
		newObject = (s2Object*)(pool->memory + pool->freeList * pool->objectSize);
		newObject->index = pool->freeList;
		newObject->revision += 1;
		pool->freeList = newObject->next;
		newObject->next = newObject->index;
		pool->count += 1;
		return newObject;
	}
	else
	{
		int32_t oldCapacity = pool->capacity;
		int32_t newCapacity = oldCapacity + oldCapacity / 2;
		newCapacity = newCapacity > 2 ? newCapacity : 2;
		pool->capacity = newCapacity;
		char* newMemory = (char*)s2Alloc(pool->capacity * pool->objectSize);
		memcpy(newMemory, pool->memory, oldCapacity * pool->objectSize);
		s2Free(pool->memory, oldCapacity * pool->objectSize);
		pool->memory = newMemory;

		newObject = (s2Object*)(pool->memory + oldCapacity * pool->objectSize);
		newObject->index = oldCapacity;
		newObject->revision = 0;
		newObject->next = newObject->index;

		pool->freeList = oldCapacity + 1;
		for (int32_t i = oldCapacity + 1; i < newCapacity - 1; ++i)
		{
			s2Object* object = (s2Object*)(pool->memory + i * pool->objectSize);
			object->index = i;
			object->next = i + 1;
			object->revision = 0;
		}

		s2Object* object = (s2Object*)(pool->memory + (newCapacity - 1) * pool->objectSize);
		object->index = newCapacity - 1;
		object->next = S2_NULL_INDEX;
		object->revision = 0;

		pool->count += 1;

#if S2_VALIDATE
		s2ValidatePool(pool);
#endif

		return newObject;
	}
}

void s2FreeObject(s2Pool* pool, s2Object* object)
{
	S2_ASSERT(pool->memory <= (char*)object && (char*)object < pool->memory + pool->capacity * pool->objectSize);
	S2_ASSERT(object->index == object->next);
	S2_ASSERT(object->index < pool->capacity);

	object->next = pool->freeList;
	pool->freeList = object->index;
	pool->count -= 1;
}
