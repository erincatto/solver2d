// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "stack_allocator.h"

#include "allocate.h"
#include "array.h"
#include "core.h"

#include <stdbool.h>

typedef struct s2StackEntry
{
	char* data;
	const char* name;
	int32_t size;
	bool usedMalloc;
} s2StackEntry;

// This is a stack allocator used for fast per step allocations.
// You must nest allocate/free pairs. The code will S2_ASSERT
// if you try to interleave multiple allocate/free pairs.
// Unlike a scratch allocator, this lets me use the heap if the allocator
// space is insufficient.
typedef struct s2StackAllocator
{
	char* data;
	int32_t capacity;
	int32_t index;

	int32_t allocation;
	int32_t maxAllocation;

	s2StackEntry* entries;
} s2StackAllocator;

s2StackAllocator* s2CreateStackAllocator(int32_t capacity)
{
	S2_ASSERT(capacity >= 0);
	s2StackAllocator* allocator = s2Alloc(sizeof(s2StackAllocator));
	allocator->capacity = capacity;
	allocator->data = s2Alloc(capacity);
	allocator->allocation = 0;
	allocator->maxAllocation = 0;
	allocator->index = 0;
	allocator->entries = s2CreateArray(sizeof(s2StackEntry), 32);
	return allocator;
}

void s2DestroyStackAllocator(s2StackAllocator* allocator)
{
	s2DestroyArray(allocator->entries, sizeof(s2StackEntry));
	s2Free(allocator->data, allocator->capacity);
	s2Free(allocator, sizeof(s2StackAllocator));
}

void* s2AllocateStackItem(s2StackAllocator* alloc, int32_t size, const char* name)
{
	s2StackEntry entry;
	entry.size = size;
	entry.name = name;
	if (alloc->index + size > alloc->capacity)
	{
		// fall back to the heap (undesirable)
		entry.data = (char*)s2Alloc(size);
		entry.usedMalloc = true;
	}
	else
	{
		entry.data = alloc->data + alloc->index;
		entry.usedMalloc = false;
		alloc->index += size;
	}

	alloc->allocation += size;
	if (alloc->allocation > alloc->maxAllocation)
	{
		alloc->maxAllocation = alloc->allocation;
	}

	s2Array_Push(alloc->entries, entry);

	memset(entry.data, 0, size);

	return entry.data;
}

void s2FreeStackItem(s2StackAllocator* alloc, void* mem)
{
	int32_t entryCount = s2Array(alloc->entries).count;
	S2_ASSERT(entryCount > 0);
	s2StackEntry* entry = alloc->entries + (entryCount - 1);
	S2_ASSERT(mem == entry->data);
	if (entry->usedMalloc)
	{
		s2Free(mem, entry->size);
	}
	else
	{
		alloc->index -= entry->size;
	}
	alloc->allocation -= entry->size;
	s2Array_Pop(alloc->entries);
}

void s2GrowStack(s2StackAllocator* alloc)
{
	// Stack must not be in use
	S2_ASSERT(alloc->allocation == 0);

	if (alloc->maxAllocation > alloc->capacity)
	{
		s2Free(alloc->data, alloc->capacity);
		alloc->capacity = alloc->maxAllocation + alloc->maxAllocation / 2;
		alloc->data = s2Alloc(alloc->capacity);
	}
}

int32_t s2GetStackCapacity(s2StackAllocator* alloc)
{
	return alloc->capacity;
}

int32_t s2GetStackAllocation(s2StackAllocator* alloc)
{
	return alloc->allocation;
}

int32_t s2GetMaxStackAllocation(s2StackAllocator* alloc)
{
	return alloc->maxAllocation;
}
