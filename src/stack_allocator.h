// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include <stdint.h>

typedef struct s2StackAllocator s2StackAllocator;

s2StackAllocator* s2CreateStackAllocator(int32_t capacity);
void s2DestroyStackAllocator(s2StackAllocator* allocator);

void* s2AllocateStackItem(s2StackAllocator* alloc, int32_t size, const char* name);
void s2FreeStackItem(s2StackAllocator* alloc, void* mem);

// Grow the stack based on usage
void s2GrowStack(s2StackAllocator* alloc);

int32_t s2GetStackCapacity(s2StackAllocator* alloc);
int32_t s2GetStackAllocation(s2StackAllocator* alloc);
int32_t s2GetMaxStackAllocation(s2StackAllocator* alloc);
