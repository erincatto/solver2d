// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver2d/types.h"

typedef struct s2BlockAllocator s2BlockAllocator;

/// Create an allocator suitable for allocating and freeing small objects quickly.
/// Does not return memory to the heap.
s2BlockAllocator* s2CreateBlockAllocator(void);

/// Destroy a block alloctor instance
void s2DestroyBlockAllocator(s2BlockAllocator* allocator);

/// Allocate memory. This will use malloc if the size is larger than s2_maxBlockSize.
void* s2AllocBlock(s2BlockAllocator* allocator, int32_t size);

/// Free memory. This will use free if the size is larger than s2_maxBlockSize.
void s2FreeBlock(s2BlockAllocator* allocator,void* p, int32_t size);
