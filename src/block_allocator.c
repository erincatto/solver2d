// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "block_allocator.h"
#include "allocate.h"
#include "core.h"

#include <limits.h>
#include <string.h>

#define s2_chunkSize (16 * 1024)
#define s2_maxBlockSize 640
#define s2_chunkArrayIncrement 128

// These are the supported object sizes. Actual allocations are rounded up the next size.
static const int32_t s2_blockSizes[] = {
	16,	 // 0
	32,	 // 1
	64,	 // 2
	96,	 // 3
	128, // 4
	160, // 5
	192, // 6
	224, // 7
	256, // 8
	320, // 9
	384, // 10
	448, // 11
	512, // 12
	640, // 13
};

#define s2_blockSizeCount S2_ARRAY_COUNT(s2_blockSizes)

// This maps an arbitrary allocation size to a suitable slot in s2_blockSizes.
typedef struct s2SizeMap
{
	uint8_t values[s2_maxBlockSize + 1];
} s2SizeMap;

static s2SizeMap s2_sizeMap;

static void s2InitializeSizeMap(void)
{
	int32_t j = 0;
	s2_sizeMap.values[0] = 0;
	for (int32_t i = 1; i <= s2_maxBlockSize; ++i)
	{
		S2_ASSERT(j < s2_blockSizeCount);
		if (i <= s2_blockSizes[j])
		{
			s2_sizeMap.values[i] = (uint8_t)j;
		}
		else
		{
			++j;
			s2_sizeMap.values[i] = (uint8_t)j;
		}
	}
}

static bool s2_sizeMapInitialized = false;

typedef struct s2Chunk
{
	int32_t blockSize;
	struct s2Block* blocks;
} s2Chunk;

typedef struct s2Block
{
	struct s2Block* next;
} s2Block;

// This is a small object allocator used for allocating small objects that persist for more than one time step.
// See: http://www.codeproject.com/useritems/Small_Block_Allocator.asp
typedef struct s2BlockAllocator
{
	s2Chunk* chunks;
	int32_t chunkCount;
	int32_t chunkSpace;

	s2Block* freeLists[s2_blockSizeCount];
} s2BlockAllocator;

s2BlockAllocator* s2CreateBlockAllocator(void)
{
	if (s2_sizeMapInitialized == false)
	{
		s2InitializeSizeMap();
		s2_sizeMapInitialized = true;
	}

	_Static_assert(s2_blockSizeCount < UCHAR_MAX, "block size too large");

	s2BlockAllocator* allocator = (s2BlockAllocator*)s2Alloc(sizeof(s2BlockAllocator));
	allocator->chunkSpace = s2_chunkArrayIncrement;
	allocator->chunkCount = 0;
	allocator->chunks = (s2Chunk*)s2Alloc(allocator->chunkSpace * sizeof(s2Chunk));

	memset(allocator->chunks, 0, allocator->chunkSpace * sizeof(s2Chunk));
	memset(allocator->freeLists, 0, sizeof(allocator->freeLists));

	return allocator;
}

void s2DestroyBlockAllocator(s2BlockAllocator* allocator)
{
	for (int32_t i = 0; i < allocator->chunkCount; ++i)
	{
		s2Free(allocator->chunks[i].blocks, s2_chunkSize);
	}

	s2Free(allocator->chunks, allocator->chunkSpace * sizeof(s2Chunk));
	s2Free(allocator, sizeof(s2BlockAllocator));
}

void* s2AllocBlock(s2BlockAllocator* allocator, int32_t size)
{
	if (size == 0)
	{
		return NULL;
	}

	S2_ASSERT(0 < size);

	if (size > s2_maxBlockSize)
	{
		return s2Alloc(size);
	}

	int32_t index = s2_sizeMap.values[size];
	S2_ASSERT(0 <= index && index < s2_blockSizeCount);

	if (allocator->freeLists[index])
	{
		s2Block* block = allocator->freeLists[index];
		allocator->freeLists[index] = block->next;
		return block;
	}
	else
	{
		if (allocator->chunkCount == allocator->chunkSpace)
		{
			s2Chunk* oldChunks = allocator->chunks;
			int32_t oldSize = allocator->chunkSpace * sizeof(s2Chunk);
			allocator->chunkSpace += s2_chunkArrayIncrement;
			allocator->chunks = (s2Chunk*)s2Alloc(allocator->chunkSpace * sizeof(s2Chunk));
			memcpy(allocator->chunks, oldChunks, allocator->chunkCount * sizeof(s2Chunk));
			memset(allocator->chunks + allocator->chunkCount, 0, s2_chunkArrayIncrement * sizeof(s2Chunk));
			s2Free(oldChunks, oldSize);
		}

		s2Chunk* chunk = allocator->chunks + allocator->chunkCount;
		chunk->blocks = (s2Block*)s2Alloc(s2_chunkSize);
#if defined(_DEBUG)
		memset(chunk->blocks, 0xcd, s2_chunkSize);
#endif
		int32_t blockSize = s2_blockSizes[index];
		chunk->blockSize = blockSize;
		int32_t blockCount = s2_chunkSize / blockSize;
		S2_ASSERT(blockCount * blockSize <= s2_chunkSize);
		for (int32_t i = 0; i < blockCount - 1; ++i)
		{
			s2Block* block = (s2Block*)((int8_t*)chunk->blocks + blockSize * i);
			s2Block* next = (s2Block*)((int8_t*)chunk->blocks + blockSize * (i + 1));
			block->next = next;
		}
		s2Block* last = (s2Block*)((int8_t*)chunk->blocks + blockSize * (blockCount - 1));
		last->next = NULL;

		allocator->freeLists[index] = chunk->blocks->next;
		++allocator->chunkCount;

		return chunk->blocks;
	}
}

void s2FreeBlock(s2BlockAllocator* allocator, void* p, int32_t size)
{
	if (size == 0)
	{
		return;
	}

	S2_ASSERT(0 < size);

	if (size > s2_maxBlockSize)
	{
		s2Free(p, size);
		return;
	}

	int32_t index = s2_sizeMap.values[size];
	S2_ASSERT(0 <= index && index < s2_blockSizeCount);

#if defined(_DEBUG)
	// Verify the memory address and size is valid.
	int32_t blockSize = s2_blockSizes[index];
	bool found = false;
	for (int32_t i = 0; i < allocator->chunkCount; ++i)
	{
		s2Chunk* chunk = allocator->chunks + i;
		if (chunk->blockSize != blockSize)
		{
			S2_ASSERT((int8_t*)p + blockSize <= (int8_t*)chunk->blocks || (int8_t*)chunk->blocks + s2_chunkSize <= (int8_t*)p);
		}
		else
		{
			if ((int8_t*)chunk->blocks <= (int8_t*)p && (int8_t*)p + blockSize <= (int8_t*)chunk->blocks + s2_chunkSize)
			{
				found = true;
			}
		}
	}

	S2_ASSERT(found);

	memset(p, 0xfd, blockSize);
#endif

	s2Block* block = (s2Block*)p;
	block->next = allocator->freeLists[index];
	allocator->freeLists[index] = block;
}
