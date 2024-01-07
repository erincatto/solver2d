// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "core.h"

#include <stdbool.h>
#include <stdint.h>

// Bit set provides fast operations on large arrays of bits
typedef struct s2BitSet
{
	uint64_t* bits;
	uint32_t wordCapacity;
	uint32_t wordCount;
} s2BitSet;

s2BitSet s2CreateBitSet(uint32_t bitCapacity);
void s2DestroyBitSet(s2BitSet* bitSet);
void s2SetBitCountAndClear(s2BitSet* bitset, uint32_t bitCount);
void s2InPlaceUnion(s2BitSet* setA, const s2BitSet* setB);
void s2GrowBitSet(s2BitSet* set, uint32_t wordCount);

static inline void s2SetBit(s2BitSet* bitSet, uint32_t bitIndex)
{
	uint32_t wordIndex = bitIndex / 64;
	// TODO_ERIN support growing
	S2_ASSERT(wordIndex < bitSet->wordCount);
	bitSet->bits[wordIndex] |= ((uint64_t)1 << bitIndex % 64);
}

static inline void s2SetBitGrow(s2BitSet* bitSet, uint32_t bitIndex)
{
	uint32_t wordIndex = bitIndex / 64;
	if (wordIndex >= bitSet->wordCount)
	{
		s2GrowBitSet(bitSet, wordIndex + 1);
	}
	bitSet->bits[wordIndex] |= ((uint64_t)1 << bitIndex % 64);
}

static inline void s2ClearBit(s2BitSet* bitSet, uint32_t bitIndex)
{
	uint32_t wordIndex = bitIndex / 64;
	if (wordIndex >= bitSet->wordCount)
	{
		return;
	}
	bitSet->bits[wordIndex] &= ~((uint64_t)1 << bitIndex % 64);
}

static inline bool s2GetBit(const s2BitSet* bitSet, uint32_t bitIndex)
{
	uint32_t wordIndex = bitIndex / 64;
	if (wordIndex >= bitSet->wordCount)
	{
		return false;
	}
	return (bitSet->bits[wordIndex] & ((uint64_t)1 << bitIndex % 64)) != 0;
}

#if defined(_MSC_VER) && !defined(__clang__)
#include <intrin.h>

// https://en.wikipedia.org/wiki/Find_first_set
static inline uint32_t s2CTZ(uint64_t word)
{
	unsigned long index;

#ifdef _WIN64
	_BitScanForward64(&index, word);
#else
	// 32-bit fall back
	if ((uint32_t)word != 0)
	{
		_BitScanForward(&index, (uint32_t)word);
	}
	else
	{
		_BitScanForward(&index, (uint32_t)(word >> 32));
		index += 32;
	}
#endif

	return index;
}

#else

static inline uint32_t s2CTZ(uint64_t word)
{
	return __builtin_ctzll(word);
}

#endif
