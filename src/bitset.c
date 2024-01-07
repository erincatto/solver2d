// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "bitset.h"

#include "allocate.h"

#include <string.h>

s2BitSet s2CreateBitSet(uint32_t bitCapacity)
{
	s2BitSet bitSet = {0};

	bitSet.wordCapacity = (bitCapacity + sizeof(uint64_t) * 8 - 1) / (sizeof(uint64_t) * 8);
	bitSet.wordCount = 0;
	bitSet.bits = s2Alloc(bitSet.wordCapacity * sizeof(uint64_t));
	memset(bitSet.bits, 0, bitSet.wordCapacity * sizeof(uint64_t));
	return bitSet;
}

void s2DestroyBitSet(s2BitSet* bitSet)
{
	s2Free(bitSet->bits, bitSet->wordCapacity * sizeof(uint64_t));
	bitSet->wordCapacity = 0;
	bitSet->wordCount = 0;
	bitSet->bits = NULL;
}

void s2SetBitCountAndClear(s2BitSet* bitSet, uint32_t bitCount)
{
	uint32_t wordCount = (bitCount + sizeof(uint64_t) * 8 - 1) / (sizeof(uint64_t) * 8);
	if (bitSet->wordCapacity < wordCount)
	{
		s2DestroyBitSet(bitSet);
		uint32_t newBitCapacity = bitCount + (bitCount >> 1);
		*bitSet = s2CreateBitSet(newBitCapacity);
	}

	bitSet->wordCount = wordCount;
	memset(bitSet->bits, 0, bitSet->wordCount * sizeof(uint64_t));
}

void s2GrowBitSet(s2BitSet* bitSet, uint32_t wordCount)
{
	S2_ASSERT(wordCount > bitSet->wordCount);
	if (wordCount > bitSet->wordCapacity)
	{
		uint32_t oldCapacity = bitSet->wordCapacity;
		bitSet->wordCapacity = wordCount + wordCount / 2;
		uint64_t* newBits = s2Alloc(bitSet->wordCapacity * sizeof(uint64_t));
		memset(newBits, 0, bitSet->wordCapacity * sizeof(uint64_t));
		memcpy(newBits, bitSet->bits, bitSet->wordCount * sizeof(uint64_t));
		s2Free(bitSet->bits, oldCapacity * sizeof(uint64_t));
		bitSet->bits = newBits;
	}

	bitSet->wordCount = wordCount;
}

void s2InPlaceUnion(s2BitSet* setA, const s2BitSet* setB)
{
	S2_ASSERT(setA->wordCount == setB->wordCount);
	uint32_t wordCount = setA->wordCount;
	for (uint32_t i = 0; i < wordCount; ++i)
	{
		setA->bits[i] |= setB->bits[i];
	}
}
