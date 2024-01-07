// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "array.h"
#include "table.h"

#include "solver2d/dynamic_tree.h"

typedef struct s2Shape s2Shape;
typedef struct s2MovePair s2MovePair;
typedef struct s2MoveResult s2MoveResult;
typedef struct s2StackAllocator s2StackAllocator;
typedef struct s2World s2World;

// Store the proxy type in the lower 4 bits of the proxy key. This leaves 28 bits for the id.
#define S2_PROXY_TYPE(KEY) ((s2BodyType)((KEY)&0xF))
#define S2_PROXY_ID(KEY) ((KEY) >> 4)
#define S2_PROXY_KEY(ID, TYPE) (((ID) << 4) | (TYPE))

/// The broad-phase is used for computing pairs and performing volume queries and ray casts.
/// This broad-phase does not persist pairs. Instead, this reports potentially new pairs.
/// It is up to the client to consume the new pairs and to track subsequent overlap.
typedef struct s2BroadPhase
{
	s2DynamicTree trees[s2_bodyTypeCount];
	int32_t proxyCount;

	// The move set and array are used to track shapes that have moved significantly
	// and need a pair query for new contacts. The array has a deterministic order.
	// TODO_ERIN perhaps just a move set?
	s2Set moveSet;
	int32_t* moveArray;

	// These are the results from the pair query and are used to create new contacts
	// in deterministic order.
	s2MoveResult* moveResults;
	s2MovePair* movePairs;
	int32_t movePairCapacity;
	int movePairIndex;

	s2Set pairSet;

} s2BroadPhase;

void s2CreateBroadPhase(s2BroadPhase* bp);
void s2DestroyBroadPhase(s2BroadPhase* bp);
int32_t s2BroadPhase_CreateProxy(s2BroadPhase* bp, s2BodyType bodyType, s2Box aabb, uint32_t categoryBits, int32_t shapeIndex);
void s2BroadPhase_DestroyProxy(s2BroadPhase* bp, int32_t proxyKey);

void s2BroadPhase_MoveProxy(s2BroadPhase* bp, int32_t proxyKey, s2Box aabb);
void s2BroadPhase_EnlargeProxy(s2BroadPhase* bp, int32_t proxyKey, s2Box aabb);

void s2BroadPhase_RebuildTrees(s2BroadPhase* bp);

int32_t s2BroadPhase_GetShapeIndex(s2BroadPhase* bp, int32_t proxyKey);

void s2UpdateBroadPhasePairs(s2World* world);
bool s2BroadPhase_TestOverlap(const s2BroadPhase* bp, int32_t proxyKeyA, int32_t proxyKeyB);

void s2ValidateBroadphase(const s2BroadPhase* bp);
void s2ValidateNoEnlarged(const s2BroadPhase* bp);

// Warning: this must be called in deterministic order
static inline void s2BufferMove(s2BroadPhase* bp, int32_t proxyKey)
{
	bool alreadyAdded = s2AddKey(&bp->moveSet, proxyKey);
	if (alreadyAdded == false)
	{
		s2Array_Push(bp->moveArray, proxyKey);
	}
}
