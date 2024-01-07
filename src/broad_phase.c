// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#define _CRT_SECURE_NO_WARNINGS

#include "broad_phase.h"

#include "allocate.h"
#include "array.h"
#include "body.h"
#include "contact.h"
#include "core.h"
#include "shape.h"
#include "stack_allocator.h"
#include "world.h"

#include "solver2d/aabb.h"
#include "solver2d/timer.h"

#include <stdbool.h>
#include <string.h>

//#include <stdio.h>

//static FILE* s_file = NULL;

void s2CreateBroadPhase(s2BroadPhase* bp)
{
	//if (s_file == NULL)
	//{
	//	s_file = fopen("pairs01.txt", "a");
	//	fprintf(s_file, "============\n\n");
	//}

	bp->proxyCount = 0;

	// TODO_ERIN initial size in s2WorldDef?
	bp->moveSet = s2CreateSet(16);
	bp->moveArray = s2CreateArray(sizeof(int), 16);

	bp->moveResults = NULL;
	bp->movePairs = NULL;
	bp->movePairCapacity = 0;
	bp->movePairIndex = 0;

	// TODO_ERIN initial size from s2WorldDef
	bp->pairSet = s2CreateSet(32);

	for (int i = 0; i < s2_bodyTypeCount; ++i)
	{
		bp->trees[i] = s2DynamicTree_Create();
	}
}

void s2DestroyBroadPhase(s2BroadPhase* bp)
{
	for (int i = 0; i < s2_bodyTypeCount; ++i)
	{
		s2DynamicTree_Destroy(bp->trees + i);
	}

	s2DestroySet(&bp->moveSet);
	s2DestroyArray(bp->moveArray, sizeof(int));
	s2DestroySet(&bp->pairSet);

	memset(bp, 0, sizeof(s2BroadPhase));

	//if (s_file != NULL)
	//{
	//	fclose(s_file);
	//	s_file = NULL;
	//}
}

static inline void s2UnBufferMove(s2BroadPhase* bp, int proxyKey)
{
	bool found = s2RemoveKey(&bp->moveSet, proxyKey);

	if (found)
	{
		// Purge from move buffer. Linear search.
		int count = s2Array(bp->moveArray).count;
		for (int i = 0; i < count; ++i)
		{
			if (bp->moveArray[i] == proxyKey)
			{
				s2Array_RemoveSwap(bp->moveArray, i);
				break;
			}
		}
	}
}

int s2BroadPhase_CreateProxy(s2BroadPhase* bp, s2BodyType bodyType, s2Box aabb, uint32_t categoryBits, int shapeIndex)
{
	S2_ASSERT(0 <= bodyType && bodyType < s2_bodyTypeCount);
	int proxyId = s2DynamicTree_CreateProxy(bp->trees + bodyType, aabb, categoryBits, shapeIndex);
	int proxyKey = S2_PROXY_KEY(proxyId, bodyType);
	if (bodyType != s2_staticBody)
	{
		s2BufferMove(bp, proxyKey);
	}
	return proxyKey;
}

void s2BroadPhase_DestroyProxy(s2BroadPhase* bp, int proxyKey)
{
	S2_ASSERT(s2Array(bp->moveArray).count == (int)bp->moveSet.count);
	s2UnBufferMove(bp, proxyKey);

	--bp->proxyCount;

	int typeIndex = S2_PROXY_TYPE(proxyKey);
	int proxyId = S2_PROXY_ID(proxyKey);

	S2_ASSERT(0 <= typeIndex && typeIndex <= s2_bodyTypeCount);
	s2DynamicTree_DestroyProxy(bp->trees + typeIndex, proxyId);
}

void s2BroadPhase_MoveProxy(s2BroadPhase* bp, int proxyKey, s2Box aabb)
{
	s2BodyType bodyType = S2_PROXY_TYPE(proxyKey);
	int proxyId = S2_PROXY_ID(proxyKey);

	s2DynamicTree_MoveProxy(bp->trees + bodyType, proxyId, aabb);
	if (bodyType != s2_staticBody)
	{
		s2BufferMove(bp, proxyKey);
	}
}

void s2BroadPhase_EnlargeProxy(s2BroadPhase* bp, int proxyKey, s2Box aabb)
{
	int typeIndex = S2_PROXY_TYPE(proxyKey);
	int proxyId = S2_PROXY_ID(proxyKey);

	S2_ASSERT(typeIndex == s2_dynamicBody || typeIndex == s2_kinematicBody);

	s2DynamicTree_EnlargeProxy(bp->trees + typeIndex, proxyId, aabb);
	s2BufferMove(bp, proxyKey);
}

typedef struct s2MovePair
{
	int shapeIndexA;
	int shapeIndexB;
	s2MovePair* next;
	bool heap;
} s2MovePair;

typedef struct s2MoveResult
{
	s2MovePair* pairList;
} s2MoveResult;

typedef struct s2QueryPairContext
{
	s2World* world;
	s2MoveResult* moveResult;
	s2BodyType queryTreeType;
	int queryProxyKey;
	int queryShapeIndex;
} s2QueryPairContext;

// This is called from s2DynamicTree::Query when we are gathering pairs.
static bool s2PairQueryCallback(int proxyId, int shapeIndex, void* context)
{
	s2QueryPairContext* queryContext = context;
	s2BroadPhase* bp = &queryContext->world->broadPhase;

	int proxyKey = S2_PROXY_KEY(proxyId, queryContext->queryTreeType);

	// A proxy cannot form a pair with itself.
	if (proxyKey == queryContext->queryProxyKey)
	{
		return true;
	}

	bool moved = s2ContainsKey(&bp->moveSet, proxyKey);
	if (moved && proxyKey > queryContext->queryProxyKey)
	{
		// Both proxies are moving. Avoid duplicate pairs.
		return true;
	}

	uint64_t pairKey = S2_SHAPE_PAIR_KEY(shapeIndex, queryContext->queryShapeIndex);
	if (s2ContainsKey(&bp->pairSet, pairKey))
	{
		// contact exists
		return true;
	}

	int shapeIndexA, shapeIndexB;
	if (proxyKey < queryContext->queryProxyKey)
	{
		shapeIndexA = shapeIndex;
		shapeIndexB = queryContext->queryShapeIndex;
	}
	else
	{
		shapeIndexA = queryContext->queryShapeIndex;
		shapeIndexB = shapeIndex;
	}

	s2World* world = queryContext->world;

	S2_ASSERT(0 <= shapeIndexA && shapeIndexA < world->shapePool.capacity);
	S2_ASSERT(0 <= shapeIndexB && shapeIndexB < world->shapePool.capacity);

	s2Shape* shapeA = world->shapes + shapeIndexA;
	s2Shape* shapeB = world->shapes + shapeIndexB;

	// Are the fixtures on the same body?
	if (shapeA->bodyIndex == shapeB->bodyIndex)
	{
		return true;
	}

	int bodyIndexA = shapeA->bodyIndex;
	int bodyIndexB = shapeB->bodyIndex;
	s2Body* bodyA = world->bodies + bodyIndexA;
	s2Body* bodyB = world->bodies + bodyIndexB;

	// Does a joint override collision? Is at least one body dynamic?
	// TODO_ERIN this could be a hash set
	if (s2ShouldBodiesCollide(world, bodyA, bodyB) == false)
	{
		return true;
	}

	int pairIndex = bp->movePairIndex;
	bp->movePairIndex += 1;
	
	s2MovePair* pair;
	if (pairIndex < bp->movePairCapacity)
	{
		pair = bp->movePairs + pairIndex;
		pair->heap = false;
	}
	else
	{
		pair = malloc(sizeof(s2MovePair));
		pair->heap = true;
	}

	pair->shapeIndexA = shapeIndexA;
	pair->shapeIndexB = shapeIndexB;
	pair->next = queryContext->moveResult->pairList;
	queryContext->moveResult->pairList = pair;

	// continue the query
	return true;
}

void s2FindPairs(int count, s2World* world)
{
	s2BroadPhase* bp = &world->broadPhase;

	s2QueryPairContext queryContext;
	queryContext.world = world;

	for (int i = 0; i < count; ++i)
	{
		// Initialize move result for this moved proxy
		queryContext.moveResult = bp->moveResults + i;
		queryContext.moveResult->pairList = NULL;

		int proxyKey = bp->moveArray[i];
		if (proxyKey == S2_NULL_INDEX)
		{
			// proxy was destroyed after it moved
			continue;
		}

		int proxyType = S2_PROXY_TYPE(proxyKey);
		int proxyId = S2_PROXY_ID(proxyKey);
		queryContext.queryProxyKey = proxyKey;

		const s2DynamicTree* baseTree = bp->trees + proxyType;

		// We have to query the tree with the fat AABB so that
		// we don't fail to create a contact that may touch later.
		s2Box fatAABB = s2DynamicTree_GetAABB(baseTree, proxyId);
		queryContext.queryShapeIndex = s2DynamicTree_GetUserData(baseTree, proxyId);

		// Query trees
		if (proxyType == s2_dynamicBody)
		{
			queryContext.queryTreeType = s2_dynamicBody;
			s2DynamicTree_Query(bp->trees + s2_dynamicBody, fatAABB, s2PairQueryCallback, &queryContext);
			queryContext.queryTreeType = s2_kinematicBody;
			s2DynamicTree_Query(bp->trees + s2_kinematicBody, fatAABB, s2PairQueryCallback, &queryContext);
			queryContext.queryTreeType = s2_staticBody;
			s2DynamicTree_Query(bp->trees + s2_staticBody, fatAABB, s2PairQueryCallback, &queryContext);
		}
		else if (proxyType == s2_kinematicBody)
		{
			queryContext.queryTreeType = s2_dynamicBody;
			s2DynamicTree_Query(bp->trees + s2_dynamicBody, fatAABB, s2PairQueryCallback, &queryContext);
		}
	}
}

void s2UpdateBroadPhasePairs(s2World* world)
{
	s2BroadPhase* bp = &world->broadPhase;

	int moveCount = s2Array(bp->moveArray).count;
	S2_ASSERT(moveCount == (int)bp->moveSet.count);

	if (moveCount == 0)
	{
		return;
	}

	s2StackAllocator* alloc = world->stackAllocator;

	bp->moveResults = s2AllocateStackItem(alloc, moveCount * sizeof(s2MoveResult), "move results");
	bp->movePairCapacity = 16 * moveCount;
	bp->movePairs = s2AllocateStackItem(alloc, bp->movePairCapacity * sizeof(s2MovePair), "move pairs");
	bp->movePairIndex = 0;

	s2FindPairs(moveCount, world);

	s2Shape* shapes = world->shapes;

	for (int i = 0; i < moveCount; ++i)
	{
		s2MoveResult* result = bp->moveResults + i;
		s2MovePair* pair = result->pairList;
		while (pair != NULL)
		{
			int shapeIndexA = pair->shapeIndexA;
			int shapeIndexB = pair->shapeIndexB;

			S2_ASSERT(0 <= shapeIndexA && shapeIndexA < world->shapePool.capacity);
			S2_ASSERT(0 <= shapeIndexB && shapeIndexB < world->shapePool.capacity);

			s2CreateContact(world, shapes + shapeIndexA, shapes + shapeIndexB);

			if (pair->heap)
			{
				s2MovePair* temp = pair;
				pair = pair->next;
				s2Free(temp, sizeof(s2MovePair));
			}
			else
			{
				pair = pair->next;
			}
		}
	}

	// Reset move buffer
	s2Array_Clear(bp->moveArray);
	s2ClearSet(&bp->moveSet);

	s2FreeStackItem(alloc, bp->movePairs);
	bp->movePairs = NULL;
	s2FreeStackItem(alloc, bp->moveResults);
	bp->moveResults = NULL;
}

bool s2BroadPhase_TestOverlap(const s2BroadPhase* bp, int proxyKeyA, int proxyKeyB)
{
	int typeIndexA = S2_PROXY_TYPE(proxyKeyA);
	int proxyIdA = S2_PROXY_ID(proxyKeyA);
	int typeIndexB = S2_PROXY_TYPE(proxyKeyB);
	int proxyIdB = S2_PROXY_ID(proxyKeyB);

	s2Box aabbA = s2DynamicTree_GetAABB(bp->trees + typeIndexA, proxyIdA);
	s2Box aabbB = s2DynamicTree_GetAABB(bp->trees + typeIndexB, proxyIdB);
	return s2AABB_Overlaps(aabbA, aabbB);
}

void s2BroadPhase_RebuildTrees(s2BroadPhase* bp)
{
	s2DynamicTree_Rebuild(bp->trees + s2_dynamicBody, false);
	s2DynamicTree_Rebuild(bp->trees + s2_kinematicBody, false);
}

int s2BroadPhase_GetShapeIndex(s2BroadPhase* bp, int proxyKey)
{
	int typeIndex = S2_PROXY_TYPE(proxyKey);
	int proxyId = S2_PROXY_ID(proxyKey);

	return s2DynamicTree_GetUserData(bp->trees + typeIndex, proxyId);
}

void s2ValidateBroadphase(const s2BroadPhase* bp)
{
	s2DynamicTree_Validate(bp->trees + s2_dynamicBody);
	s2DynamicTree_Validate(bp->trees + s2_kinematicBody);

	// TODO_ERIN validate every shape AABB is contained in tree AABB
}

void s2ValidateNoEnlarged(const s2BroadPhase* bp)
{
#if S2_VALIDATE == 1
	for (int j = 0; j < s2_bodyTypeCount; ++j)
	{
		const s2DynamicTree* tree = bp->trees + j;
		int capacity = tree->nodeCapacity;
		const s2TreeNode* nodes = tree->nodes;
		for (int i = 0; i < capacity; ++i)
		{
			const s2TreeNode* node = nodes + i;
			if (node->height < 0)
			{
				continue;
			}

			if (node->enlarged == true)
			{
				capacity += 0;
			}

			S2_ASSERT(node->enlarged == false);
		}
	}
#else
	S2_MAYBE_UNUSED(bp);
#endif
}
