// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver2d/constants.h"
#include "solver2d/types.h"

#define s2_defaultCategoryBits (0x00000001)
#define s2_defaultMaskBits (0xFFFFFFFF)

// A node in the dynamic tree. The client does not interact with this directly.
// 16 + 16 + 8 + pad(8)
typedef struct s2TreeNode
{
	s2Box aabb; // 16

	// Category bits for collision filtering
	uint32_t categoryBits; // 4

	union
	{
		int32_t parent;
		int32_t next;
	}; // 4

	int32_t child1; // 4
	int32_t child2; // 4

	// TODO_ERIN could be union with child index
	int32_t userData; // 4

	// leaf = 0, free node = -1
	int16_t height; // 2

	bool enlarged; // 1

	char pad[9];
} s2TreeNode;

/// A dynamic AABB tree broad-phase, inspired by Nathanael Presson's btDbvt.
/// A dynamic tree arranges data in a binary tree to accelerate
/// queries such as volume queries and ray casts. Leaf nodes are proxies
/// with an AABB. These are used to hold a user collision object, such as a reference to a s2Shape.
///
/// Nodes are pooled and relocatable, so I use node indices rather than pointers.
typedef struct s2DynamicTree
{
	s2TreeNode* nodes;

	int32_t root;
	int32_t nodeCount;
	int32_t nodeCapacity;
	int32_t freeList;
	int32_t proxyCount;

	int32_t* leafIndices;
	s2Box* leafBoxes;
	s2Vec2* leafCenters;
	int32_t* binIndices;
	int32_t rebuildCapacity;
} s2DynamicTree;

#ifdef __cplusplus
extern "C"
{
#endif

/// Constructing the tree initializes the node pool.
s2DynamicTree s2DynamicTree_Create(void);

/// Destroy the tree, freeing the node pool.
void s2DynamicTree_Destroy(s2DynamicTree* tree);

/// Create a proxy. Provide a tight fitting AABB and a userData value.
int32_t s2DynamicTree_CreateProxy(s2DynamicTree* tree, s2Box aabb, uint32_t categoryBits, int32_t userData);

/// Destroy a proxy. This asserts if the id is invalid.
void s2DynamicTree_DestroyProxy(s2DynamicTree* tree, int32_t proxyId);

// Clone one tree to another, reusing storage in the outTree if possible
void s2DynamicTree_Clone(s2DynamicTree* outTree, const s2DynamicTree* inTree);

/// Move a proxy to a new AABB by removing and reinserting into the tree.
void s2DynamicTree_MoveProxy(s2DynamicTree* tree, int32_t proxyId, s2Box aabb);

/// Enlarge a proxy and enlarge ancestors as necessary.
void s2DynamicTree_EnlargeProxy(s2DynamicTree* tree, int32_t proxyId, s2Box aabb);

/// This function receives proxies found in the AABB query.
/// @return true if the query should continue
typedef bool s2TreeQueryCallbackFcn(int32_t proxyId, int32_t userData, void* context);

/// Query an AABB for overlapping proxies. The callback class
/// is called for each proxy that overlaps the supplied AABB.
void s2DynamicTree_QueryFiltered(const s2DynamicTree* tree, s2Box aabb, uint32_t maskBits, s2TreeQueryCallbackFcn* callback,
								 void* context);

/// Query an AABB for overlapping proxies. The callback class
/// is called for each proxy that overlaps the supplied AABB.
void s2DynamicTree_Query(const s2DynamicTree* tree, s2Box aabb, s2TreeQueryCallbackFcn* callback, void* context);

/// This function receives clipped raycast input for a proxy. The function
/// returns the new ray fraction.
/// - return a value of 0 to terminate the ray cast
/// - return a value less than input->maxFraction to clip the ray
/// - return a value of input->maxFraction to continue the ray cast without clipping
typedef float s2TreeRayCastCallbackFcn(const s2RayCastInput* input, int32_t proxyId, int32_t userData, void* context);

/// Ray-cast against the proxies in the tree. This relies on the callback
/// to perform a exact ray-cast in the case were the proxy contains a shape.
/// The callback also performs the any collision filtering. This has performance
/// roughly equal to k * log(n), where k is the number of collisions and n is the
/// number of proxies in the tree.
/// @param input the ray-cast input data. The ray extends from p1 to p1 + maxFraction * (p2 - p1).
/// @param callback a callback class that is called for each proxy that is hit by the ray.
void s2DynamicTree_RayCast(const s2DynamicTree* tree, const s2RayCastInput* input, uint32_t maskBits, s2TreeRayCastCallbackFcn* callback,
						   void* context);

/// Validate this tree. For testing.
void s2DynamicTree_Validate(const s2DynamicTree* tree);

/// Compute the height of the binary tree in O(N) time. Should not be
/// called often.
int32_t s2DynamicTree_GetHeight(const s2DynamicTree* tree);

/// Get the maximum balance of the tree. The balance is the difference in height of the two children of a node.
int32_t s2DynamicTree_GetMaxBalance(const s2DynamicTree* tree);

/// Get the ratio of the sum of the node areas to the root area.
float s2DynamicTree_GetAreaRatio(const s2DynamicTree* tree);

/// Build an optimal tree. Very expensive. For testing.
void s2DynamicTree_RebuildBottomUp(s2DynamicTree* tree);

/// Get the number of proxies created
int32_t s2DynamicTree_GetProxyCount(const s2DynamicTree* tree);

/// Rebuild the tree while retaining subtrees that haven't changed. Returns the number of boxes sorted.
int32_t s2DynamicTree_Rebuild(s2DynamicTree* tree, bool fullBuild);

/// Shift the world origin. Useful for large worlds.
/// The shift formula is: position -= newOrigin
/// @param newOrigin the new origin with respect to the old origin
void s2DynamicTree_ShiftOrigin(s2DynamicTree* tree, s2Vec2 newOrigin);

#ifdef __cplusplus
}
#endif

/// Get proxy user data.
/// @return the proxy user data or 0 if the id is invalid.
static inline int32_t s2DynamicTree_GetUserData(const s2DynamicTree* tree, int32_t proxyId)
{
	return tree->nodes[proxyId].userData;
}

static inline s2Box s2DynamicTree_GetAABB(const s2DynamicTree* tree, int32_t proxyId)
{
	return tree->nodes[proxyId].aabb;
}

static inline uint32_t s2DynamicTree_GetCategoryBits(const s2DynamicTree* tree, int32_t proxyId)
{
	return tree->nodes[proxyId].categoryBits;
}
