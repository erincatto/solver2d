// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "solver2d/dynamic_tree.h"

#include "allocate.h"
#include "array.h"
#include "core.h"

#include "solver2d/aabb.h"
#include "solver2d/constants.h"

#include <float.h>
#include <string.h>

// TODO_ERIN
// - try incrementally sorting internal nodes by height for better cache efficiency during depth first traversal.

static s2TreeNode s2_defaultTreeNode = {{{0.0f, 0.0f}, {0.0f, 0.0f}}, 0, {S2_NULL_INDEX}, S2_NULL_INDEX, S2_NULL_INDEX, -1, -2, false,
										{0, 0, 0, 0, 0, 0, 0, 0, 0}};

static inline bool s2IsLeaf(const s2TreeNode* node)
{
	return node->height == 0;
}

s2DynamicTree s2DynamicTree_Create(void)
{
	_Static_assert((sizeof(s2TreeNode) & 0xF) == 0, "tree node size not a multiple of 16");

	s2DynamicTree tree;
	tree.root = S2_NULL_INDEX;

	tree.nodeCapacity = 16;
	tree.nodeCount = 0;
	tree.nodes = (s2TreeNode*)s2Alloc(tree.nodeCapacity * sizeof(s2TreeNode));
	memset(tree.nodes, 0, tree.nodeCapacity * sizeof(s2TreeNode));

	// Build a linked list for the free list.
	for (int32_t i = 0; i < tree.nodeCapacity - 1; ++i)
	{
		tree.nodes[i].next = i + 1;
		tree.nodes[i].height = -1;
	}
	tree.nodes[tree.nodeCapacity - 1].next = S2_NULL_INDEX;
	tree.nodes[tree.nodeCapacity - 1].height = -1;
	tree.freeList = 0;

	tree.proxyCount = 0;

	tree.leafIndices = NULL;
	tree.leafBoxes = NULL;
	tree.leafCenters = NULL;
	tree.binIndices = NULL;
	tree.rebuildCapacity = 0;

	return tree;
}

void s2DynamicTree_Destroy(s2DynamicTree* tree)
{
	s2Free(tree->nodes, tree->nodeCapacity * sizeof(s2TreeNode));
	s2Free(tree->leafIndices, tree->rebuildCapacity * sizeof(int32_t));
	s2Free(tree->leafBoxes, tree->rebuildCapacity * sizeof(s2Box));
	s2Free(tree->leafCenters, tree->rebuildCapacity * sizeof(s2Vec2));
	s2Free(tree->binIndices, tree->rebuildCapacity * sizeof(int32_t));

	memset(tree, 0, sizeof(s2DynamicTree));
}

void s2DynamicTree_Clone(s2DynamicTree* outTree, const s2DynamicTree* inTree)
{
	if (outTree->nodeCapacity < inTree->nodeCapacity)
	{
		s2Free(outTree->nodes, outTree->nodeCapacity * sizeof(s2TreeNode));
		outTree->nodeCapacity = inTree->nodeCapacity;
		outTree->nodes = (s2TreeNode*)s2Alloc(outTree->nodeCapacity * sizeof(s2TreeNode));
	}

	memcpy(outTree->nodes, inTree->nodes, inTree->nodeCapacity * sizeof(s2TreeNode));
	outTree->root = inTree->root;
	outTree->nodeCount = inTree->nodeCount;
	outTree->freeList = inTree->freeList;
	outTree->proxyCount = inTree->proxyCount;

	// Hook up free list.
	// TODO_ERIN make this optional?
	// TODO_ERIN perhaps find tail of existing free list and append
	int32_t inCapacity = inTree->nodeCapacity;
	int32_t outCapacity = outTree->nodeCapacity;
	if (outCapacity > inCapacity)
	{
		for (int32_t i = inCapacity; i < outCapacity - 1; ++i)
		{
			outTree->nodes[i].next = i + 1;
			outTree->nodes[i].height = -1;
		}
		outTree->nodes[outCapacity - 1].next = outTree->freeList;
		outTree->nodes[outCapacity - 1].height = -1;
		outTree->freeList = inCapacity;
	}
}

// Allocate a node from the pool. Grow the pool if necessary.
static int32_t s2AllocateNode(s2DynamicTree* tree)
{
	// Expand the node pool as needed.
	if (tree->freeList == S2_NULL_INDEX)
	{
		S2_ASSERT(tree->nodeCount == tree->nodeCapacity);

		// The free list is empty. Rebuild a bigger pool.
		s2TreeNode* oldNodes = tree->nodes;
		int32_t oldCapcity = tree->nodeCapacity;
		tree->nodeCapacity += oldCapcity >> 1;
		tree->nodes = (s2TreeNode*)s2Alloc(tree->nodeCapacity * sizeof(s2TreeNode));
		memcpy(tree->nodes, oldNodes, tree->nodeCount * sizeof(s2TreeNode));
		s2Free(oldNodes, oldCapcity * sizeof(s2TreeNode));

		// Build a linked list for the free list. The parent
		// pointer becomes the "next" pointer.
		for (int32_t i = tree->nodeCount; i < tree->nodeCapacity - 1; ++i)
		{
			tree->nodes[i].next = i + 1;
			tree->nodes[i].height = -1;
		}
		tree->nodes[tree->nodeCapacity - 1].next = S2_NULL_INDEX;
		tree->nodes[tree->nodeCapacity - 1].height = -1;
		tree->freeList = tree->nodeCount;
	}

	// Peel a node off the free list.
	int32_t nodeIndex = tree->freeList;
	s2TreeNode* node = tree->nodes + nodeIndex;
	tree->freeList = node->next;
	*node = s2_defaultTreeNode;
	++tree->nodeCount;
	return nodeIndex;
}

// Return a node to the pool.
static void s2FreeNode(s2DynamicTree* tree, int32_t nodeId)
{
	S2_ASSERT(0 <= nodeId && nodeId < tree->nodeCapacity);
	S2_ASSERT(0 < tree->nodeCount);
	tree->nodes[nodeId].next = tree->freeList;
	tree->nodes[nodeId].height = -1;
	tree->freeList = nodeId;
	--tree->nodeCount;
}

// Greedy algorithm for sibling selection using the SAH
// We have three nodes A-(B,C) and want to add a leaf D, there are three choices.
// 1: make a new parent for A and D : E-(A-(B,C), D)
// 2: associate D with B
//   a: B is a leaf : A-(E-(B,D), C)
//   b: B is an internal node: A-(B{D},C)
// 3: associate D with C
//   a: C is a leaf : A-(B, E-(C,D))
//   b: C is an internal node: A-(B, C{D})
// All of these have a clear cost except when B or C is an internal node. Hence we need to be greedy.

// The cost for cases 1, 2a, and 3a can be computed using the sibling cost formula.
// cost of sibling H = area(union(H, D)) + increased are of ancestors

// Suppose B (or C) is an internal node, then the lowest cost would be one of two cases:
// case1: D becomes a sibling of B
// case2: D becomes a descendant of B along with a new internal node of area(D).
static int32_t s2FindBestSibling(const s2DynamicTree* tree, s2Box boxD)
{
	s2Vec2 centerD = s2AABB_Center(boxD);
	float areaD = s2AABB_Perimeter(boxD);

	const s2TreeNode* nodes = tree->nodes;
	int32_t rootIndex = tree->root;

	s2Box rootBox = nodes[rootIndex].aabb;

	// Area of current node
	float areaBase = s2AABB_Perimeter(rootBox);

	// Area of inflated node
	float directCost = s2AABB_Perimeter(s2AABB_Union(rootBox, boxD));
	float inheritedCost = 0.0f;

	int32_t bestSibling = rootIndex;
	float bestCost = directCost;

	// Descend the tree from root, following a single greedy path.
	int32_t index = rootIndex;
	while (nodes[index].height > 0)
	{
		int32_t child1 = nodes[index].child1;
		int32_t child2 = nodes[index].child2;

		// Cost of creating a new parent for this node and the new leaf
		float cost = directCost + inheritedCost;

		// Sometimes there are multiple identical costs within tolerance.
		// This breaks the ties using the centroid distance.
		if (cost < bestCost)
		{
			bestSibling = index;
			bestCost = cost;
		}

		// Inheritance cost seen by children
		inheritedCost += directCost - areaBase;

		bool leaf1 = nodes[child1].height == 0;
		bool leaf2 = nodes[child2].height == 0;

		// Cost of descending into child 1
		float lowerCost1 = FLT_MAX;
		s2Box box1 = nodes[child1].aabb;
		float directCost1 = s2AABB_Perimeter(s2AABB_Union(box1, boxD));
		float area1 = 0.0f;
		if (leaf1)
		{
			// Child 1 is a leaf
			// Cost of creating new node and increasing area of node P
			float cost1 = directCost1 + inheritedCost;

			// Need this here due to while condition above
			if (cost1 < bestCost)
			{
				bestSibling = child1;
				bestCost = cost1;
			}
		}
		else
		{
			// Child 1 is an internal node
			area1 = s2AABB_Perimeter(box1);

			// Lower bound cost of inserting under child 1.
			lowerCost1 = inheritedCost + directCost1 + S2_MIN(areaD - area1, 0.0f);
		}

		// Cost of descending into child 2
		float lowerCost2 = FLT_MAX;
		s2Box box2 = nodes[child2].aabb;
		float directCost2 = s2AABB_Perimeter(s2AABB_Union(box2, boxD));
		float area2 = 0.0f;
		if (leaf2)
		{
			// Child 2 is a leaf
			// Cost of creating new node and increasing area of node P
			float cost2 = directCost2 + inheritedCost;

			// Need this here due to while condition above
			if (cost2 < bestCost)
			{
				bestSibling = child2;
				bestCost = cost2;
			}
		}
		else
		{
			// Child 2 is an internal node
			area2 = s2AABB_Perimeter(box2);

			// Lower bound cost of inserting under child 2. This is not the cost
			// of child 2, it is the best we can hope for under child 2.
			lowerCost2 = inheritedCost + directCost2 + S2_MIN(areaD - area2, 0.0f);
		}

		if (leaf1 && leaf2)
		{
			break;
		}

		// Can the cost possibly be decreased?
		if (bestCost <= lowerCost1 && bestCost <= lowerCost2)
		{
			break;
		}

		if (lowerCost1 == lowerCost2 && leaf1 == false)
		{
			S2_ASSERT(lowerCost1 < FLT_MAX);
			S2_ASSERT(lowerCost2 < FLT_MAX);

			// No clear choice based on lower bound surface area. This can happen when both
			// children fully contain D. Fall back to node distance.
			s2Vec2 d1 = s2Sub(s2AABB_Center(box1), centerD);
			s2Vec2 d2 = s2Sub(s2AABB_Center(box2), centerD);
			lowerCost1 = s2LengthSquared(d1);
			lowerCost2 = s2LengthSquared(d2);
		}

		// Descend
		if (lowerCost1 < lowerCost2 && leaf1 == false)
		{
			index = child1;
			areaBase = area1;
			directCost = directCost1;
		}
		else
		{
			index = child2;
			areaBase = area2;
			directCost = directCost2;
		}

		S2_ASSERT(nodes[index].height > 0);
	}

	return bestSibling;
}

enum s2RotateType
{
	s2_rotateNone,
	s2_rotateBF,
	s2_rotateBG,
	s2_rotateCD,
	s2_rotateCE
};

// Perform a left or right rotation if node A is imbalanced.
// Returns the new root index.
static void s2RotateNodes(s2DynamicTree* tree, int32_t iA)
{
	S2_ASSERT(iA != S2_NULL_INDEX);

	s2TreeNode* nodes = tree->nodes;

	s2TreeNode* A = nodes + iA;
	if (A->height < 2)
	{
		return;
	}

	int32_t iB = A->child1;
	int32_t iC = A->child2;
	S2_ASSERT(0 <= iB && iB < tree->nodeCapacity);
	S2_ASSERT(0 <= iC && iC < tree->nodeCapacity);

	s2TreeNode* B = nodes + iB;
	s2TreeNode* C = nodes + iC;

	if (B->height == 0)
	{
		// B is a leaf and C is internal
		S2_ASSERT(C->height > 0);

		int32_t iF = C->child1;
		int32_t iG = C->child2;
		s2TreeNode* F = nodes + iF;
		s2TreeNode* G = nodes + iG;
		S2_ASSERT(0 <= iF && iF < tree->nodeCapacity);
		S2_ASSERT(0 <= iG && iG < tree->nodeCapacity);

		// Base cost
		float costBase = s2AABB_Perimeter(C->aabb);

		// Cost of swapping B and F
		s2Box aabbBG = s2AABB_Union(B->aabb, G->aabb);
		float costBF = s2AABB_Perimeter(aabbBG);

		// Cost of swapping B and G
		s2Box aabbBF = s2AABB_Union(B->aabb, F->aabb);
		float costBG = s2AABB_Perimeter(aabbBF);

		if (costBase < costBF && costBase < costBG)
		{
			// Rotation does not improve cost
			return;
		}

		if (costBF < costBG)
		{
			// Swap B and F
			A->child1 = iF;
			C->child1 = iB;

			B->parent = iC;
			F->parent = iA;

			C->aabb = aabbBG;

			C->height = 1 + S2_MAX(B->height, G->height);
			A->height = 1 + S2_MAX(C->height, F->height);
			C->categoryBits = B->categoryBits | G->categoryBits;
			A->categoryBits = C->categoryBits | F->categoryBits;
			C->enlarged = B->enlarged || G->enlarged;
			A->enlarged = C->enlarged || F->enlarged;
		}
		else
		{
			// Swap B and G
			A->child1 = iG;
			C->child2 = iB;

			B->parent = iC;
			G->parent = iA;

			C->aabb = aabbBF;

			C->height = 1 + S2_MAX(B->height, F->height);
			A->height = 1 + S2_MAX(C->height, G->height);
			C->categoryBits = B->categoryBits | F->categoryBits;
			A->categoryBits = C->categoryBits | G->categoryBits;
			C->enlarged = B->enlarged || F->enlarged;
			A->enlarged = C->enlarged || G->enlarged;
		}
	}
	else if (C->height == 0)
	{
		// C is a leaf and B is internal
		S2_ASSERT(B->height > 0);

		int iD = B->child1;
		int iE = B->child2;
		s2TreeNode* D = nodes + iD;
		s2TreeNode* E = nodes + iE;
		S2_ASSERT(0 <= iD && iD < tree->nodeCapacity);
		S2_ASSERT(0 <= iE && iE < tree->nodeCapacity);

		// Base cost
		float costBase = s2AABB_Perimeter(B->aabb);

		// Cost of swapping C and D
		s2Box aabbCE = s2AABB_Union(C->aabb, E->aabb);
		float costCD = s2AABB_Perimeter(aabbCE);

		// Cost of swapping C and E
		s2Box aabbCD = s2AABB_Union(C->aabb, D->aabb);
		float costCE = s2AABB_Perimeter(aabbCD);

		if (costBase < costCD && costBase < costCE)
		{
			// Rotation does not improve cost
			return;
		}

		if (costCD < costCE)
		{
			// Swap C and D
			A->child2 = iD;
			B->child1 = iC;

			C->parent = iB;
			D->parent = iA;

			B->aabb = aabbCE;

			B->height = 1 + S2_MAX(C->height, E->height);
			A->height = 1 + S2_MAX(B->height, D->height);
			B->categoryBits = C->categoryBits | E->categoryBits;
			A->categoryBits = B->categoryBits | D->categoryBits;
			B->enlarged = C->enlarged || E->enlarged;
			A->enlarged = B->enlarged || D->enlarged;
		}
		else
		{
			// Swap C and E
			A->child2 = iE;
			B->child2 = iC;

			C->parent = iB;
			E->parent = iA;

			B->aabb = aabbCD;
			B->height = 1 + S2_MAX(C->height, D->height);
			A->height = 1 + S2_MAX(B->height, E->height);
			B->categoryBits = C->categoryBits | D->categoryBits;
			A->categoryBits = B->categoryBits | E->categoryBits;
			B->enlarged = C->enlarged || D->enlarged;
			A->enlarged = B->enlarged || E->enlarged;
		}
	}
	else
	{
		int iD = B->child1;
		int iE = B->child2;
		int iF = C->child1;
		int iG = C->child2;

		s2TreeNode* D = nodes + iD;
		s2TreeNode* E = nodes + iE;
		s2TreeNode* F = nodes + iF;
		s2TreeNode* G = nodes + iG;

		S2_ASSERT(0 <= iD && iD < tree->nodeCapacity);
		S2_ASSERT(0 <= iE && iE < tree->nodeCapacity);
		S2_ASSERT(0 <= iF && iF < tree->nodeCapacity);
		S2_ASSERT(0 <= iG && iG < tree->nodeCapacity);

		// Base cost
		float areaB = s2AABB_Perimeter(B->aabb);
		float areaC = s2AABB_Perimeter(C->aabb);
		float costBase = areaB + areaC;
		enum s2RotateType bestRotation = s2_rotateNone;
		float bestCost = costBase;

		// Cost of swapping B and F
		s2Box aabbBG = s2AABB_Union(B->aabb, G->aabb);
		float costBF = areaB + s2AABB_Perimeter(aabbBG);
		if (costBF < bestCost)
		{
			bestRotation = s2_rotateBF;
			bestCost = costBF;
		}

		// Cost of swapping B and G
		s2Box aabbBF = s2AABB_Union(B->aabb, F->aabb);
		float costBG = areaB + s2AABB_Perimeter(aabbBF);
		if (costBG < bestCost)
		{
			bestRotation = s2_rotateBG;
			bestCost = costBG;
		}

		// Cost of swapping C and D
		s2Box aabbCE = s2AABB_Union(C->aabb, E->aabb);
		float costCD = areaC + s2AABB_Perimeter(aabbCE);
		if (costCD < bestCost)
		{
			bestRotation = s2_rotateCD;
			bestCost = costCD;
		}

		// Cost of swapping C and E
		s2Box aabbCD = s2AABB_Union(C->aabb, D->aabb);
		float costCE = areaC + s2AABB_Perimeter(aabbCD);
		if (costCE < bestCost)
		{
			bestRotation = s2_rotateCE;
			bestCost = costCE;
		}

		switch (bestRotation)
		{
		case s2_rotateNone:
			break;

		case s2_rotateBF:
			A->child1 = iF;
			C->child1 = iB;

			B->parent = iC;
			F->parent = iA;

			C->aabb = aabbBG;
			C->height = 1 + S2_MAX(B->height, G->height);
			A->height = 1 + S2_MAX(C->height, F->height);
			C->categoryBits = B->categoryBits | G->categoryBits;
			A->categoryBits = C->categoryBits | F->categoryBits;
			C->enlarged = B->enlarged || G->enlarged;
			A->enlarged = C->enlarged || F->enlarged;
			break;

		case s2_rotateBG:
			A->child1 = iG;
			C->child2 = iB;

			B->parent = iC;
			G->parent = iA;

			C->aabb = aabbBF;
			C->height = 1 + S2_MAX(B->height, F->height);
			A->height = 1 + S2_MAX(C->height, G->height);
			C->categoryBits = B->categoryBits | F->categoryBits;
			A->categoryBits = C->categoryBits | G->categoryBits;
			C->enlarged = B->enlarged || F->enlarged;
			A->enlarged = C->enlarged || G->enlarged;
			break;

		case s2_rotateCD:
			A->child2 = iD;
			B->child1 = iC;

			C->parent = iB;
			D->parent = iA;

			B->aabb = aabbCE;
			B->height = 1 + S2_MAX(C->height, E->height);
			A->height = 1 + S2_MAX(B->height, D->height);
			B->categoryBits = C->categoryBits | E->categoryBits;
			A->categoryBits = B->categoryBits | D->categoryBits;
			B->enlarged = C->enlarged || E->enlarged;
			A->enlarged = B->enlarged || D->enlarged;
			break;

		case s2_rotateCE:
			A->child2 = iE;
			B->child2 = iC;

			C->parent = iB;
			E->parent = iA;

			B->aabb = aabbCD;
			B->height = 1 + S2_MAX(C->height, D->height);
			A->height = 1 + S2_MAX(B->height, E->height);
			B->categoryBits = C->categoryBits | D->categoryBits;
			A->categoryBits = B->categoryBits | E->categoryBits;
			B->enlarged = C->enlarged || D->enlarged;
			A->enlarged = B->enlarged || E->enlarged;
			break;

		default:
			S2_ASSERT(false);
			break;
		}
	}
}

static void s2InsertLeaf(s2DynamicTree* tree, int32_t leaf, bool shouldRotate)
{
	if (tree->root == S2_NULL_INDEX)
	{
		tree->root = leaf;
		tree->nodes[tree->root].parent = S2_NULL_INDEX;
		return;
	}

	// Stage 1: find the best sibling for this node
	s2Box leafAABB = tree->nodes[leaf].aabb;
	int32_t sibling = s2FindBestSibling(tree, leafAABB);

	// Stage 2: create a new parent for the leaf and sibling
	int32_t oldParent = tree->nodes[sibling].parent;
	int32_t newParent = s2AllocateNode(tree);

	// warning: node pointer can change after allocation
	s2TreeNode* nodes = tree->nodes;
	nodes[newParent].parent = oldParent;
	nodes[newParent].userData = -1;
	nodes[newParent].aabb = s2AABB_Union(leafAABB, nodes[sibling].aabb);
	nodes[newParent].categoryBits = nodes[leaf].categoryBits | nodes[sibling].categoryBits;
	nodes[newParent].height = nodes[sibling].height + 1;

	if (oldParent != S2_NULL_INDEX)
	{
		// The sibling was not the root.
		if (nodes[oldParent].child1 == sibling)
		{
			nodes[oldParent].child1 = newParent;
		}
		else
		{
			nodes[oldParent].child2 = newParent;
		}

		nodes[newParent].child1 = sibling;
		nodes[newParent].child2 = leaf;
		nodes[sibling].parent = newParent;
		nodes[leaf].parent = newParent;
	}
	else
	{
		// The sibling was the root.
		nodes[newParent].child1 = sibling;
		nodes[newParent].child2 = leaf;
		nodes[sibling].parent = newParent;
		nodes[leaf].parent = newParent;
		tree->root = newParent;
	}

	// Stage 3: walk back up the tree fixing heights and AABBs
	int32_t index = nodes[leaf].parent;
	while (index != S2_NULL_INDEX)
	{
		int32_t child1 = nodes[index].child1;
		int32_t child2 = nodes[index].child2;

		S2_ASSERT(child1 != S2_NULL_INDEX);
		S2_ASSERT(child2 != S2_NULL_INDEX);

		nodes[index].aabb = s2AABB_Union(nodes[child1].aabb, nodes[child2].aabb);
		nodes[index].categoryBits = nodes[child1].categoryBits | nodes[child2].categoryBits;
		nodes[index].height = 1 + S2_MAX(nodes[child1].height, nodes[child2].height);
		nodes[index].enlarged = nodes[child1].enlarged || nodes[child2].enlarged;

		if (shouldRotate)
		{
			s2RotateNodes(tree, index);
		}

		index = nodes[index].parent;
	}
}

static void s2RemoveLeaf(s2DynamicTree* tree, int32_t leaf)
{
	if (leaf == tree->root)
	{
		tree->root = S2_NULL_INDEX;
		return;
	}

	s2TreeNode* nodes = tree->nodes;

	int32_t parent = nodes[leaf].parent;
	int32_t grandParent = nodes[parent].parent;
	int32_t sibling;
	if (nodes[parent].child1 == leaf)
	{
		sibling = nodes[parent].child2;
	}
	else
	{
		sibling = nodes[parent].child1;
	}

	if (grandParent != S2_NULL_INDEX)
	{
		// Destroy parent and connect sibling to grandParent.
		if (nodes[grandParent].child1 == parent)
		{
			nodes[grandParent].child1 = sibling;
		}
		else
		{
			nodes[grandParent].child2 = sibling;
		}
		nodes[sibling].parent = grandParent;
		s2FreeNode(tree, parent);

		// Adjust ancestor bounds.
		int32_t index = grandParent;
		while (index != S2_NULL_INDEX)
		{
			s2TreeNode* node = nodes + index;
			s2TreeNode* child1 = nodes + node->child1;
			s2TreeNode* child2 = nodes + node->child2;

			// Fast union using SSE
			//__m128 aabb1 = _mm_load_ps(&child1->aabb.lowerBound.x);
			//__m128 aabs2 = _mm_load_ps(&child2->aabb.lowerBound.x);
			//__m128 lower = _mm_min_ps(aabb1, aabs2);
			//__m128 upper = _mm_max_ps(aabb1, aabs2);
			//__m128 aabb = _mm_shuffle_ps(lower, upper, _MM_SHUFFLE(3, 2, 1, 0));
			//_mm_store_ps(&node->aabb.lowerBound.x, aabb);

			node->aabb = s2AABB_Union(child1->aabb, child2->aabb);
			node->categoryBits = child1->categoryBits | child2->categoryBits;
			node->height = 1 + S2_MAX(child1->height, child2->height);

			index = node->parent;
		}
	}
	else
	{
		tree->root = sibling;
		tree->nodes[sibling].parent = S2_NULL_INDEX;
		s2FreeNode(tree, parent);
	}
}

// Create a proxy in the tree as a leaf node. We return the index of the node instead of a pointer so that we can grow
// the node pool.
int32_t s2DynamicTree_CreateProxy(s2DynamicTree* tree, s2Box aabb, uint32_t categoryBits, int32_t userData)
{
	int32_t proxyId = s2AllocateNode(tree);
	s2TreeNode* node = tree->nodes + proxyId;

	node->aabb = aabb;
	node->userData = userData;
	node->categoryBits = categoryBits;
	node->height = 0;

	bool shouldRotate = true;
	s2InsertLeaf(tree, proxyId, shouldRotate);

	tree->proxyCount += 1;

	return proxyId;
}

void s2DynamicTree_DestroyProxy(s2DynamicTree* tree, int32_t proxyId)
{
	S2_ASSERT(0 <= proxyId && proxyId < tree->nodeCapacity);
	S2_ASSERT(s2IsLeaf(tree->nodes + proxyId));

	s2RemoveLeaf(tree, proxyId);
	s2FreeNode(tree, proxyId);

	S2_ASSERT(tree->proxyCount > 0);
	tree->proxyCount -= 1;
}

int32_t s2DynamicTree_GetProxyCount(const s2DynamicTree* tree)
{
	return tree->proxyCount;
}

void s2DynamicTree_MoveProxy(s2DynamicTree* tree, int32_t proxyId, s2Box aabb)
{
	S2_ASSERT(0 <= proxyId && proxyId < tree->nodeCapacity);
	S2_ASSERT(s2IsLeaf(tree->nodes + proxyId));

	s2RemoveLeaf(tree, proxyId);

	tree->nodes[proxyId].aabb = aabb;

	bool shouldRotate = false;
	s2InsertLeaf(tree, proxyId, shouldRotate);
}

void s2DynamicTree_EnlargeProxy(s2DynamicTree* tree, int32_t proxyId, s2Box aabb)
{
	s2TreeNode* nodes = tree->nodes;

	S2_ASSERT(0 <= proxyId && proxyId < tree->nodeCapacity);
	S2_ASSERT(s2IsLeaf(tree->nodes + proxyId));

	// Caller must ensure this
	S2_ASSERT(s2AABB_Contains(nodes[proxyId].aabb, aabb) == false);

	nodes[proxyId].aabb = aabb;

	int32_t parentIndex = nodes[proxyId].parent;
	while (parentIndex != S2_NULL_INDEX)
	{
		bool changed = s2AABB_Enlarge(&nodes[parentIndex].aabb, aabb);
		nodes[parentIndex].enlarged = true;
		parentIndex = nodes[parentIndex].parent;

		if (changed == false)
		{
			break;
		}
	}

	while (parentIndex != S2_NULL_INDEX)
	{
		if (nodes[parentIndex].enlarged == true)
		{
			// early out because this ancestor was previously ascended and marked as enlarged
			break;
		}

		nodes[parentIndex].enlarged = true;
		parentIndex = nodes[parentIndex].parent;
	}
}

int32_t s2DynamicTree_GetHeight(const s2DynamicTree* tree)
{
	if (tree->root == S2_NULL_INDEX)
	{
		return 0;
	}

	return tree->nodes[tree->root].height;
}

float s2DynamicTree_GetAreaRatio(const s2DynamicTree* tree)
{
	if (tree->root == S2_NULL_INDEX)
	{
		return 0.0f;
	}

	const s2TreeNode* root = tree->nodes + tree->root;
	float rootArea = s2AABB_Perimeter(root->aabb);

	float totalArea = 0.0f;
	for (int32_t i = 0; i < tree->nodeCapacity; ++i)
	{
		const s2TreeNode* node = tree->nodes + i;
		if (node->height < 0 || s2IsLeaf(node) || i == tree->root)
		{
			// Free node in pool
			continue;
		}

		totalArea += s2AABB_Perimeter(node->aabb);
	}

	return totalArea / rootArea;
}

// Compute the height of a sub-tree.
static int32_t s2ComputeHeight(const s2DynamicTree* tree, int32_t nodeId)
{
	S2_ASSERT(0 <= nodeId && nodeId < tree->nodeCapacity);
	s2TreeNode* node = tree->nodes + nodeId;

	if (s2IsLeaf(node))
	{
		return 0;
	}

	int32_t height1 = s2ComputeHeight(tree, node->child1);
	int32_t height2 = s2ComputeHeight(tree, node->child2);
	return 1 + S2_MAX(height1, height2);
}

int32_t s2DynamicTree_ComputeHeight(const s2DynamicTree* tree)
{
	int32_t height = s2ComputeHeight(tree, tree->root);
	return height;
}

#if S2_VALIDATE
static void s2ValidateStructure(const s2DynamicTree* tree, int32_t index)
{
	if (index == S2_NULL_INDEX)
	{
		return;
	}

	if (index == tree->root)
	{
		S2_ASSERT(tree->nodes[index].parent == S2_NULL_INDEX);
	}

	const s2TreeNode* node = tree->nodes + index;

	int32_t child1 = node->child1;
	int32_t child2 = node->child2;

	if (s2IsLeaf(node))
	{
		S2_ASSERT(child1 == S2_NULL_INDEX);
		S2_ASSERT(child2 == S2_NULL_INDEX);
		S2_ASSERT(node->height == 0);
		return;
	}

	S2_ASSERT(0 <= child1 && child1 < tree->nodeCapacity);
	S2_ASSERT(0 <= child2 && child2 < tree->nodeCapacity);

	S2_ASSERT(tree->nodes[child1].parent == index);
	S2_ASSERT(tree->nodes[child2].parent == index);

	if (tree->nodes[child1].enlarged || tree->nodes[child2].enlarged)
	{
		S2_ASSERT(node->enlarged == true);
	}

	s2ValidateStructure(tree, child1);
	s2ValidateStructure(tree, child2);
}

static void s2ValidateMetrics(const s2DynamicTree* tree, int32_t index)
{
	if (index == S2_NULL_INDEX)
	{
		return;
	}

	const s2TreeNode* node = tree->nodes + index;

	int32_t child1 = node->child1;
	int32_t child2 = node->child2;

	if (s2IsLeaf(node))
	{
		S2_ASSERT(child1 == S2_NULL_INDEX);
		S2_ASSERT(child2 == S2_NULL_INDEX);
		S2_ASSERT(node->height == 0);
		return;
	}

	S2_ASSERT(0 <= child1 && child1 < tree->nodeCapacity);
	S2_ASSERT(0 <= child2 && child2 < tree->nodeCapacity);

	int32_t height1 = tree->nodes[child1].height;
	int32_t height2 = tree->nodes[child2].height;
	int32_t height;
	height = 1 + S2_MAX(height1, height2);
	S2_ASSERT(node->height == height);

	// s2Box aabb = s2AABB_Union(tree->nodes[child1].aabb, tree->nodes[child2].aabb);

	S2_ASSERT(s2AABB_Contains(node->aabb, tree->nodes[child1].aabb));
	S2_ASSERT(s2AABB_Contains(node->aabb, tree->nodes[child2].aabb));

	// S2_ASSERT(aabb.lowerBound.x == node->aabb.lowerBound.x);
	// S2_ASSERT(aabb.lowerBound.y == node->aabb.lowerBound.y);
	// S2_ASSERT(aabb.upperBound.x == node->aabb.upperBound.x);
	// S2_ASSERT(aabb.upperBound.y == node->aabb.upperBound.y);

	uint32_t categoryBits = tree->nodes[child1].categoryBits | tree->nodes[child2].categoryBits;
	S2_ASSERT(node->categoryBits == categoryBits);

	s2ValidateMetrics(tree, child1);
	s2ValidateMetrics(tree, child2);
}
#endif

void s2DynamicTree_Validate(const s2DynamicTree* tree)
{
#if S2_VALIDATE
	if (tree->root == S2_NULL_INDEX)
	{
		return;
	}

	s2ValidateStructure(tree, tree->root);
	s2ValidateMetrics(tree, tree->root);

	int32_t freeCount = 0;
	int32_t freeIndex = tree->freeList;
	while (freeIndex != S2_NULL_INDEX)
	{
		S2_ASSERT(0 <= freeIndex && freeIndex < tree->nodeCapacity);
		freeIndex = tree->nodes[freeIndex].next;
		++freeCount;
	}

	int32_t height = s2DynamicTree_GetHeight(tree);
	int32_t computedHeight = s2DynamicTree_ComputeHeight(tree);
	S2_ASSERT(height == computedHeight);

	S2_ASSERT(tree->nodeCount + freeCount == tree->nodeCapacity);
#else
	S2_MAYBE_UNUSED(tree);
#endif
}

int32_t s2DynamicTree_GetMaxBalance(const s2DynamicTree* tree)
{
	int32_t maxBalance = 0;
	for (int32_t i = 0; i < tree->nodeCapacity; ++i)
	{
		const s2TreeNode* node = tree->nodes + i;
		if (node->height <= 1)
		{
			continue;
		}

		S2_ASSERT(s2IsLeaf(node) == false);

		int32_t child1 = node->child1;
		int32_t child2 = node->child2;
		int32_t balance = S2_ABS(tree->nodes[child2].height - tree->nodes[child1].height);
		maxBalance = S2_MAX(maxBalance, balance);
	}

	return maxBalance;
}

void s2DynamicTree_RebuildBottomUp(s2DynamicTree* tree)
{
	int32_t* nodes = (int32_t*)s2Alloc(tree->nodeCount * sizeof(int32_t));
	int32_t count = 0;

	// Build array of leaves. Free the rest.
	for (int32_t i = 0; i < tree->nodeCapacity; ++i)
	{
		if (tree->nodes[i].height < 0)
		{
			// free node in pool
			continue;
		}

		if (s2IsLeaf(tree->nodes + i))
		{
			tree->nodes[i].parent = S2_NULL_INDEX;
			nodes[count] = i;
			++count;
		}
		else
		{
			s2FreeNode(tree, i);
		}
	}

	while (count > 1)
	{
		float minCost = FLT_MAX;
		int32_t iMin = -1, jMin = -1;
		for (int32_t i = 0; i < count; ++i)
		{
			s2Box aabbi = tree->nodes[nodes[i]].aabb;

			for (int32_t j = i + 1; j < count; ++j)
			{
				s2Box aabbj = tree->nodes[nodes[j]].aabb;
				s2Box b = s2AABB_Union(aabbi, aabbj);
				float cost = s2AABB_Perimeter(b);
				if (cost < minCost)
				{
					iMin = i;
					jMin = j;
					minCost = cost;
				}
			}
		}

		int32_t index1 = nodes[iMin];
		int32_t index2 = nodes[jMin];
		s2TreeNode* child1 = tree->nodes + index1;
		s2TreeNode* child2 = tree->nodes + index2;

		int32_t parentIndex = s2AllocateNode(tree);
		s2TreeNode* parent = tree->nodes + parentIndex;
		parent->child1 = index1;
		parent->child2 = index2;
		parent->aabb = s2AABB_Union(child1->aabb, child2->aabb);
		parent->categoryBits = child1->categoryBits | child2->categoryBits;
		parent->height = 1 + S2_MAX(child1->height, child2->height);
		parent->parent = S2_NULL_INDEX;

		child1->parent = parentIndex;
		child2->parent = parentIndex;

		nodes[jMin] = nodes[count - 1];
		nodes[iMin] = parentIndex;
		--count;
	}

	tree->root = nodes[0];
	s2Free(nodes, tree->nodeCount * sizeof(s2TreeNode));

	s2DynamicTree_Validate(tree);
}

void s2DynamicTree_ShiftOrigin(s2DynamicTree* tree, s2Vec2 newOrigin)
{
	// Build array of leaves. Free the rest.
	for (int32_t i = 0; i < tree->nodeCapacity; ++i)
	{
		s2TreeNode* n = tree->nodes + i;
		n->aabb.lowerBound.x -= newOrigin.x;
		n->aabb.lowerBound.y -= newOrigin.y;
		n->aabb.upperBound.x -= newOrigin.x;
		n->aabb.upperBound.y -= newOrigin.y;
	}
}

#define s2_treeStackSize 256

void s2DynamicTree_QueryFiltered(const s2DynamicTree* tree, s2Box aabb, uint32_t maskBits, s2TreeQueryCallbackFcn* callback, void* context)
{
	int32_t stack[s2_treeStackSize];
	int32_t stackCount = 0;
	stack[stackCount++] = tree->root;

	while (stackCount > 0)
	{
		int32_t nodeId = stack[--stackCount];
		if (nodeId == S2_NULL_INDEX)
		{
			continue;
		}

		const s2TreeNode* node = tree->nodes + nodeId;

		if (s2AABB_Overlaps(node->aabb, aabb) && (node->categoryBits & maskBits) != 0)
		{
			if (s2IsLeaf(node))
			{
				// callback to user code with proxy id
				bool proceed = callback(nodeId, node->userData, context);
				if (proceed == false)
				{
					return;
				}
			}
			else
			{
				S2_ASSERT(stackCount <= s2_treeStackSize - 2);
				// TODO log this?

				if (stackCount <= s2_treeStackSize - 2)
				{
					stack[stackCount++] = node->child1;
					stack[stackCount++] = node->child2;
				}
			}
		}
	}
}

void s2DynamicTree_Query(const s2DynamicTree* tree, s2Box aabb, s2TreeQueryCallbackFcn* callback, void* context)
{
	int32_t stack[s2_treeStackSize];
	int32_t stackCount = 0;
	stack[stackCount++] = tree->root;

	while (stackCount > 0)
	{
		int32_t nodeId = stack[--stackCount];
		if (nodeId == S2_NULL_INDEX)
		{
			continue;
		}

		const s2TreeNode* node = tree->nodes + nodeId;

		if (s2AABB_Overlaps(node->aabb, aabb))
		{
			if (s2IsLeaf(node))
			{
				// callback to user code with proxy id
				bool proceed = callback(nodeId, node->userData, context);
				if (proceed == false)
				{
					return;
				}
			}
			else
			{
				S2_ASSERT(stackCount <= s2_treeStackSize - 2);
				// TODO log this?

				if (stackCount <= s2_treeStackSize - 2)
				{
					stack[stackCount++] = node->child1;
					stack[stackCount++] = node->child2;
				}
			}
		}
	}
}

void s2DynamicTree_RayCast(const s2DynamicTree* tree, const s2RayCastInput* input, uint32_t maskBits, s2TreeRayCastCallbackFcn* callback,
						   void* context)
{
	s2Vec2 p1 = input->p1;
	s2Vec2 p2 = input->p2;

	s2Vec2 r = s2Normalize(s2Sub(p2, p1));

	// v is perpendicular to the segment.
	s2Vec2 v = s2CrossSV(1.0f, r);
	s2Vec2 abs_v = s2Abs(v);

	// Separating axis for segment (Gino, p80).
	// |dot(v, p1 - c)| > dot(|v|, h)

	float maxFraction = input->maxFraction;

	// Build a bounding box for the segment.
	s2Box segmentAABB;
	{
		// t is the endpoint of the ray
		s2Vec2 t = s2MulAdd(p1, maxFraction, s2Sub(p2, p1));

		segmentAABB.lowerBound = s2Min(p1, t);
		segmentAABB.upperBound = s2Max(p1, t);
	}

	int32_t stack[s2_treeStackSize];
	int32_t stackCount = 0;
	stack[stackCount++] = tree->root;

	while (stackCount > 0)
	{
		int32_t nodeId = stack[--stackCount];
		if (nodeId == S2_NULL_INDEX)
		{
			continue;
		}

		const s2TreeNode* node = tree->nodes + nodeId;
		if (s2AABB_Overlaps(node->aabb, segmentAABB) == false || (node->categoryBits & maskBits) == 0)
		{
			continue;
		}

		// Separating axis for segment (Gino, p80).
		// |dot(v, p1 - c)| > dot(|v|, h)
		// radius extension is added to the node in this case
		s2Vec2 c = s2AABB_Center(node->aabb);
		s2Vec2 h = s2AABB_Extents(node->aabb);
		float term1 = S2_ABS(s2Dot(v, s2Sub(p1, c)));
		float term2 = s2Dot(abs_v, h);
		if (term2 < term1)
		{
			continue;
		}

		if (s2IsLeaf(node))
		{
			s2RayCastInput subInput;
			subInput.p1 = input->p1;
			subInput.p2 = input->p2;
			subInput.maxFraction = maxFraction;

			float value = callback(&subInput, nodeId, node->userData, context);
			S2_ASSERT(value >= 0.0f);

			if (value == 0.0f)
			{
				// The client has terminated the ray cast.
				return;
			}

			if (value < maxFraction)
			{
				// Update segment bounding box.
				maxFraction = value;
				s2Vec2 t = s2MulAdd(p1, maxFraction, s2Sub(p2, p1));
				segmentAABB.lowerBound = s2Min(p1, t);
				segmentAABB.upperBound = s2Max(p1, t);
			}
		}
		else
		{
			S2_ASSERT(stackCount <= s2_treeStackSize - 2);
			// TODO log this?

			if (stackCount <= s2_treeStackSize - 2)
			{
				// TODO_ERIN just put one node on the stack, continue on a child node
				// TODO_ERIN test ordering children by nearest to ray origin
				stack[stackCount++] = node->child1;
				stack[stackCount++] = node->child2;
			}
		}
	}
}

// Median split == 0, Surface area heurstic == 1
#define S2_TREE_HEURISTIC 0

#if S2_TREE_HEURISTIC == 0

// Median split heuristic
static int32_t s2PartitionMid(int32_t* indices, s2Vec2* centers, int32_t count)
{
	// Handle trivial case
	if (count <= 2)
	{
		return count / 2;
	}

	// TODO_ERIN SIMD?
	s2Vec2 lowerBound = centers[0];
	s2Vec2 upperBound = centers[0];

	for (int32_t i = 1; i < count; ++i)
	{
		lowerBound = s2Min(lowerBound, centers[i]);
		upperBound = s2Max(upperBound, centers[i]);
	}

	s2Vec2 d = s2Sub(upperBound, lowerBound);
	s2Vec2 c = {0.5f * (lowerBound.x + upperBound.x), 0.5f * (lowerBound.y + upperBound.y)};

	// Partition longest axis using the Hoare partition scheme
	// https://en.wikipedia.org/wiki/Quicksort
	// https://nicholasvadivelu.com/2021/01/11/array-partition/
	int32_t i1 = 0, i2 = count;
	if (d.x > d.y)
	{
		float pivot = c.x;

		while (i1 < i2)
		{
			while (i1 < i2 && centers[i1].x < pivot)
			{
				i1 += 1;
			};

			while (i1 < i2 && centers[i2 - 1].x >= pivot)
			{
				i2 -= 1;
			};

			if (i1 < i2)
			{
				// Swap indices
				{
					int32_t temp = indices[i1];
					indices[i1] = indices[i2 - 1];
					indices[i2 - 1] = temp;
				}

				// Swap centers
				{
					s2Vec2 temp = centers[i1];
					centers[i1] = centers[i2 - 1];
					centers[i2 - 1] = temp;
				}

				i1 += 1;
				i2 -= 1;
			}
		}
	}
	else
	{
		float pivot = c.y;

		while (i1 < i2)
		{
			while (i1 < i2 && centers[i1].y < pivot)
			{
				i1 += 1;
			};

			while (i1 < i2 && centers[i2 - 1].y >= pivot)
			{
				i2 -= 1;
			};

			if (i1 < i2)
			{
				// Swap indices
				{
					int32_t temp = indices[i1];
					indices[i1] = indices[i2 - 1];
					indices[i2 - 1] = temp;
				}

				// Swap centers
				{
					s2Vec2 temp = centers[i1];
					centers[i1] = centers[i2 - 1];
					centers[i2 - 1] = temp;
				}

				i1 += 1;
				i2 -= 1;
			}
		}
	}
	S2_ASSERT(i1 == i2);

	if (i1 > 0 && i1 < count)
	{
		return i1;
	}
	else
	{
		return count / 2;
	}
}

#else

#define S2_BIN_COUNT 8

typedef struct s2TreeBin
{
	s2Box aabb;
	int32_t count;
} s2TreeBin;

typedef struct s2TreePlane
{
	s2Box leftAABB;
	s2Box rightAABB;
	int32_t leftCount;
	int32_t rightCount;
} s2TreePlane;

// "On Fast Construction of SAH-based Bounding Volume Hierarchies" by Ingo Wald
// Returns the left child count
static int32_t s2PartitionSAH(int32_t* indices, int32_t* binIndices, s2Box* boxes, int32_t count)
{
	S2_ASSERT(count > 0);

	s2TreeBin bins[S2_BIN_COUNT];
	s2TreePlane planes[S2_BIN_COUNT - 1];

	s2Vec2 center = s2AABB_Center(boxes[0]);
	s2Box centroidAABB;
	centroidAABB.lowerBound = center;
	centroidAABB.upperBound = center;

	for (int32_t i = 1; i < count; ++i)
	{
		center = s2AABB_Center(boxes[i]);
		centroidAABB.lowerBound = s2Min(centroidAABB.lowerBound, center);
		centroidAABB.upperBound = s2Max(centroidAABB.upperBound, center);
	}

	s2Vec2 d = s2Sub(centroidAABB.upperBound, centroidAABB.lowerBound);

	// Find longest axis
	int32_t axisIndex;
	float invD;
	if (d.x > d.y)
	{
		axisIndex = 0;
		invD = d.x;
	}
	else
	{
		axisIndex = 1;
		invD = d.y;
	}

	invD = invD > 0.0f ? 1.0f / invD : 0.0f;

	// Initialize bin bounds and count
	for (int32_t i = 0; i < S2_BIN_COUNT; ++i)
	{
		bins[i].aabb.lowerBound = (s2Vec2){FLT_MAX, FLT_MAX};
		bins[i].aabb.upperBound = (s2Vec2){-FLT_MAX, -FLT_MAX};
		bins[i].count = 0;
	}

	// Assign boxes to bins and compute bin boxes
	// TODO_ERIN optimize
	float binCount = S2_BIN_COUNT;
	float lowerBoundArray[2] = {centroidAABB.lowerBound.x, centroidAABB.lowerBound.y};
	float minC = lowerBoundArray[axisIndex];
	for (int32_t i = 0; i < count; ++i)
	{
		s2Vec2 c = s2AABB_Center(boxes[i]);
		float cArray[2] = {c.x, c.y};
		int32_t binIndex = (int32_t)(binCount * (cArray[axisIndex] - minC) * invD);
		binIndex = S2_CLAMP(binIndex, 0, S2_BIN_COUNT - 1);
		binIndices[i] = binIndex;
		bins[binIndex].count += 1;
		bins[binIndex].aabb = s2AABB_Union(bins[binIndex].aabb, boxes[i]);
	}

	int32_t planeCount = S2_BIN_COUNT - 1;

	// Prepare all the left planes, candidates for left child
	planes[0].leftCount = bins[0].count;
	planes[0].leftAABB = bins[0].aabb;
	for (int32_t i = 1; i < planeCount; ++i)
	{
		planes[i].leftCount = planes[i - 1].leftCount + bins[i].count;
		planes[i].leftAABB = s2AABB_Union(planes[i - 1].leftAABB, bins[i].aabb);
	}

	// Prepare all the right planes, candidates for right child
	planes[planeCount - 1].rightCount = bins[planeCount].count;
	planes[planeCount - 1].rightAABB = bins[planeCount].aabb;
	for (int32_t i = planeCount - 2; i >= 0; --i)
	{
		planes[i].rightCount = planes[i + 1].rightCount + bins[i + 1].count;
		planes[i].rightAABB = s2AABB_Union(planes[i + 1].rightAABB, bins[i + 1].aabb);
	}

	// Find best split to minimize SAH
	float minCost = FLT_MAX;
	int32_t bestPlane = 0;
	for (int32_t i = 0; i < planeCount; ++i)
	{
		float leftArea = s2AABB_Perimeter(planes[i].leftAABB);
		float rightArea = s2AABB_Perimeter(planes[i].rightAABB);
		int32_t leftCount = planes[i].leftCount;
		int32_t rightCount = planes[i].rightCount;

		float cost = leftCount * leftArea + rightCount * rightArea;
		if (cost < minCost)
		{
			bestPlane = i;
			minCost = cost;
		}
	}

	// Partition node indices and boxes using the Hoare partition scheme
	// https://en.wikipedia.org/wiki/Quicksort
	// https://nicholasvadivelu.com/2021/01/11/array-partition/
	int32_t i1 = 0, i2 = count;
	while (i1 < i2)
	{
		while (i1 < i2 && binIndices[i1] < bestPlane)
		{
			i1 += 1;
		};

		while (i1 < i2 && binIndices[i2 - 1] >= bestPlane)
		{
			i2 -= 1;
		};

		if (i1 < i2)
		{
			// Swap indices
			{
				int32_t temp = indices[i1];
				indices[i1] = indices[i2 - 1];
				indices[i2 - 1] = temp;
			}

			// Swap boxes
			{
				s2Box temp = boxes[i1];
				boxes[i1] = boxes[i2 - 1];
				boxes[i2 - 1] = temp;
			}

			i1 += 1;
			i2 -= 1;
		}
	}
	S2_ASSERT(i1 == i2);

	if (i1 > 0 && i1 < count)
	{
		return i1;
	}
	else
	{
		return count / 2;
	}
}

#endif

// Temporary data used to track the rebuild of a tree node
struct s2RebuildItem
{
	int32_t nodeIndex;
	int32_t childCount;

	// Leaf indices
	int32_t startIndex;
	int32_t splitIndex;
	int32_t endIndex;
};

// Returns root node index
static int32_t s2BuildTree(s2DynamicTree* tree, int32_t leafCount)
{
	s2TreeNode* nodes = tree->nodes;
	int32_t* leafIndices = tree->leafIndices;

	if (leafCount == 1)
	{
		nodes[leafIndices[0]].parent = S2_NULL_INDEX;
		return leafIndices[0];
	}

#if S2_TREE_HEURISTIC == 0
	s2Vec2* leafCenters = tree->leafCenters;
#else
	s2Box* leafBoxes = tree->leafBoxes;
	int32_t* binIndices = tree->binIndices;
#endif

	struct s2RebuildItem stack[s2_treeStackSize];
	int32_t top = 0;

	stack[0].nodeIndex = s2AllocateNode(tree);
	stack[0].childCount = -1;
	stack[0].startIndex = 0;
	stack[0].endIndex = leafCount;
#if S2_TREE_HEURISTIC == 0
	stack[0].splitIndex = s2PartitionMid(leafIndices, leafCenters, leafCount);
#else
	stack[0].splitIndex = s2PartitionSAH(leafIndices, binIndices, leafBoxes, leafCount);
#endif

	while (true)
	{
		struct s2RebuildItem* item = stack + top;

		item->childCount += 1;

		if (item->childCount == 2)
		{
			// This internal node has both children established

			if (top == 0)
			{
				// all done
				break;
			}

			struct s2RebuildItem* parentItem = stack + (top - 1);
			s2TreeNode* parentNode = nodes + parentItem->nodeIndex;

			if (parentItem->childCount == 0)
			{
				S2_ASSERT(parentNode->child1 == S2_NULL_INDEX);
				parentNode->child1 = item->nodeIndex;
			}
			else
			{
				S2_ASSERT(parentItem->childCount == 1);
				S2_ASSERT(parentNode->child2 == S2_NULL_INDEX);
				parentNode->child2 = item->nodeIndex;
			}

			s2TreeNode* node = nodes + item->nodeIndex;

			S2_ASSERT(node->parent == S2_NULL_INDEX);
			node->parent = parentItem->nodeIndex;

			S2_ASSERT(node->child1 != S2_NULL_INDEX);
			S2_ASSERT(node->child2 != S2_NULL_INDEX);
			s2TreeNode* child1 = nodes + node->child1;
			s2TreeNode* child2 = nodes + node->child2;

			node->aabb = s2AABB_Union(child1->aabb, child2->aabb);
			node->height = 1 + S2_MAX(child1->height, child2->height);
			node->categoryBits = child1->categoryBits | child2->categoryBits;

			// Pop stack
			top -= 1;
		}
		else
		{
			int32_t startIndex, endIndex;
			if (item->childCount == 0)
			{
				startIndex = item->startIndex;
				endIndex = item->splitIndex;
			}
			else
			{
				S2_ASSERT(item->childCount == 1);
				startIndex = item->splitIndex;
				endIndex = item->endIndex;
			}

			int32_t count = endIndex - startIndex;

			if (count == 1)
			{
				int32_t childIndex = leafIndices[startIndex];
				s2TreeNode* node = nodes + item->nodeIndex;

				if (item->childCount == 0)
				{
					S2_ASSERT(node->child1 == S2_NULL_INDEX);
					node->child1 = childIndex;
				}
				else
				{
					S2_ASSERT(item->childCount == 1);
					S2_ASSERT(node->child2 == S2_NULL_INDEX);
					node->child2 = childIndex;
				}

				s2TreeNode* childNode = nodes + childIndex;
				S2_ASSERT(childNode->parent == S2_NULL_INDEX);
				childNode->parent = item->nodeIndex;
			}
			else
			{
				S2_ASSERT(count > 0);
				S2_ASSERT(top < s2_treeStackSize);

				top += 1;
				struct s2RebuildItem* newItem = stack + top;
				newItem->nodeIndex = s2AllocateNode(tree);
				newItem->childCount = -1;
				newItem->startIndex = startIndex;
				newItem->endIndex = endIndex;
#if S2_TREE_HEURISTIC == 0
				newItem->splitIndex = s2PartitionMid(leafIndices + startIndex, leafCenters + startIndex, count);
#else
				newItem->splitIndex = s2PartitionSAH(leafIndices + startIndex, binIndices + startIndex, leafBoxes + startIndex, count);
#endif
				newItem->splitIndex += startIndex;
			}
		}
	}

	s2TreeNode* rootNode = nodes + stack[0].nodeIndex;
	S2_ASSERT(rootNode->parent == S2_NULL_INDEX);
	S2_ASSERT(rootNode->child1 != S2_NULL_INDEX);
	S2_ASSERT(rootNode->child2 != S2_NULL_INDEX);

	s2TreeNode* child1 = nodes + rootNode->child1;
	s2TreeNode* child2 = nodes + rootNode->child2;

	rootNode->aabb = s2AABB_Union(child1->aabb, child2->aabb);
	rootNode->height = 1 + S2_MAX(child1->height, child2->height);
	rootNode->categoryBits = child1->categoryBits | child2->categoryBits;

	return stack[0].nodeIndex;
}

// Not safe to access tree during this operation because it may grow
int32_t s2DynamicTree_Rebuild(s2DynamicTree* tree, bool fullBuild)
{
	int32_t proxyCount = tree->proxyCount;
	if (proxyCount == 0)
	{
		return 0;
	}

	// Ensure capacity for rebuild space
	if (proxyCount > tree->rebuildCapacity)
	{
		int32_t newCapacity = proxyCount + proxyCount / 2;

		s2Free(tree->leafIndices, tree->rebuildCapacity * sizeof(int32_t));
		tree->leafIndices = s2Alloc(newCapacity * sizeof(int32_t));

#if S2_TREE_HEURISTIC == 0
		s2Free(tree->leafCenters, tree->rebuildCapacity * sizeof(s2Vec2));
		tree->leafCenters = s2Alloc(newCapacity * sizeof(s2Vec2));
#else
		s2Free(tree->leafBoxes, tree->rebuildCapacity * sizeof(s2Box));
		tree->leafBoxes = s2Alloc(newCapacity * sizeof(s2Box));
		s2Free(tree->binIndices, tree->rebuildCapacity * sizeof(int32_t));
		tree->binIndices = s2Alloc(newCapacity * sizeof(int32_t));
#endif
		tree->rebuildCapacity = newCapacity;
	}

	int32_t leafCount = 0;
	int32_t stack[s2_treeStackSize];
	int32_t stackCount = 0;

	int32_t nodeIndex = tree->root;
	s2TreeNode* nodes = tree->nodes;
	s2TreeNode* node = nodes + nodeIndex;

	// These are the nodes that get sorted to rebuild the tree.
	// I'm using indices because the node pool may grow during the build.
	int32_t* leafIndices = tree->leafIndices;

#if S2_TREE_HEURISTIC == 0
	s2Vec2* leafCenters = tree->leafCenters;
#else
	s2Box* leafBoxes = tree->leafBoxes;
#endif

	// Gather all proxy nodes that have grown and all internal nodes that haven't grown. Both are
	// considered leaves in the tree rebuild.
	// Free all internal nodes that have grown.
	while (true)
	{
		if (node->height == 0 || (node->enlarged == false && fullBuild == false))
		{
			leafIndices[leafCount] = nodeIndex;
#if S2_TREE_HEURISTIC == 0
			leafCenters[leafCount] = s2AABB_Center(node->aabb);
#else
			leafBoxes[leafCount] = node->aabb;
#endif
			leafCount += 1;

			// Detach
			node->parent = S2_NULL_INDEX;
		}
		else
		{
			int32_t doomedNodeIndex = nodeIndex;

			// Handle children
			nodeIndex = node->child1;

			S2_ASSERT(stackCount < s2_treeStackSize);
			if (stackCount < s2_treeStackSize)
			{
				stack[stackCount++] = node->child2;
			}

			node = nodes + nodeIndex;

			// Remove doomed node
			s2FreeNode(tree, doomedNodeIndex);

			continue;
		}

		if (stackCount == 0)
		{
			break;
		}

		nodeIndex = stack[--stackCount];
		node = nodes + nodeIndex;
	}

	#if S2_VALIDATE == 1
	int32_t capacity = tree->nodeCapacity;
	for (int32_t i = 0; i < capacity; ++i)
	{
		if (nodes[i].height >= 0)
		{
			S2_ASSERT(nodes[i].enlarged == false);
		}
	}
	#endif

	S2_ASSERT(leafCount <= proxyCount);

	tree->root = s2BuildTree(tree, leafCount);

	s2DynamicTree_Validate(tree);

	return leafCount;
}
