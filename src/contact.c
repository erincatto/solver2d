// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "contact.h"

#include "array.h"
#include "block_allocator.h"
#include "body.h"
#include "core.h"
#include "shape.h"
#include "table.h"
#include "world.h"

#include "solver2d/distance.h"
#include "solver2d/manifold.h"
#include "solver2d/timer.h"

#include <float.h>
#include <math.h>

// Contacts and determinism
// A deterministic simulation requires contacts to exist in the same order in s2Island no matter the thread count.
// The order must reproduce from run to run. This is necessary because the Gauss-Seidel constraint solver is order dependent.
//
// Creation:
// - Contacts are created using results from s2UpdateBroadPhasePairs
// - These results are ordered according to the order of the broad-phase move array
// - The move array is ordered according to the shape creation order using a bitset.
// - The island/shape/body order is determined by creation order
// - Logically contacts are only created for awake bodies, so they are immediately added to the awake contact array (serially)
//
// Island linking:
// - The awake contact array is built from the body-contact graph for all awake bodies in awake islands.
// - Awake contacts are solved in parallel and they generate contact state changes.
// - These state changes may link islands together using union find.
// - The state changes are ordered using a bit array that encompasses all contacts
// - As long as contacts are created in deterministic order, island link order is deterministic.
// - This keeps the order of contacts in islands deterministic

// Friction mixing law. The idea is to allow either fixture to drive the friction to zero.
// For example, anything slides on ice.
static inline float s2MixFriction(float friction1, float friction2)
{
	return sqrtf(friction1 * friction2);
}

// Restitution mixing law. The idea is allow for anything to bounce off an inelastic surface.
// For example, a superball bounces on anything.
static inline float s2MixRestitution(float restitution1, float restitution2)
{
	return restitution1 > restitution2 ? restitution1 : restitution2;
}

typedef s2Manifold s2ManifoldFcn(const s2Shape* shapeA, s2Transform xfA, const s2Shape* shapeB, s2Transform xfB, float maxDistance,
								 s2DistanceCache* cache);

struct s2ContactRegister
{
	s2ManifoldFcn* fcn;
	bool primary;
};

static struct s2ContactRegister s_registers[s2_shapeTypeCount][s2_shapeTypeCount];
static bool s_initialized = false;

static s2Manifold s2CircleManifold(const s2Shape* shapeA, s2Transform xfA, const s2Shape* shapeB, s2Transform xfB, float maxDistance,
								   s2DistanceCache* cache)
{
	S2_MAYBE_UNUSED(cache);
	return s2CollideCircles(&shapeA->circle, xfA, &shapeB->circle, xfB, maxDistance);
}

static s2Manifold s2CapsuleAndCircleManifold(const s2Shape* shapeA, s2Transform xfA, const s2Shape* shapeB, s2Transform xfB, float maxDistance,
								   s2DistanceCache* cache)
{
	S2_MAYBE_UNUSED(cache);
	return s2CollideCapsuleAndCircle(&shapeA->capsule, xfA, &shapeB->circle, xfB, maxDistance);
}

static s2Manifold s2CapsuleManifold(const s2Shape* shapeA, s2Transform xfA, const s2Shape* shapeB, s2Transform xfB, float maxDistance,
								   s2DistanceCache* cache)
{
	return s2CollideCapsules(&shapeA->capsule, xfA, &shapeB->capsule, xfB, maxDistance, cache);
}

static s2Manifold s2PolygonAndCircleManifold(const s2Shape* shapeA, s2Transform xfA, const s2Shape* shapeB, s2Transform xfB,
											 float maxDistance, s2DistanceCache* cache)
{
	S2_MAYBE_UNUSED(cache);
	return s2CollidePolygonAndCircle(&shapeA->polygon, xfA, &shapeB->circle, xfB, maxDistance);
}

static s2Manifold s2PolygonAndCapsuleManifold(const s2Shape* shapeA, s2Transform xfA, const s2Shape* shapeB, s2Transform xfB,
											 float maxDistance, s2DistanceCache* cache)
{
	return s2CollidePolygonAndCapsule(&shapeA->polygon, xfA, &shapeB->capsule, xfB, maxDistance, cache);
}

static s2Manifold s2PolygonManifold(const s2Shape* shapeA, s2Transform xfA, const s2Shape* shapeB, s2Transform xfB, float maxDistance,
									s2DistanceCache* cache)
{
	return s2CollidePolygons(&shapeA->polygon, xfA, &shapeB->polygon, xfB, maxDistance, cache);
}

static s2Manifold s2SegmentAndCircleManifold(const s2Shape* shapeA, s2Transform xfA, const s2Shape* shapeB, s2Transform xfB, float maxDistance,
								   s2DistanceCache* cache)
{
	S2_MAYBE_UNUSED(cache);
	return s2CollideSegmentAndCircle(&shapeA->segment, xfA, &shapeB->circle, xfB, maxDistance);
}

static s2Manifold s2SegmentAndCapsuleManifold(const s2Shape* shapeA, s2Transform xfA, const s2Shape* shapeB, s2Transform xfB, float maxDistance,
								   s2DistanceCache* cache)
{
	return s2CollideSegmentAndCapsule(&shapeA->segment, xfA, &shapeB->capsule, xfB, maxDistance, cache);
}

static s2Manifold s2SegmentAndPolygonManifold(const s2Shape* shapeA, s2Transform xfA, const s2Shape* shapeB, s2Transform xfB, float maxDistance,
								   s2DistanceCache* cache)
{
	return s2CollideSegmentAndPolygon(&shapeA->segment, xfA, &shapeB->polygon, xfB, maxDistance, cache);
}

static void s2AddType(s2ManifoldFcn* fcn, enum s2ShapeType type1, enum s2ShapeType type2)
{
	S2_ASSERT(0 <= type1 && type1 < s2_shapeTypeCount);
	S2_ASSERT(0 <= type2 && type2 < s2_shapeTypeCount);

	s_registers[type1][type2].fcn = fcn;
	s_registers[type1][type2].primary = true;

	if (type1 != type2)
	{
		s_registers[type2][type1].fcn = fcn;
		s_registers[type2][type1].primary = false;
	}
}

void s2InitializeContactRegisters(void)
{
	if (s_initialized == false)
	{
		s2AddType(s2CircleManifold, s2_circleShape, s2_circleShape);
		s2AddType(s2CapsuleAndCircleManifold, s2_capsuleShape, s2_circleShape);
		s2AddType(s2CapsuleManifold, s2_capsuleShape, s2_capsuleShape);
		s2AddType(s2PolygonAndCircleManifold, s2_polygonShape, s2_circleShape);
		s2AddType(s2PolygonAndCapsuleManifold, s2_polygonShape, s2_capsuleShape);
		s2AddType(s2PolygonManifold, s2_polygonShape, s2_polygonShape);
		s2AddType(s2SegmentAndCircleManifold, s2_segmentShape, s2_circleShape);
		s2AddType(s2SegmentAndCapsuleManifold, s2_segmentShape, s2_capsuleShape);
		s2AddType(s2SegmentAndPolygonManifold, s2_segmentShape, s2_polygonShape);
		s_initialized = true;
	}
}

void s2CreateContact(s2World* world, s2Shape* shapeA, s2Shape* shapeB)
{
	s2ShapeType type1 = shapeA->type;
	s2ShapeType type2 = shapeB->type;

	S2_ASSERT(0 <= type1 && type1 < s2_shapeTypeCount);
	S2_ASSERT(0 <= type2 && type2 < s2_shapeTypeCount);

	if (s_registers[type1][type2].fcn == NULL)
	{
		// For example, no segment vs segment collision
		return;
	}

	if (s_registers[type1][type2].primary == false)
	{
		// flip order
		s2CreateContact(world, shapeB, shapeA);
		return;
	}

	s2Contact* contact = (s2Contact*)s2AllocObject(&world->contactPool);
	world->contacts = (s2Contact*)world->contactPool.memory;

	int32_t contactIndex = contact->object.index;

	contact->shapeIndexA = shapeA->object.index;
	contact->shapeIndexB = shapeB->object.index;
	contact->cache = s2_emptyDistanceCache;
	contact->manifold = s2_emptyManifold;
	contact->friction = s2MixFriction(shapeA->friction, shapeB->friction);
	contact->restitution = s2MixRestitution(shapeA->restitution, shapeB->restitution);

	s2Body* bodyA = world->bodies + shapeA->bodyIndex;
	s2Body* bodyB = world->bodies + shapeB->bodyIndex;

	// Connect to body A
	{
		contact->edges[0].bodyIndex = shapeA->bodyIndex;
		contact->edges[0].prevKey = S2_NULL_INDEX;
		contact->edges[0].nextKey = bodyA->contactList;

		int32_t keyA = (contactIndex << 1) | 0;
		if (bodyA->contactList != S2_NULL_INDEX)
		{
			s2Contact* contactA = world->contacts + (bodyA->contactList >> 1);
			s2ContactEdge* edgeA = contactA->edges + (bodyA->contactList & 1);
			edgeA->prevKey = keyA;
		}
		bodyA->contactList = keyA;
		bodyA->contactCount += 1;
	}

	// Connect to body B
	{
		contact->edges[1].bodyIndex = shapeB->bodyIndex;
		contact->edges[1].prevKey = S2_NULL_INDEX;
		contact->edges[1].nextKey = bodyB->contactList;

		int32_t keyB = (contactIndex << 1) | 1;
		if (bodyB->contactList != S2_NULL_INDEX)
		{
			s2Contact* contactB = world->contacts + (bodyB->contactList >> 1);
			s2ContactEdge* edgeB = contactB->edges + (bodyB->contactList & 1);
			edgeB->prevKey = keyB;
		}
		bodyB->contactList = keyB;
		bodyB->contactCount += 1;
	}

	// Add to pair set for fast lookup
	uint64_t pairKey = S2_SHAPE_PAIR_KEY(contact->shapeIndexA, contact->shapeIndexB);
	s2AddKey(&world->broadPhase.pairSet, pairKey);
}

void s2DestroyContact(s2World* world, s2Contact* contact)
{
	// Remove pair from set
	uint64_t pairKey = S2_SHAPE_PAIR_KEY(contact->shapeIndexA, contact->shapeIndexB);
	s2RemoveKey(&world->broadPhase.pairSet, pairKey);

	s2ContactEdge* edgeA = contact->edges + 0;
	s2ContactEdge* edgeB = contact->edges + 1;

	s2Body* bodyA = world->bodies + edgeA->bodyIndex;
	s2Body* bodyB = world->bodies + edgeB->bodyIndex;

	// Remove from body A
	if (edgeA->prevKey != S2_NULL_INDEX)
	{
		s2Contact* prevContact = world->contacts + (edgeA->prevKey >> 1);
		s2ContactEdge* prevEdge = prevContact->edges + (edgeA->prevKey & 1);
		prevEdge->nextKey = edgeA->nextKey;
	}

	if (edgeA->nextKey != S2_NULL_INDEX)
	{
		s2Contact* nextContact = world->contacts + (edgeA->nextKey >> 1);
		s2ContactEdge* nextEdge = nextContact->edges + (edgeA->nextKey & 1);
		nextEdge->prevKey = edgeA->prevKey;
	}

	int32_t edgeKeyA = (contact->object.index << 1) | 0;
	if (bodyA->contactList == edgeKeyA)
	{
		bodyA->contactList = edgeA->nextKey;
	}

	bodyA->contactCount -= 1;

	// Remove from body B
	if (edgeB->prevKey != S2_NULL_INDEX)
	{
		s2Contact* prevContact = world->contacts + (edgeB->prevKey >> 1);
		s2ContactEdge* prevEdge = prevContact->edges + (edgeB->prevKey & 1);
		prevEdge->nextKey = edgeB->nextKey;
	}

	if (edgeB->nextKey != S2_NULL_INDEX)
	{
		s2Contact* nextContact = world->contacts + (edgeB->nextKey >> 1);
		s2ContactEdge* nextEdge = nextContact->edges + (edgeB->nextKey & 1);
		nextEdge->prevKey = edgeB->prevKey;
	}

	int32_t contactIndex = contact->object.index;

	int32_t edgeKeyB = (contactIndex << 1) | 1;
	if (bodyB->contactList == edgeKeyB)
	{
		bodyB->contactList = edgeB->nextKey;
	}

	bodyB->contactCount -= 1;

	s2FreeObject(&world->contactPool, &contact->object);
}

// Update the contact manifold and touching status.
// Note: do not assume the fixture AABBs are overlapping or are valid.
void s2UpdateContact(s2World* world, s2Contact* contact, s2Shape* shapeA, s2Body* bodyA, s2Shape* shapeB, s2Body* bodyB)
{
	s2Manifold oldManifold = contact->manifold;

	S2_ASSERT(shapeA->object.index == contact->shapeIndexA);
	S2_ASSERT(shapeB->object.index == contact->shapeIndexB);

	bool touching = false;
	contact->manifold.pointCount = 0;

	// bool wasTouching = (contact->flags & s2_contactTouchingFlag) == s2_contactTouchingFlag;

	s2ManifoldFcn* fcn = s_registers[shapeA->type][shapeB->type].fcn;

	float maxDistance = s2_speculativeDistance;
	s2Transform transformA = {bodyA->origin, bodyA->rot};
	s2Transform transformB = {bodyB->origin, bodyB->rot};
	contact->manifold = fcn(shapeA, transformA, shapeB, transformB, maxDistance, &contact->cache);

	touching = contact->manifold.pointCount > 0;

	contact->manifold.frictionPersisted = true;
	
	if (contact->manifold.pointCount != oldManifold.pointCount)
	{
		contact->manifold.frictionPersisted = false;
	}

	// TODO_ERIN testing
	contact->manifold.constraintIndex = oldManifold.constraintIndex;

	// Match old contact ids to new contact ids and copy the
	// stored impulses to warm start the solver.
	for (int32_t i = 0; i < contact->manifold.pointCount; ++i)
	{
		s2ManifoldPoint* mp2 = contact->manifold.points + i;
		mp2->normalImpulse = 0.0f;
		mp2->tangentImpulse = 0.0f;
		mp2->persisted = false;
		uint16_t id2 = mp2->id;

		for (int32_t j = 0; j < oldManifold.pointCount; ++j)
		{
			s2ManifoldPoint* mp1 = oldManifold.points + j;

			if (mp1->id == id2)
			{
				mp2->localNormalA = mp1->localNormalA;
				mp2->localNormalB = mp1->localNormalB;
				mp2->localAnchorA = mp1->localAnchorA;
				mp2->localAnchorB = mp1->localAnchorB;

				mp2->normalImpulse = mp1->normalImpulse;
				mp2->tangentImpulse = mp1->tangentImpulse;
				mp2->persisted = true;
				break;
			}
		}

		if (mp2->persisted == false)
		{
			contact->manifold.frictionPersisted = false;
		}
	}
}
