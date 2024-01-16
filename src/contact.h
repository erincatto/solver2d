// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "pool.h"

#include "solver2d/distance.h"
#include "solver2d/manifold.h"
#include "solver2d/types.h"

typedef struct s2Body s2Body;
typedef struct s2Shape s2Shape;
typedef struct s2World s2World;

// A contact edge is used to connect bodies and contacts together
// in a contact graph where each body is a node and each contact
// is an edge. A contact edge belongs to a doubly linked list
// maintained in each attached body. Each contact has two contact
// edges, one for each attached body.
typedef struct s2ContactEdge
{
	int32_t bodyIndex;
	int32_t prevKey;
	int32_t nextKey;
} s2ContactEdge;

// Flags stored in s2Contact::flags
enum s2ContactFlags
{
	// This contact no longer has overlapping AABBs
	s2_contactDisjoint = 0x00000020,

	// This contact started touching
	s2_contactStartedTouching = 0x00000040,

	// This contact stopped touching
	s2_contactStoppedTouching = 0x00000080,
};

/// The class manages contact between two shapes. A contact exists for each overlapping
/// AABB in the broad-phase (except if filtered). Therefore a contact object may exist
/// that has no contact points.
typedef struct s2Contact
{
	s2Object object;

	uint32_t flags;

	s2ContactEdge edges[2];

	int32_t shapeIndexA;
	int32_t shapeIndexB;

	s2DistanceCache cache;
	s2Manifold manifold;

	// Mixed friction and restitution
	float friction;
	float restitution;
} s2Contact;

void s2InitializeContactRegisters(void);

void s2CreateContact(s2World* world, s2Shape* shapeA, s2Shape* shapeB);
void s2DestroyContact(s2World* world, s2Contact* contact);

void s2UpdateContact(s2World* world, s2Contact* contact, s2Shape* shapeA, s2Body* bodyA, s2Shape* shapeB, s2Body* bodyB);

static inline bool s2ShouldShapesCollide(s2Filter filterA, s2Filter filterB)
{
	if (filterA.groupIndex == filterB.groupIndex && filterA.groupIndex != 0)
	{
		return filterA.groupIndex > 0;
	}

	bool collide = (filterA.maskBits & filterB.categoryBits) != 0 && (filterA.categoryBits & filterB.maskBits) != 0;
	return collide;
}
