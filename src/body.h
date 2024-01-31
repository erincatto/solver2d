// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "pool.h"

#include "solver2d/distance.h"
#include "solver2d/id.h"
#include "solver2d/math.h"

typedef struct s2Polygon s2Polygon;
typedef struct s2World s2World;

// A rigid body
typedef struct s2Body
{
	s2Object object;

	enum s2BodyType type;

	// the body origin (not center of mass)
	s2Vec2 origin;

	// center of mass position in world
	s2Vec2 position;

	// delta position for the whole time step
	s2Vec2 deltaPosition;

	// delta position at the beginning of each sub-step
	s2Vec2 deltaPosition0;

	// rotation
	s2Rot rot0;
	s2Rot rot;

	// location of center of mass relative to the body origin
	s2Vec2 localCenter;

	s2Vec2 linearVelocity;
	float angularVelocity;

	s2Vec2 linearVelocity0;
	float angularVelocity0;

	s2Vec2 force;
	float torque;

	int32_t shapeList;

	// This is a key: [jointIndex:31, edgeIndex:1]
	int32_t jointList;
	int32_t jointCount;

	int32_t contactList;
	int32_t contactCount;

	float mass, invMass;

	// Rotational inertia about the center of mass.
	float I, invI;

	float linearDamping;
	float angularDamping;
	float gravityScale;

	void* userData;
	int16_t world;

	bool enlargeAABB;
} s2Body;

bool s2ShouldBodiesCollide(s2World* world, s2Body* bodyA, s2Body* bodyB);

s2ShapeId s2CreatePolygonShape(s2BodyId bodyId, const s2ShapeDef* def, const s2Polygon* polygon);
void s2Body_DestroyShape(s2ShapeId shapeId);

#define S2_TRANSFORM(body)                                                                                                       \
	(s2Transform)                                                                                                                \
	{                                                                                                                            \
		body->origin, body->rot                                                                                                  \
	}
