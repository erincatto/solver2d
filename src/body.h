// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver2d/distance.h"
#include "solver2d/id.h"
#include "solver2d/math.h"

#include "pool.h"

typedef struct s2Polygon s2Polygon;
typedef struct s2World s2World;

// A rigid body
typedef struct s2Body
{
	s2Object object;

	enum s2BodyType type;

	// the body origin transform (not center of mass)
	s2Transform transform;
	
	// center of mass position in world
	s2Vec2 position0;
	s2Vec2 position;

	// rotation in radians
	float angle0;
	float angle;

	// location of center of mass relative to the body origin
	s2Vec2 localCenter;

	s2Vec2 linearVelocity;
	float angularVelocity;

	s2Vec2 deltaPosition;
	float deltaAngle;

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

	float sleepTime;

	void* userData;
	int16_t world;

	bool isMarked;
	bool enlargeAABB;
} s2Body;

bool s2ShouldBodiesCollide(s2World* world, s2Body* bodyA, s2Body* bodyB);

s2ShapeId s2Body_CreatePolygon(s2BodyId bodyId, const s2ShapeDef* def, const s2Polygon* polygon);
void s2Body_DestroyShape(s2ShapeId shapeId);
