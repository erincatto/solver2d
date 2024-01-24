// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "body.h"

#include "array.h"
#include "block_allocator.h"
#include "contact.h"
#include "core.h"
#include "joint.h"
#include "shape.h"
#include "world.h"

#include "solver2d/aabb.h"
#include "solver2d/id.h"

s2BodyId s2CreateBody(s2WorldId worldId, const s2BodyDef* def)
{
	s2World* world = s2GetWorldFromId(worldId);

	s2Body* b = (s2Body*)s2AllocObject(&world->bodyPool);
	world->bodies = (s2Body*)world->bodyPool.memory;

	S2_ASSERT(0 <= def->type && def->type < s2_bodyTypeCount);
	S2_ASSERT(s2IsValidVec2(def->position));
	S2_ASSERT(s2IsValid(def->angle));
	S2_ASSERT(s2IsValidVec2(def->linearVelocity));
	S2_ASSERT(s2IsValid(def->angularVelocity));

	b->type = def->type;
	b->origin = def->position;
	b->position0 = def->position;
	b->position = def->position;
	b->rot0 = s2MakeRot(def->angle);
	b->rot = b->rot0;
	b->localCenter = s2Vec2_zero;
	b->linearVelocity = def->linearVelocity;
	b->angularVelocity = def->angularVelocity;
	b->linearVelocity0 = def->linearVelocity;
	b->angularVelocity0 = def->angularVelocity;
	b->force = s2Vec2_zero;
	b->torque = 0.0f;
	b->shapeList = S2_NULL_INDEX;
	b->jointList = S2_NULL_INDEX;
	b->jointCount = 0;
	b->contactList = S2_NULL_INDEX;
	b->contactCount = 0;
	b->mass = 0.0f;
	b->invMass = 0.0f;
	b->I = 0.0f;
	b->invI = 0.0f;
	b->linearDamping = def->linearDamping;
	b->angularDamping = def->angularDamping;
	b->gravityScale = def->gravityScale;

	b->userData = def->userData;
	b->world = worldId.index;
	b->enlargeAABB = false;

	s2BodyId id = {b->object.index, worldId.index, b->object.revision};
	return id;
}

s2Body* s2GetBody(s2World* world, s2BodyId id)
{
	S2_ASSERT(0 <= id.index && id.index < world->bodyPool.capacity);
	s2Body* body = world->bodies + id.index;
	S2_ASSERT(s2ObjectValid(&body->object));
	S2_ASSERT(id.revision == body->object.revision);
	return body;
}

void s2DestroyBody(s2BodyId bodyId)
{
	s2World* world = s2GetWorldFromIndex(bodyId.world);

	S2_ASSERT(0 <= bodyId.index && bodyId.index < world->bodyPool.capacity);

	s2Body* body = world->bodies + bodyId.index;

	// User must destroy joints before destroying bodies
	S2_ASSERT(body->jointList == S2_NULL_INDEX && body->jointCount == 0);

	// Destroy the attached contacts
	int32_t edgeKey = body->contactList;
	while (edgeKey != S2_NULL_INDEX)
	{
		int32_t contactIndex = edgeKey >> 1;
		int32_t edgeIndex = edgeKey & 1;

		int32_t twinKey = edgeKey ^ 1;
		int32_t twinIndex = twinKey & 1;

		s2Contact* contact = world->contacts + contactIndex;

		s2ContactEdge* twin = contact->edges + twinIndex;

		// Remove contact from other body's doubly linked list
		if (twin->prevKey != S2_NULL_INDEX)
		{
			s2Contact* prevContact = world->contacts + (twin->prevKey >> 1);
			s2ContactEdge* prevEdge = prevContact->edges + (twin->prevKey & 1);
			prevEdge->nextKey = twin->nextKey;
		}

		if (twin->nextKey != S2_NULL_INDEX)
		{
			s2Contact* nextContact = world->contacts + (twin->nextKey >> 1);
			s2ContactEdge* nextEdge = nextContact->edges + (twin->nextKey & 1);
			nextEdge->prevKey = twin->prevKey;
		}

		// Check other body's list head
		s2Body* other = world->bodies + twin->bodyIndex;
		if (other->contactList == twinKey)
		{
			other->contactList = twin->nextKey;
		}

		S2_ASSERT(other->contactCount > 0);
		other->contactCount -= 1;

		// Remove pair from set
		uint64_t pairKey = S2_SHAPE_PAIR_KEY(contact->shapeIndexA, contact->shapeIndexB);
		s2RemoveKey(&world->broadPhase.pairSet, pairKey);

		s2ContactEdge* edge = contact->edges + edgeIndex;
		edgeKey = edge->nextKey;

		// Free contact
		s2FreeObject(&world->contactPool, &contact->object);
	}

	// Delete the attached shapes. This destroys broad-phase proxies.
	int32_t shapeIndex = body->shapeList;
	while (shapeIndex != S2_NULL_INDEX)
	{
		s2Shape* shape = world->shapes + shapeIndex;
		shapeIndex = shape->nextShapeIndex;

		// The broad-phase proxies only exist if the body is enabled
		s2Shape_DestroyProxy(shape, &world->broadPhase);

		s2FreeObject(&world->shapePool, &shape->object);
	}

	// Free body
	s2FreeObject(&world->bodyPool, &body->object);
}

static void s2ComputeMass(s2World* w, s2Body* b)
{
	// Compute mass data from shapes. Each shape has its own density.
	b->mass = 0.0f;
	b->invMass = 0.0f;
	b->I = 0.0f;
	b->invI = 0.0f;
	b->localCenter = s2Vec2_zero;

	// Static and kinematic bodies have zero mass.
	if (b->type == s2_staticBody || b->type == s2_kinematicBody)
	{
		b->position = b->origin;
		return;
	}

	S2_ASSERT(b->type == s2_dynamicBody);

	// Accumulate mass over all shapes.
	s2Vec2 localCenter = s2Vec2_zero;
	int32_t shapeIndex = b->shapeList;
	while (shapeIndex != S2_NULL_INDEX)
	{
		const s2Shape* s = w->shapes + shapeIndex;
		shapeIndex = s->nextShapeIndex;

		if (s->density == 0.0f)
		{
			continue;
		}

		s2MassData massData = s2Shape_ComputeMass(s);

		b->mass += massData.mass;
		localCenter = s2MulAdd(localCenter, massData.mass, massData.center);
		b->I += massData.I;
	}

	// Compute center of mass.
	if (b->mass > 0.0f)
	{
		b->invMass = 1.0f / b->mass;
		localCenter = s2MulSV(b->invMass, localCenter);
	}

	if (b->I > 0.0f)
	{
		// Center the inertia about the center of mass.
		b->I -= b->mass * s2Dot(localCenter, localCenter);
		S2_ASSERT(b->I > 0.0f);
		b->invI = 1.0f / b->I;
	}
	else
	{
		b->I = 0.0f;
		b->invI = 0.0f;
	}

	// Move center of mass.
	s2Vec2 oldCenter = b->position;
	b->localCenter = localCenter;
	b->position = s2Add(s2RotateVector(b->rot, b->localCenter), b->origin);

	// Update center of mass velocity.
	s2Vec2 deltaLinear = s2CrossSV(b->angularVelocity, s2Sub(b->position, oldCenter));
	b->linearVelocity = s2Add(b->linearVelocity, deltaLinear);
}

static s2ShapeId s2CreateShape(s2BodyId bodyId, const s2ShapeDef* def, const void* geometry, s2ShapeType shapeType)
{
	s2World* w = s2GetWorldFromIndex(bodyId.world);

	S2_ASSERT(0 <= bodyId.index && bodyId.index < w->bodyPool.capacity);

	s2Body* body = w->bodies + bodyId.index;

	s2Shape* shape = (s2Shape*)s2AllocObject(&w->shapePool);
	w->shapes = (s2Shape*)w->shapePool.memory;

	S2_ASSERT(s2IsValid(def->density) && def->density >= 0.0f);
	S2_ASSERT(s2IsValid(def->friction) && def->friction >= 0.0f);
	S2_ASSERT(s2IsValid(def->restitution) && def->restitution >= 0.0f);

	switch (shapeType)
	{
		case s2_capsuleShape:
			shape->capsule = *(const s2Capsule*)geometry;
			break;

		case s2_circleShape:
			shape->circle = *(const s2Circle*)geometry;
			break;

		case s2_polygonShape:
			shape->polygon = *(const s2Polygon*)geometry;
			break;

		case s2_segmentShape:
			shape->segment = *(const s2Segment*)geometry;
			break;

		default:
			S2_ASSERT(false);
			break;
	}

	shape->bodyIndex = body->object.index;
	shape->type = shapeType;
	shape->density = def->density;
	shape->friction = def->friction;
	shape->restitution = def->restitution;
	shape->userData = def->userData;
	shape->filter = def->filter;
	shape->enlargedAABB = false;

	s2Shape_CreateProxy(shape, &w->broadPhase, body->type, (s2Transform){body->origin, body->rot});

	// Add to shape linked list
	shape->nextShapeIndex = body->shapeList;
	body->shapeList = shape->object.index;

	if (shape->density)
	{
		s2ComputeMass(w, body);
	}

	s2ShapeId id = {shape->object.index, bodyId.world, shape->object.revision};
	return id;
}

s2ShapeId s2CreateCircleShape(s2BodyId bodyId, const s2ShapeDef* def, const s2Circle* circle)
{
	return s2CreateShape(bodyId, def, circle, s2_circleShape);
}

s2ShapeId s2CreatePolygonShape(s2BodyId bodyId, const s2ShapeDef* def, const s2Polygon* polygon)
{
	return s2CreateShape(bodyId, def, polygon, s2_polygonShape);
}

s2ShapeId s2CreateSegmentShape(s2BodyId bodyId, const s2ShapeDef* def, const s2Segment* segment)
{
	float lengthSqr = s2DistanceSquared(segment->point1, segment->point2);
	if (lengthSqr <= s2_linearSlop * s2_linearSlop)
	{
		S2_ASSERT(false);
		return s2_nullShapeId;
	}

	return s2CreateShape(bodyId, def, segment, s2_segmentShape);
}

s2ShapeId s2CreateCapsuleShape(s2BodyId bodyId, const s2ShapeDef* def, const s2Capsule* capsule)
{
	float lengthSqr = s2DistanceSquared(capsule->point1, capsule->point2);
	if (lengthSqr <= s2_linearSlop * s2_linearSlop)
	{
		S2_ASSERT(false);
		return s2_nullShapeId;
	}

	return s2CreateShape(bodyId, def, capsule, s2_capsuleShape);
}

s2Vec2 s2Body_GetPosition(s2BodyId bodyId)
{
	s2World* world = s2GetWorldFromIndex(bodyId.world);
	S2_ASSERT(0 <= bodyId.index && bodyId.index < world->bodyPool.capacity);
	return world->bodies[bodyId.index].origin;
}

float s2Body_GetAngle(s2BodyId bodyId)
{
	s2World* world = s2GetWorldFromIndex(bodyId.world);
	S2_ASSERT(0 <= bodyId.index && bodyId.index < world->bodyPool.capacity);
	return s2Rot_GetAngle(world->bodies[bodyId.index].rot);
}

s2Vec2 s2Body_GetLocalPoint(s2BodyId bodyId, s2Vec2 globalPoint)
{
	s2World* world = s2GetWorldFromIndex(bodyId.world);
	s2Body* body = s2GetBody(world, bodyId);
	return s2InvTransformPoint((s2Transform){body->origin, body->rot}, globalPoint);
}

void s2Body_SetLinearVelocity(s2BodyId bodyId, s2Vec2 linearVelocity)
{
	s2World* world = s2GetWorldFromIndex(bodyId.world);
	S2_ASSERT(0 <= bodyId.index && bodyId.index < world->bodyPool.capacity);
	world->bodies[bodyId.index].linearVelocity = linearVelocity;
}

void s2Body_SetAngularVelocity(s2BodyId bodyId, float angularVelocity)
{
	s2World* world = s2GetWorldFromIndex(bodyId.world);
	S2_ASSERT(0 <= bodyId.index && bodyId.index < world->bodyPool.capacity);

	world->bodies[bodyId.index].angularVelocity = angularVelocity;
}

void s2Body_ApplyForceToCenter(s2BodyId bodyId, s2Vec2 force)
{
	s2World* world = s2GetWorldFromIndex(bodyId.world);
	s2Body* body = s2GetBody(world, bodyId);
	body->force = s2Add(body->force, force);
}

void s2Body_ApplyLinearImpulse(s2BodyId bodyId, s2Vec2 impulse, s2Vec2 point)
{
	s2World* world = s2GetWorldFromIndex(bodyId.world);
	s2Body* body = s2GetBody(world, bodyId);
	if (body->type != s2_dynamicBody)
	{
		return;
	}

	body->linearVelocity = s2MulAdd(body->linearVelocity, body->invMass, impulse);
	body->angularVelocity += body->invI * s2Cross(s2Sub(point, body->position), impulse);
}

s2BodyType s2Body_GetType(s2BodyId bodyId)
{
	s2World* world = s2GetWorldFromIndex(bodyId.world);
	S2_ASSERT(0 <= bodyId.index && bodyId.index < world->bodyPool.capacity);
	return world->bodies[bodyId.index].type;
}

float s2Body_GetMass(s2BodyId bodyId)
{
	s2World* world = s2GetWorldFromIndex(bodyId.world);
	S2_ASSERT(0 <= bodyId.index && bodyId.index < world->bodyPool.capacity);
	return world->bodies[bodyId.index].mass;
}

bool s2ShouldBodiesCollide(s2World* world, s2Body* bodyA, s2Body* bodyB)
{
	int32_t jointKey;
	int32_t otherBodyIndex;
	if (bodyA->jointCount < bodyB->jointCount)
	{
		jointKey = bodyA->jointList;
		otherBodyIndex = bodyB->object.index;
	}
	else
	{
		jointKey = bodyB->jointList;
		otherBodyIndex = bodyA->object.index;
	}

	while (jointKey != S2_NULL_INDEX)
	{
		int32_t jointIndex = jointKey >> 1;
		int32_t edgeIndex = jointKey & 1;
		int32_t otherEdgeIndex = edgeIndex ^ 1;

		s2Joint* joint = world->joints + jointIndex;
		if (joint->edges[otherEdgeIndex].bodyIndex == otherBodyIndex)
		{
			return false;
		}

		jointKey = joint->edges[edgeIndex].nextKey;
	}

	return true;
}
