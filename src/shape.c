// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "shape.h"

#include "body.h"
#include "broad_phase.h"
#include "world.h"

s2Box s2Shape_ComputeAABB(const s2Shape* shape, s2Transform xf)
{
	switch (shape->type)
	{
		case s2_capsuleShape:
			return s2ComputeCapsuleAABB(&shape->capsule, xf);
		case s2_circleShape:
			return s2ComputeCircleAABB(&shape->circle, xf);
		case s2_polygonShape:
			return s2ComputePolygonAABB(&shape->polygon, xf);
		case s2_segmentShape:
			return s2ComputeSegmentAABB(&shape->segment, xf);
		default: {
			S2_ASSERT(false);
			s2Box empty = {xf.p, xf.p};
			return empty;
		}
	}
}

s2MassData s2Shape_ComputeMass(const s2Shape* shape)
{
	switch (shape->type)
	{
		case s2_capsuleShape:
			return s2ComputeCapsuleMass(&shape->capsule, shape->density);
		case s2_circleShape:
			return s2ComputeCircleMass(&shape->circle, shape->density);
		case s2_polygonShape:
			return s2ComputePolygonMass(&shape->polygon, shape->density);
		default: {
			S2_ASSERT(false);
			s2MassData data = {0};
			return data;
		}
	}
}

void s2Shape_CreateProxy(s2Shape* shape, s2BroadPhase* bp, s2BodyType type, s2Transform xf)
{
	// Create proxies in the broad-phase.
	shape->aabb = s2Shape_ComputeAABB(shape, xf);

	// Smaller margin for static bodies. Cannot be zero due to TOI tolerance.
	float margin = type == s2_staticBody ? 4.0f * s2_linearSlop : s2_aabbMargin;
	shape->fatAABB.lowerBound.x = shape->aabb.lowerBound.x - margin;
	shape->fatAABB.lowerBound.y = shape->aabb.lowerBound.y - margin;
	shape->fatAABB.upperBound.x = shape->aabb.upperBound.x + margin;
	shape->fatAABB.upperBound.y = shape->aabb.upperBound.y + margin;
	shape->enlargedAABB = false;

	shape->proxyKey = s2BroadPhase_CreateProxy(bp, type, shape->fatAABB, 1, shape->object.index);
	S2_ASSERT(S2_PROXY_TYPE(shape->proxyKey) < s2_bodyTypeCount);
}

void s2Shape_DestroyProxy(s2Shape* shape, s2BroadPhase* bp)
{
	s2BroadPhase_DestroyProxy(bp, shape->proxyKey);
	shape->proxyKey = S2_NULL_INDEX;
}

s2DistanceProxy s2Shape_MakeDistanceProxy(const s2Shape* shape)
{
	switch (shape->type)
	{
		case s2_capsuleShape:
			return s2MakeProxy(&shape->capsule.point1, 2, shape->capsule.radius);
		case s2_circleShape:
			return s2MakeProxy(&shape->circle.point, 1, shape->circle.radius);
		case s2_polygonShape:
			return s2MakeProxy(shape->polygon.vertices, shape->polygon.count, shape->polygon.radius);
		case s2_segmentShape:
			return s2MakeProxy(&shape->segment.point1, 2, 0.0f);
		default: {
			S2_ASSERT(false);
			s2DistanceProxy empty = {0};
			return empty;
		}
	}
}

s2BodyId s2Shape_GetBody(s2ShapeId shapeId)
{
	s2World* world = s2GetWorldFromIndex(shapeId.world);
	S2_ASSERT(0 <= shapeId.index && shapeId.index < world->shapePool.capacity);
	s2Shape* shape = world->shapes + shapeId.index;
	S2_ASSERT(s2ObjectValid(&shape->object));

	S2_ASSERT(0 <= shape->bodyIndex && shape->bodyIndex < world->bodyPool.capacity);
	s2Body* body = world->bodies + shape->bodyIndex;
	S2_ASSERT(s2ObjectValid(&body->object));

	s2BodyId bodyId = {body->object.index, shapeId.world, body->object.revision};
	return bodyId;
}

bool s2Shape_TestPoint(s2ShapeId shapeId, s2Vec2 point)
{
	s2World* world = s2GetWorldFromIndex(shapeId.world);
	S2_ASSERT(0 <= shapeId.index && shapeId.index < world->shapePool.capacity);
	s2Shape* shape = world->shapes + shapeId.index;
	S2_ASSERT(s2ObjectValid(&shape->object));

	S2_ASSERT(0 <= shape->bodyIndex && shape->bodyIndex < world->bodyPool.capacity);
	s2Body* body = world->bodies + shape->bodyIndex;
	S2_ASSERT(s2ObjectValid(&body->object));

	s2Vec2 localPoint = s2InvTransformPoint(body->transform, point);

	switch (shape->type)
	{
		case s2_capsuleShape:
			return s2PointInCapsule(localPoint, &shape->capsule);

		case s2_circleShape:
			return s2PointInCircle(localPoint, &shape->circle);

		case s2_polygonShape:
			return s2PointInPolygon(localPoint, &shape->polygon);

		default:
			return false;
	}
}
