// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "probe.h"

#include "body.h"
#include "contact.h"
#include "joint.h"
#include "world.h"

#include "solver2d/math.h"

#include <math.h>
#include <stdlib.h>

typedef struct BodyPose
{
	s2Vec2 center;
	s2Rot rot;
	uint16_t revision;
	bool valid;
} BodyPose;

struct Probe
{
	BodyPose* poses;
	int capacity;
};

Probe* CreateProbe(void)
{
	Probe* probe = calloc(1, sizeof(Probe));
	return probe;
}

void DestroyProbe(Probe* probe)
{
	free(probe->poses);
	free(probe);
}

void Probe_Snapshot(Probe* probe, s2WorldId worldId)
{
	s2World* world = s2GetWorldFromId(worldId);
	int capacity = world->bodyPool.capacity;
	if (capacity > probe->capacity)
	{
		probe->poses = realloc(probe->poses, capacity * sizeof(BodyPose));
		probe->capacity = capacity;
	}

	for (int i = 0; i < probe->capacity; ++i)
	{
		BodyPose* pose = probe->poses + i;
		if (i >= capacity || s2IsFree(&world->bodies[i].object))
		{
			pose->valid = false;
			continue;
		}

		s2Body* body = world->bodies + i;
		pose->center = body->position;
		pose->rot = body->rot;
		pose->revision = body->object.revision;
		pose->valid = true;
	}
}

// Narrow phase runs before the solver, so the manifold describes the pose at the start of the step.
// Carry it forward to the solved pose the same way the solvers measure current separation.
static float CurrentSeparation(const Probe* probe, const s2Body* bodyA, const s2Body* bodyB, s2Vec2 normal,
							   const s2ManifoldPoint* mp)
{
	const BodyPose* poseA = probe->poses + bodyA->object.index;
	const BodyPose* poseB = probe->poses + bodyB->object.index;
	if (bodyA->object.index >= probe->capacity || bodyB->object.index >= probe->capacity || poseA->valid == false ||
		poseB->valid == false || poseA->revision != bodyA->object.revision || poseB->revision != bodyB->object.revision)
	{
		return mp->separation;
	}

	s2Vec2 localA = s2Sub(mp->localAnchorA, bodyA->localCenter);
	s2Vec2 localB = s2Sub(mp->localAnchorB, bodyB->localCenter);

	s2Vec2 rA0 = s2RotateVector(poseA->rot, localA);
	s2Vec2 rB0 = s2RotateVector(poseB->rot, localB);
	s2Vec2 rA = s2RotateVector(bodyA->rot, localA);
	s2Vec2 rB = s2RotateVector(bodyB->rot, localB);

	s2Vec2 dcA = s2Sub(bodyA->position, poseA->center);
	s2Vec2 dcB = s2Sub(bodyB->position, poseB->center);

	s2Vec2 d = s2Add(s2Sub(dcB, dcA), s2Sub(s2Sub(rB, rA), s2Sub(rB0, rA0)));
	return mp->separation + s2Dot(d, normal);
}

ProbeStep Probe_Measure(const Probe* probe, s2WorldId worldId)
{
	s2World* world = s2GetWorldFromId(worldId);
	ProbeStep result = {0};

	s2Body* bodies = world->bodies;
	int bodyCapacity = world->bodyPool.capacity;
	for (int i = 0; i < bodyCapacity; ++i)
	{
		s2Body* body = bodies + i;
		if (s2IsFree(&body->object) || body->type != s2_dynamicBody)
		{
			continue;
		}

		result.bodyCount += 1;

		float v2 = s2LengthSquared(body->linearVelocity);
		float w = body->angularVelocity;
		result.kineticEnergy += 0.5f * (body->mass * v2 + body->I * w * w);
		result.potentialEnergy -= body->mass * body->gravityScale * s2Dot(world->gravity, body->position);

		float speed = sqrtf(v2);
		result.maxSpeed = speed > result.maxSpeed ? speed : result.maxSpeed;
		result.maxAngularSpeed = fabsf(w) > result.maxAngularSpeed ? fabsf(w) : result.maxAngularSpeed;
	}

	float penetrationSum = 0.0f;
	int contactCapacity = world->contactPool.capacity;
	for (int i = 0; i < contactCapacity; ++i)
	{
		s2Contact* contact = world->contacts + i;
		if (s2IsFree(&contact->object) || contact->manifold.pointCount == 0)
		{
			continue;
		}

		result.contactCount += 1;

		const s2Body* bodyA = bodies + contact->edges[0].bodyIndex;
		const s2Body* bodyB = bodies + contact->edges[1].bodyIndex;
		const s2Manifold* manifold = &contact->manifold;

		for (int j = 0; j < manifold->pointCount; ++j)
		{
			const s2ManifoldPoint* mp = manifold->points + j;
			float separation = CurrentSeparation(probe, bodyA, bodyB, manifold->normal, mp);
			float penetration = separation < 0.0f ? -separation : 0.0f;
			result.maxPenetration = penetration > result.maxPenetration ? penetration : result.maxPenetration;
			penetrationSum += penetration;
			result.normalImpulse += mp->normalImpulse;
			result.pointCount += 1;
		}
	}

	if (result.pointCount > 0)
	{
		result.meanPenetration = penetrationSum / result.pointCount;
	}

	float jointErrorSum = 0.0f;
	int jointCapacity = world->jointPool.capacity;
	for (int i = 0; i < jointCapacity; ++i)
	{
		s2Joint* joint = world->joints + i;
		if (s2IsFree(&joint->object) || joint->type != s2_revoluteJoint)
		{
			continue;
		}

		const s2Body* bodyA = bodies + joint->edges[0].bodyIndex;
		const s2Body* bodyB = bodies + joint->edges[1].bodyIndex;

		s2Vec2 pA = s2TransformPoint(S2_TRANSFORM(bodyA), joint->localOriginAnchorA);
		s2Vec2 pB = s2TransformPoint(S2_TRANSFORM(bodyB), joint->localOriginAnchorB);
		float error = s2Distance(pA, pB);

		result.jointCount += 1;
		jointErrorSum += error;
		result.maxJointError = error > result.maxJointError ? error : result.maxJointError;

		const s2RevoluteJoint* revolute = &joint->revoluteJoint;
		if (revolute->enableLimit)
		{
			float angle = s2RelativeAngle(bodyB->rot, bodyA->rot) - revolute->referenceAngle;
			float limitError = 0.0f;
			if (angle < revolute->lowerAngle)
			{
				limitError = revolute->lowerAngle - angle;
			}
			else if (angle > revolute->upperAngle)
			{
				limitError = angle - revolute->upperAngle;
			}

			result.maxLimitError = limitError > result.maxLimitError ? limitError : result.maxLimitError;
		}
	}

	if (result.jointCount > 0)
	{
		result.meanJointError = jointErrorSum / result.jointCount;
	}

	return result;
}

int Probe_GetBodies(s2WorldId worldId, ProbeBody* bodies, int capacity)
{
	s2World* world = s2GetWorldFromId(worldId);
	int count = 0;
	int bodyCapacity = world->bodyPool.capacity;
	for (int i = 0; i < bodyCapacity; ++i)
	{
		s2Body* body = world->bodies + i;
		if (s2IsFree(&body->object) || body->type != s2_dynamicBody)
		{
			continue;
		}

		if (count < capacity)
		{
			ProbeBody* pb = bodies + count;
			pb->index = i;
			pb->position = body->position;
			pb->angle = s2Rot_GetAngle(body->rot);
			pb->linearVelocity = body->linearVelocity;
			pb->angularVelocity = body->angularVelocity;
		}

		count += 1;
	}

	return count;
}
