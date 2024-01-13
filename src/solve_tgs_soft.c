// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "allocate.h"
#include "array.h"
#include "body.h"
#include "contact.h"
#include "core.h"
#include "joint.h"
#include "shape.h"
#include "solvers.h"
#include "stack_allocator.h"
#include "world.h"

#include "solver2d/aabb.h"

#include <stdbool.h>
#include <stdlib.h>

#define maxBaumgarteVelocity 3.0f

static void s2UpdatePositionDeltas(s2World* world, float h)
{
	s2Body* bodies = world->bodies;
	int bodyCapacity = world->bodyPool.capacity;

	for (int i = 0; i < bodyCapacity; ++i)
	{
		s2Body* body = bodies + i;
		if (s2ObjectValid(&body->object) == false)
		{
			continue;
		}

		if (body->type == s2_staticBody)
		{
			continue;
		}

		body->deltaAngle += h * body->angularVelocity;
		body->deltaPosition = s2MulAdd(body->deltaPosition, h, body->linearVelocity);
	}
}

static void s2FinalizePositions(s2World* world)
{
	s2Body* bodies = world->bodies;
	int bodyCapacity = world->bodyPool.capacity;

	for (int i = 0; i < bodyCapacity; ++i)
	{
		s2Body* body = bodies + i;
		if (s2ObjectValid(&body->object) == false)
		{
			continue;
		}

		if (body->type == s2_staticBody)
		{
			continue;
		}

		body->position = s2Add(body->position, body->deltaPosition);
		body->angle += body->deltaAngle;
	}
}

// todo same as PGS
static void s2PrepareSoftContacts(s2World* world, s2ContactConstraint* constraints, int constraintCount, float h)
{
	s2Body* bodies = world->bodies;

	for (int i = 0; i < constraintCount; ++i)
	{
		s2ContactConstraint* constraint = constraints + i;

		s2Contact* contact = constraint->contact;
		const s2Manifold* manifold = &contact->manifold;
		int pointCount = manifold->pointCount;
		S2_ASSERT(0 < pointCount && pointCount <= 2);
		int indexA = contact->edges[0].bodyIndex;
		int indexB = contact->edges[1].bodyIndex;

		constraint->indexA = indexA;
		constraint->indexB = indexB;
		constraint->normal = manifold->normal;
		constraint->friction = contact->friction;
		constraint->pointCount = pointCount;

		s2Body* bodyA = bodies + indexA;
		s2Body* bodyB = bodies + indexB;

		float mA = bodyA->invMass;
		float iA = bodyA->invI;
		float mB = bodyB->invMass;
		float iB = bodyB->invI;

		s2Vec2 cA = bodyA->position;
		s2Vec2 cB = bodyB->position;
		s2Rot qA = s2MakeRot(bodyA->angle);
		s2Rot qB = s2MakeRot(bodyB->angle);

		s2Vec2 vA = bodyA->linearVelocity;
		float wA = bodyA->angularVelocity;
		s2Vec2 vB = bodyB->linearVelocity;
		float wB = bodyB->angularVelocity;

		s2Vec2 normal = constraint->normal;
		s2Vec2 tangent = s2RightPerp(constraint->normal);

		for (int j = 0; j < pointCount; ++j)
		{
			const s2ManifoldPoint* mp = manifold->points + j;
			s2ContactConstraintPoint* cp = constraint->points + j;

			cp->normalImpulse = mp->normalImpulse;
			cp->tangentImpulse = mp->tangentImpulse;

			cp->rA = s2Sub(mp->point, cA);
			cp->rB = s2Sub(mp->point, cB);
			cp->localAnchorA = s2InvRotateVector(qA, cp->rA);
			cp->localAnchorB = s2InvRotateVector(qB, cp->rB);

			float rnA = s2Cross(cp->rA, normal);
			float rnB = s2Cross(cp->rB, normal);
			float kNormal = mA + mB + iA * rnA * rnA + iB * rnB * rnB;

			float rtA = s2Cross(cp->rA, tangent);
			float rtB = s2Cross(cp->rB, tangent);
			float kTangent = mA + mB + iA * rtA * rtA + iB * rtB * rtB;

			cp->tangentMass = kTangent > 0.0f ? 1.0f / kTangent : 0.0f;

			// Soft contact with speculation
			const float hertz = 30.0f;
			const float zeta = 1.0f;
			float omega = 2.0f * s2_pi * hertz;
			// float d = 2.0f * zeta * omega / kNormal;
			// float k = omega * omega / kNormal;

			// cp->gamma = 1.0f / (h * (d + h * k));
			// cp->gamma = 1.0f / (h * (2.0f * zeta * omega / kNormal + h * omega * omega / kNormal));
			// cp->gamma = kNormal / (h * omega * (2.0f * zeta + h * omega));

			cp->separation = mp->separation;

			// cp->bias = h * k * cp->gamma * mp->separation;
			// cp->bias = k / (d + h * k) * mp->separation;
			// cp->bias =
			//	(omega * omega / kNormal) / (2 * zeta * omega / kNormal + h * omega * omega / kNormal) * mp->separation;
			cp->biasCoefficient = omega / (2.0f * zeta + h * omega);
			// cp->gamma = 0.0f;
			// cp->bias = (0.2f / h) * mp->separation;

			cp->normalMass = kNormal > 0.0f ? 1.0f / kNormal : 0.0f;

			float c = h * omega * (2.0f * zeta + h * omega);
			cp->impulseCoefficient = 1.0f / (1.0f + c);
			cp->massCoefficient = c * cp->impulseCoefficient;

			// meff = 1.0f / kNormal * 1.0f / (1.0f + 1.0f / (h * omega * (2 * zeta + h * omega)))
			// float impulse = -cp->normalMass * (vn + bias + cp->gamma * cp->normalImpulse);
			// = -meff * mscale * (vn + bias) - imp_scale * impulse

			// Warm start
			s2Vec2 P = s2Add(s2MulSV(cp->normalImpulse, normal), s2MulSV(cp->tangentImpulse, tangent));
			wA -= iA * s2Cross(cp->rA, P);
			vA = s2MulAdd(vA, -mA, P);
			wB += iB * s2Cross(cp->rB, P);
			vB = s2MulAdd(vB, mB, P);
		}

		bodyA->linearVelocity = vA;
		bodyA->angularVelocity = wA;
		bodyB->linearVelocity = vB;
		bodyB->angularVelocity = wB;
	}
}

static void s2SolveContactVelocities_TGS_Soft(s2World* world, s2ContactConstraint* constraints, int constraintCount, float inv_dt,
											  bool useBias)
{
	s2Body* bodies = world->bodies;

	for (int i = 0; i < constraintCount; ++i)
	{
		s2ContactConstraint* constraint = constraints + i;

		s2Body* bodyA = bodies + constraint->indexA;
		s2Body* bodyB = bodies + constraint->indexB;

		float mA = bodyA->invMass;
		float iA = bodyA->invI;
		float mB = bodyB->invMass;
		float iB = bodyB->invI;
		int pointCount = constraint->pointCount;

		s2Vec2 vA = bodyA->linearVelocity;
		float wA = bodyA->angularVelocity;
		s2Vec2 vB = bodyB->linearVelocity;
		float wB = bodyB->angularVelocity;

		const s2Vec2 dpA = bodyA->deltaPosition;
		const float daA = bodyA->deltaAngle;
		const s2Vec2 dpB = bodyB->deltaPosition;
		const float daB = bodyB->deltaAngle;

		s2Vec2 normal = constraint->normal;
		s2Vec2 tangent = s2RightPerp(normal);
		float friction = constraint->friction;

		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			// Approximate change in anchor points
			s2Vec2 drA = s2CrossSV(daA, cp->rA);
			s2Vec2 drB = s2CrossSV(daB, cp->rB);

			// Compute change in separation
			s2Vec2 prA = s2Add(dpA, drA);
			s2Vec2 prB = s2Add(dpB, drB);
			float ds = s2Dot(s2Sub(prB, prA), normal);
			float s = cp->separation + ds;

			float bias = 0.0f;
			float massScale = 1.0f;
			float impulseScale = 0.0f;
			if (s > 0.0f)
			{
				// Speculative
				bias = s * inv_dt;
			}
			else if (useBias)
			{
				bias = S2_MAX(cp->biasCoefficient * s, -maxBaumgarteVelocity);
				massScale = cp->massCoefficient;
				impulseScale = cp->impulseCoefficient;
			}

			// Current anchors
			s2Vec2 rA = s2Add(cp->rA, drA);
			s2Vec2 rB = s2Add(cp->rB, drB);

			// Relative velocity at contact
			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, rB));
			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, rA));
			float vn = s2Dot(s2Sub(vrB, vrA), normal);

			// Compute normal impulse
			float impulse = -cp->normalMass * massScale * (vn + bias) - impulseScale * cp->normalImpulse;

			// Clamp the accumulated impulse
			float newImpulse = S2_MAX(cp->normalImpulse + impulse, 0.0f);
			impulse = newImpulse - cp->normalImpulse;
			cp->normalImpulse = newImpulse;

			// Apply contact impulse
			s2Vec2 P = s2MulSV(impulse, normal);
			vA = s2MulSub(vA, mA, P);
			wA -= iA * s2Cross(rA, P);

			vB = s2MulAdd(vB, mB, P);
			wB += iB * s2Cross(rB, P);
		}

		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			// Approximate anchor points
			s2Vec2 rA = s2Add(cp->rA, s2CrossSV(daA, cp->rA));
			s2Vec2 rB = s2Add(cp->rB, s2CrossSV(daB, cp->rB));

			// Relative velocity at contact
			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, rB));
			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, rA));
			float vt = s2Dot(s2Sub(vrB, vrA), tangent);

			// Compute tangent force
			float impulse = -cp->tangentMass * vt;

			// Clamp the accumulated force
			float maxFriction = friction * cp->normalImpulse;
			float newImpulse = S2_CLAMP(cp->tangentImpulse + impulse, -maxFriction, maxFriction);
			impulse = newImpulse - cp->tangentImpulse;
			cp->tangentImpulse = newImpulse;

			// Apply contact impulse
			s2Vec2 P = s2MulSV(impulse, tangent);

			vA = s2MulSub(vA, mA, P);
			wA -= iA * s2Cross(rA, P);

			vB = s2MulAdd(vB, mB, P);
			wB += iB * s2Cross(rB, P);
		}

		bodyA->linearVelocity = vA;
		bodyA->angularVelocity = wA;
		bodyB->linearVelocity = vB;
		bodyB->angularVelocity = wB;
	}
}

void s2Solve_TGS_Soft(s2World* world, s2StepContext* context)
{
	s2Contact* contacts = world->contacts;
	int contactCapacity = world->contactPool.capacity;

	s2Joint* joints = world->joints;
	int jointCapacity = world->jointPool.capacity;

	s2ContactConstraint* constraints =
		s2AllocateStackItem(world->stackAllocator, contactCapacity * sizeof(s2ContactConstraint), "constraint");

	int constraintCount = 0;
	for (int i = 0; i < contactCapacity; ++i)
	{
		s2Contact* contact = contacts + i;
		if (s2IsFree(&contact->object))
		{
			continue;
		}

		if (contact->manifold.pointCount == 0)
		{
			continue;
		}

		constraints[constraintCount].contact = contact;
		constraints[constraintCount].contact->manifold.constraintIndex = constraintCount;
		constraintCount += 1;
	}

	// Full step apply gravity
	s2IntegrateVelocities(world, context->dt);

	s2PrepareSoftContacts(world, constraints, constraintCount, context->dt);

	for (int i = 0; i < jointCapacity; ++i)
	{
		s2Joint* joint = joints + i;
		if (s2IsFree(&joint->object))
		{
			continue;
		}

		// todo soft
		s2PrepareJoint(joint, context);
	}

	int substepCount = context->velocityIterations;
	float h = context->dt / substepCount;
	float inv_h = 1.0f / h;

	for (int substep = 0; substep < substepCount; ++substep)
	{
		for (int i = 0; i < jointCapacity; ++i)
		{
			s2Joint* joint = joints + i;
			if (s2IsFree(&joint->object))
			{
				continue;
			}

			// todo soft
			s2SolveJointVelocity(joint, context);
		}

		// One constraint iteration
		bool useBias = true;
		s2SolveContactVelocities_TGS_Soft(world, constraints, constraintCount, inv_h, useBias);

		// Update position deltas
		s2UpdatePositionDeltas(world, h);
	}

	// Update positions using deltas
	s2FinalizePositions(world);

	int positionIterations = context->positionIterations;
	for (int iter = 0; iter < positionIterations; ++iter)
	{
		for (int i = 0; i < jointCapacity; ++i)
		{
			s2Joint* joint = joints + i;
			if (s2IsFree(&joint->object))
			{
				continue;
			}

			// todo relax soft
			s2SolveJointPosition(joint, context);
		}

		bool useBias = false;
		s2SolveContactVelocities_TGS_Soft(world, constraints, constraintCount, 0.0f, useBias);
	}

	s2StoreContactImpulses(constraints, constraintCount);

	s2FreeStackItem(world->stackAllocator, constraints);
}
