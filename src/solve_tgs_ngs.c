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

// Should switch to quaternions or incremental rotation updates
// cos(a + da) = cos(a)*cos(da) + sin(a)*sin(da) ~= cos(a) + da * sin(a)
// sin(a + da) = sin(a)*cos(da) + cos(a)*sin(da) ~= sin(a) + da * cos(a)

static void s2PrepareContacts(s2World* world, s2ContactConstraint* constraints, int constraintCount, bool warmStart)
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

			if (warmStart)
			{
				cp->normalImpulse = mp->normalImpulse;
				cp->tangentImpulse = mp->tangentImpulse;
			}
			else
			{
				cp->normalImpulse = 0.0f;
				cp->tangentImpulse = 0.0f;
			}

			s2Vec2 rA = s2Sub(mp->point, cA);
			s2Vec2 rB = s2Sub(mp->point, cB);
			cp->localAnchorA = s2InvRotateVector(qA, rA);
			cp->localAnchorB = s2InvRotateVector(qB, rB);
			cp->separation = mp->separation;

			float rnA = s2Cross(rA, normal);
			float rnB = s2Cross(rB, normal);
			float kNormal = mA + mB + iA * rnA * rnA + iB * rnB * rnB;
			cp->normalMass = kNormal > 0.0f ? 1.0f / kNormal : 0.0f;

			float rtA = s2Cross(rA, tangent);
			float rtB = s2Cross(rB, tangent);
			float kTangent = mA + mB + iA * rtA * rtA + iB * rtB * rtB;
			cp->tangentMass = kTangent > 0.0f ? 1.0f / kTangent : 0.0f;
		}

		bodyA->linearVelocity = vA;
		bodyA->angularVelocity = wA;
		bodyB->linearVelocity = vB;
		bodyB->angularVelocity = wB;
	}
}

static void s2SolveContacts_TGS(s2World* world, s2ContactConstraint* constraints, int constraintCount, float inv_h)
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

		s2Vec2 cA = bodyA->position;
		s2Vec2 cB = bodyB->position;
		s2Rot qA = s2MakeRot(bodyA->angle);
		s2Rot qB = s2MakeRot(bodyB->angle);

		s2Vec2 normal = constraint->normal;
		s2Vec2 tangent = s2RightPerp(normal);
		float friction = constraint->friction;

		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			// Current anchors
			s2Vec2 rA = s2RotateVector(qA, cp->localAnchorA);
			s2Vec2 rB = s2RotateVector(qB, cp->localAnchorB);

			// Compute separation
			s2Vec2 d = s2Sub(s2Add(cB, rB), s2Add(cA, rA));
			float s = s2Dot(d, normal) + cp->separation;
			float bias = 0.0f;
			if (s > 0.0f)
			{
				// Speculative
				bias = s * inv_h;
			}

			// Relative velocity at contact
			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, rB));
			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, rA));
			float vn = s2Dot(s2Sub(vrB, vrA), normal);

			// Compute normal impulse
			float impulse = -cp->normalMass * (vn + bias);

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

			// Current anchors
			s2Vec2 rA = s2RotateVector(qA, cp->localAnchorA);
			s2Vec2 rB = s2RotateVector(qB, cp->localAnchorB);

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

void s2Solve_TGS_NGS(s2World* world, s2StepContext* context)
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

	s2PrepareContacts(world, constraints, constraintCount, context->warmStart);

	for (int i = 0; i < jointCapacity; ++i)
	{
		s2Joint* joint = joints + i;
		if (s2IsFree(&joint->object))
		{
			continue;
		}

		s2PrepareJoint(joint, context);
	}

	int substepCount = context->iterations;
	float h = context->dt / substepCount;
	float inv_h = 1.0f / h;

	for (int substep = 0; substep < substepCount; ++substep)
	{
		s2IntegrateVelocities(world, h);

		if (context->warmStart)
		{
			for (int i = 0; i < jointCapacity; ++i)
			{
				s2Joint* joint = joints + i;
				if (s2IsFree(&joint->object))
				{
					continue;
				}

				s2WarmStartJoint(joint, context);
			}

			s2WarmStartContacts(world, constraints, constraintCount);
		}

		for (int i = 0; i < jointCapacity; ++i)
		{
			s2Joint* joint = joints + i;
			if (s2IsFree(&joint->object))
			{
				continue;
			}

			s2SolveJoint(joint, context, h);
		}

		s2SolveContacts_TGS(world, constraints, constraintCount, inv_h);
		
		s2IntegratePositions(world, h);

		for (int i = 0; i < jointCapacity; ++i)
		{
			s2Joint* joint = joints + i;
			if (s2IsFree(&joint->object))
			{
				continue;
			}

			s2SolveJointPosition(joint, context);
		}

		s2SolveContact_NGS(world, constraints, constraintCount);
	}

	s2StoreContactImpulses(constraints, constraintCount);

	s2FreeStackItem(world->stackAllocator, constraints);
}
