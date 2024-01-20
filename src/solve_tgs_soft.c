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

// this differes from PGS because it uses updated anchors
static void s2SolveContacts_TGS_Soft(s2World* world, s2ContactConstraint* constraints, int constraintCount, float inv_h,
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

		s2Vec2 cA = bodyA->position;
		s2Vec2 cB = bodyB->position;
		s2Rot qA = bodyA->rot;
		s2Rot qB = bodyB->rot;

		s2Vec2 normal = constraint->normal;
		s2Vec2 tangent = s2RightPerp(normal);
		float friction = constraint->friction;

		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			// Current anchor points
			s2Vec2 rA = s2RotateVector(qA, cp->localAnchorA);
			s2Vec2 rB = s2RotateVector(qB, cp->localAnchorB);

			// Compute change in separation
			s2Vec2 d = s2Sub(s2Add(cB, rB), s2Add(cA, rA));
			float s = s2Dot(d, normal) + cp->separation;

			float bias = 0.0f;
			float massScale = 1.0f;
			float impulseScale = 0.0f;
			if (s > 0.0f)
			{
				// Speculative
				bias = s * inv_h;
			}
			else if (useBias)
			{
				bias = S2_MAX(cp->biasCoefficient * s, -s2_maxBaumgarteVelocity);
				massScale = cp->massCoefficient;
				impulseScale = cp->impulseCoefficient;
			}

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

			// Current anchor points
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

// Warm starting and relaxing in the substep loop
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

	int substepCount = context->iterations;
	float h = context->dt / substepCount;
	float inv_h = 1.0f / h;

	float contactHertz = S2_MIN(s2_contactHertz, 0.25f * inv_h);
	float jointHertz = S2_MIN(s2_jointHertz, 0.125f * inv_h);

	// Prepare
	s2PrepareContacts_Soft(world, constraints, constraintCount, context, h, contactHertz);

	for (int i = 0; i < jointCapacity; ++i)
	{
		s2Joint* joint = joints + i;
		if (s2IsFree(&joint->object))
		{
			continue;
		}

		bool warmStart = true;
		s2PrepareJoint_Soft(joint, context, h, jointHertz, warmStart);
	}

	// Solve
	for (int substep = 0; substep < substepCount; ++substep)
	{
		// Integrate gravity and forces
		s2IntegrateVelocities(world, h);

		// Apply warm starting
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

		// Solve velocities using position bias
		bool useBias = true;
		for (int i = 0; i < jointCapacity; ++i)
		{
			s2Joint* joint = joints + i;
			if (s2IsFree(&joint->object))
			{
				continue;
			}

			s2SolveJoint_Soft(joint, context, h, inv_h, useBias);
		}

		s2SolveContacts_TGS_Soft(world, constraints, constraintCount, inv_h, useBias);

		// Integrate positions using biased velocities
		s2IntegratePositions(world, h);

		// Relax biased velocities and impulses.
		// Relaxing the impulses reduces warm starting overshoot.
		useBias = false;
		for (int i = 0; i < jointCapacity; ++i)
		{
			s2Joint* joint = joints + i;
			if (s2IsFree(&joint->object))
			{
				continue;
			}

			s2SolveJoint_Soft(joint, context, h, inv_h, useBias);
		}
		
		s2SolveContacts_TGS_Soft(world, constraints, constraintCount, inv_h, useBias);
	}

	// Store results
	s2StoreContactImpulses(constraints, constraintCount);

	// no separate relax stage!

	s2FreeStackItem(world->stackAllocator, constraints);
}
