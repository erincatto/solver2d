// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "body.h"
#include "contact.h"
#include "joint.h"
#include "solvers.h"
#include "stack_allocator.h"
#include "world.h"

#include "solver2d/aabb.h"

#include <stdbool.h>
#include <stdlib.h>

static void s2WarmStartContacts_Fixed(s2World* world, s2ContactConstraint* constraints, int constraintCount)
{
	s2Body* bodies = world->bodies;

	for (int i = 0; i < constraintCount; ++i)
	{
		s2ContactConstraint* constraint = constraints + i;

		int pointCount = constraint->pointCount;
		S2_ASSERT(0 < pointCount && pointCount <= 2);

		s2Body* bodyA = bodies + constraint->indexA;
		s2Body* bodyB = bodies + constraint->indexB;

		float mA = bodyA->invMass;
		float iA = bodyA->invI;
		float mB = bodyB->invMass;
		float iB = bodyB->invI;

		s2Vec2 vA = bodyA->linearVelocity;
		float wA = bodyA->angularVelocity;
		s2Vec2 vB = bodyB->linearVelocity;
		float wB = bodyB->angularVelocity;

		s2Vec2 normal = constraint->normal;
		s2Vec2 tangent = s2RightPerp(normal);

		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			// fixed anchors
			s2Vec2 rA = cp->rA0;
			s2Vec2 rB = cp->rB0;

			s2Vec2 P = s2Add(s2MulSV(cp->normalImpulse, normal), s2MulSV(cp->tangentImpulse, tangent));
			wA -= iA * s2Cross(rA, P);
			vA = s2MulAdd(vA, -mA, P);
			wB += iB * s2Cross(rB, P);
			vB = s2MulAdd(vB, mB, P);
		}

		bodyA->linearVelocity = vA;
		bodyA->angularVelocity = wA;
		bodyB->linearVelocity = vB;
		bodyB->angularVelocity = wB;
	}
}

// Uses fixed contact offsets which may be better for fast moving wheels
static void s2SolveContacts_TGS_Fixed(s2World* world, s2ContactConstraint* constraints, int constraintCount, float inv_h,
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

		s2Vec2 dcA = bodyA->deltaPosition;
		s2Rot qA = bodyA->rot;
		s2Vec2 dcB = bodyB->deltaPosition;
		s2Rot qB = bodyB->rot;

		s2Vec2 normal = constraint->normal;
		s2Vec2 tangent = s2RightPerp(normal);
		float friction = constraint->friction;

		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			// compute current separation (current anchors)
			s2Vec2 ds = s2Add(s2Sub(dcB, dcA), s2Sub(s2RotateVector(qB, cp->localAnchorB), s2RotateVector(qA, cp->localAnchorA)));
			float s = s2Dot(ds, normal) + cp->adjustedSeparation;

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
			// fixed anchors
			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, cp->rB0));
			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, cp->rA0));
			float vn = s2Dot(s2Sub(vrB, vrA), normal);
			
			// Compute normal impulse
			float impulse = -cp->normalMass * massScale * (vn + bias) - impulseScale * cp->normalImpulse;

			// Clamp the accumulated impulse
			float newImpulse = S2_MAX(cp->normalImpulse + impulse, 0.0f);
			impulse = newImpulse - cp->normalImpulse;
			cp->normalImpulse = newImpulse;

			// Apply contact impulse
			// fixed anchors
			s2Vec2 P = s2MulSV(impulse, normal);
			vA = s2MulSub(vA, mA, P);
			wA -= iA * s2Cross(cp->rA0, P);
			vB = s2MulAdd(vB, mB, P);
			wB += iB * s2Cross(cp->rB0, P);
		}

		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			// Relative velocity at contact
			// fixed anchors
			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, cp->rA0));
			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, cp->rB0));
			float vt = s2Dot(s2Sub(vrB, vrA), tangent);

			// Compute tangent force
			float impulse = -cp->tangentMass * vt;

			// Clamp the accumulated force
			float maxFriction = friction * cp->normalImpulse;
			float newImpulse = S2_CLAMP(cp->tangentImpulse + impulse, -maxFriction, maxFriction);
			impulse = newImpulse - cp->tangentImpulse;
			cp->tangentImpulse = newImpulse;

			// Apply contact impulse
			// fixed anchors
			s2Vec2 P = s2MulSV(impulse, tangent);
			vA = s2MulSub(vA, mA, P);
			wA -= iA * s2Cross(cp->rA0, P);
			vB = s2MulAdd(vB, mB, P);
			wB += iB * s2Cross(cp->rB0, P);
		}

		bodyA->linearVelocity = vA;
		bodyA->angularVelocity = wA;
		bodyB->linearVelocity = vB;
		bodyB->angularVelocity = wB;
	}
}

// This solver is almost identical to TGS_Soft except that contact anchors are not updated with sub-steps.
// This can improve rolling where the anchor should stay roughly fixed in direction.
// Like TGS_Soft this uses warm starting and relaxing in the substep loop
void s2Solve_SoftStep(s2World* world, s2StepContext* context)
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
	float h = context->h;
	float inv_h = context->inv_h;

	float contactHertz = S2_MIN(s2_contactHertz, 0.25f * inv_h);
	float jointHertz = S2_MIN(s2_jointHertz, 0.125f * inv_h);

	// Loops
	// body: 1 + 2 * substepCount
	// constraint: 2 + 2 * substepCount

	// Prepare
	// constraint loop
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
	// body 2 * substepCount
	// constraint 2 * substepCount (merge warm starting)
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

			s2WarmStartContacts_Fixed(world, constraints, constraintCount);
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

		s2SolveContacts_TGS_Fixed(world, constraints, constraintCount, inv_h, useBias);

		// Integrate positions using biased velocities
		s2IntegratePositions(world, h);

		// Relax biased velocities and impulses.
		// Relaxing the impulses reduces warm starting overshoot.
		if (context->extraIterations > 0)
		{
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
		
			s2SolveContacts_TGS_Fixed(world, constraints, constraintCount, inv_h, useBias);
		}
	}

	// Finalize body position
	// body loop
	s2FinalizePositions(world);

	// Store results
	// constraint loop
	s2StoreContactImpulses(constraints, constraintCount);

	s2FreeStackItem(world->stackAllocator, constraints);
}
