// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "allocate.h"
#include "array.h"
#include "body.h"
#include "contact.h"
#include "core.h"
#include "joint.h"
#include "solvers.h"
#include "stack_allocator.h"
#include "world.h"

#include "solver2d/aabb.h"

#include <stdbool.h>
#include <stdlib.h>

static void s2SolveContacts_PGS(s2World* world, s2ContactConstraint* constraints, int constraintCount, float inv_dt)
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

		s2Vec2 normal = constraint->normal;
		s2Vec2 tangent = s2CrossVS(normal, 1.0f);
		float friction = constraint->friction;

		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			// static anchors
			s2Vec2 rA = cp->rA0;
			s2Vec2 rB = cp->rB0;

			// Relative velocity at contact
			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, rB));
			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, rA));
			float vn = s2Dot(s2Sub(vrB, vrA), normal);

			// Compute normal impulse
			float impulse = -cp->normalMass * (vn + cp->biasCoefficient * cp->separation * inv_dt);

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

			// static anchors
			s2Vec2 rA = cp->rA0;
			s2Vec2 rB = cp->rB0;

			// Relative velocity at contact
			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, rB));
			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, rA));
			float vt = s2Dot(s2Sub(vrB, vrA), tangent);

			// Compute tangent force
			float lambda = cp->tangentMass * (-vt);

			// Clamp the accumulated force
			float maxFriction = friction * cp->normalImpulse;
			float newImpulse = S2_CLAMP(cp->tangentImpulse + lambda, -maxFriction, maxFriction);
			lambda = newImpulse - cp->tangentImpulse;
			cp->tangentImpulse = newImpulse;

			// Apply contact impulse
			s2Vec2 P = s2MulSV(lambda, tangent);

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

void s2Solve_PGS_NGS(s2World* world, s2StepContext* context)
{
	s2Contact* contacts = world->contacts;
	int contactCapacity = world->contactPool.capacity;

	s2Joint* joints = world->joints;
	int jointCapacity = world->jointPool.capacity;

	s2ContactConstraint* constraints = s2AllocateStackItem(world->stackAllocator, contactCapacity * sizeof(s2ContactConstraint), "constraint");
	
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

	int velocityIterations = context->iterations;
	int positionIterations = context->extraIterations;
	float h = context->dt;
	float inv_h = context->inv_dt;

	// Loops: body 3, constraint 2 + vel iter + pos iter

	// body loop
	s2IntegrateVelocities(world, h);

	// constraint loop
	s2PrepareContacts_PGS(world, constraints, constraintCount, context->warmStart);

	if (context->warmStart)
	{
		s2WarmStartContacts(world, constraints, constraintCount);
	}

	for (int i = 0; i < jointCapacity; ++i)
	{
		s2Joint* joint = joints + i;
		if (s2IsFree(&joint->object))
		{
			continue;
		}
		s2PrepareJoint(joint, context, context->warmStart);

		if (context->warmStart)
		{
			s2WarmStartJoint(joint, context);
		}
	}

	// constraint loop * velocityIterations
	for (int iter = 0; iter < velocityIterations; ++iter)
	{
		for (int i = 0; i < jointCapacity; ++i)
		{
			s2Joint* joint = joints + i;
			if (s2IsFree(&joint->object))
			{
				continue;
			}

			s2SolveJoint(joint, context, h);
		}

		s2SolveContacts_PGS(world, constraints, constraintCount, inv_h);
	}

	// body loop
	s2IntegratePositions(world, h);

	// constraint loop
	s2StoreContactImpulses(constraints, constraintCount);

	// constraint loop * positionIterations
	for (int iter = 0; iter < positionIterations; ++iter)
	{
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

	// body loop
	s2FinalizePositions(world);

	s2FreeStackItem(world->stackAllocator, constraints);
}
