// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "allocate.h"
#include "body.h"
#include "contact.h"
#include "core.h"
#include "joint.h"
#include "solvers.h"
#include "stack_allocator.h"
#include "world.h"

#include <stdbool.h>

// This uses fixed anchors
static void s2SolveContacts_Jacobi_Soft(s2World* world, s2ContactConstraint* constraints, int constraintCount, float inv_h,
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

		s2Vec2 normal = constraint->normal;
		s2Vec2 tangent = s2RightPerp(normal);
		float friction = constraint->friction;

		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			float bias = 0.0f;
			float massScale = 1.0f;
			float impulseScale = 0.0f;
			if (cp->separation > 0.0f)
			{
				// Speculative
				bias = cp->separation * inv_h;
			}
			else if (useBias)
			{
				bias = S2_MAX(cp->biasCoefficient * cp->separation, -s2_maxBaumgarteVelocity);
				massScale = cp->massCoefficient;
				impulseScale = cp->impulseCoefficient;
			}

			// static anchors
			s2Vec2 rA = cp->rA0;
			s2Vec2 rB = cp->rB0;

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

			// static anchors
			s2Vec2 rA = cp->rA0;
			s2Vec2 rB = cp->rB0;

			// Relative velocity at contact
			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, rB));
			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, rA));
			s2Vec2 dv = s2Sub(vrB, vrA);

			// Compute tangent force
			float vt = s2Dot(dv, tangent);
			float lambda = cp->tangentMass * (-vt);

			// Clamp the accumulated force
			float maxFriction = friction * cp->normalImpulse;
			float newImpulse = S2_CLAMP(cp->tangentImpulse + lambda, -maxFriction, maxFriction);
			lambda = newImpulse - cp->tangentImpulse;
			cp->tangentImpulse = newImpulse;

			// Apply contact impulse
			s2Vec2 P = s2MulSV(lambda, tangent);

			vA = s2MulSub(vA, mA, P);
			wA -= iA * s2Cross(rA, P);

			vB = s2MulAdd(vB, mB, P);
			wB += iB * s2Cross(rB, P);
		}

		bodyA->dv = s2Add(bodyA->dv, s2Sub(vA, bodyA->linearVelocity));
		bodyA->dw += wA - bodyA->angularVelocity;

		bodyB->dv = s2Add(bodyB->dv, s2Sub(vB, bodyB->linearVelocity));
		bodyB->dw += wB - bodyB->angularVelocity;
	}
}

void s2Solve_Jacobi(s2World* world, s2StepContext* context)
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

	int velocityIterations = context->iterations;
	int positionIterations = context->extraIterations;
	float h = context->dt;
	float inv_h = context->inv_dt;

	float contactHertz = S2_MIN(s2_contactHertz, 0.333f * inv_h);
	float jointHertz = S2_MIN(s2_jointHertz, 0.5f * inv_h);

	// Loops: body 2, constraint 2 + iterations

	// Reset delta velocities for Jacobi
	s2Body* bodies = world->bodies;
	int bodyCapacity = world->bodyPool.capacity;
	for (int i = 0; i < bodyCapacity; ++i)
	{
		s2Body* body = bodies + i;
		if (s2IsFree(&body->object))
		{
			continue;
		}

		body->dv = s2Vec2_zero;
		body->dw = 0.0f;
	}

	// body loop
	s2IntegrateVelocities(world, h);

	// constraint loop
	s2PrepareContacts_Soft(world, constraints, constraintCount, context, h, contactHertz);

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

		s2PrepareJoint_Soft(joint, context, h, jointHertz, context->warmStart);

		if (context->warmStart)
		{
			s2WarmStartJoint(joint, context);
		}
	}

	// constraint loop * velocityIterations
	bool useBias = true;
	for (int iteration = 0; iteration < velocityIterations; ++iteration)
	{
		for (int i = 0; i < jointCapacity; ++i)
		{
			s2Joint* joint = joints + i;
			if (s2IsFree(&joint->object))
			{
				continue;
			}

			s2SolveJoint_Soft(joint, context, h, inv_h, useBias);
		}

		s2SolveContacts_Jacobi_Soft(world, constraints, constraintCount, inv_h, useBias);

		for (int i = 0; i < bodyCapacity; ++i)
		{
			s2Body* body = bodies + i;
			if (s2IsFree(&body->object))
			{
				continue;
			}

			body->linearVelocity = s2Add(body->linearVelocity, body->dv);
			body->angularVelocity += body->dw;
			body->dv = s2Vec2_zero;
			body->dw = 0.0f;
		}
	}

	// Update positions from velocity
	// body loop
	s2IntegratePositions(world, h);

	// Relax
	// constraint loop * positionIterations
	useBias = false;
	for (int iteration = 0; iteration < positionIterations; ++iteration)
	{
		for (int i = 0; i < jointCapacity; ++i)
		{
			s2Joint* joint = joints + i;
			if (s2IsFree(&joint->object))
			{
				continue;
			}

			s2SolveJoint_Soft(joint, context, h, inv_h, useBias);
		}

		s2SolveContacts_Jacobi_Soft(world, constraints, constraintCount, inv_h, useBias);

		for (int i = 0; i < bodyCapacity; ++i)
		{
			s2Body* body = bodies + i;
			if (s2IsFree(&body->object))
			{
				continue;
			}

			body->linearVelocity = s2Add(body->linearVelocity, body->dv);
			body->angularVelocity += body->dw;
			body->dv = s2Vec2_zero;
			body->dw = 0.0f;
		}
	}

	// body loop
	s2FinalizePositions(world);

	// constraint loop
	s2StoreContactImpulses(constraints, constraintCount);

	s2FreeStackItem(world->stackAllocator, constraints);
}
