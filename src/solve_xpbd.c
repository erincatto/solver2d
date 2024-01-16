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

#include <assert.h>
#include <stdbool.h>

static void s2PrepareContacts_XPBD(s2World* world, s2ContactConstraint* constraints, int constraintCount)
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

		s2Vec2 normal = constraint->normal;
		s2Vec2 tangent = s2RightPerp(normal);

		for (int j = 0; j < pointCount; ++j)
		{
			const s2ManifoldPoint* mp = manifold->points + j;
			s2ContactConstraintPoint* cp = constraint->points + j;

			cp->normalImpulse = 0.0f;
			cp->tangentImpulse = 0.0f;

			s2Vec2 rA = s2Sub(mp->point, cA);
			s2Vec2 rB = s2Sub(mp->point, cB);
			cp->localAnchorA = s2InvRotateVector(qA, rA);
			cp->localAnchorB = s2InvRotateVector(qB, rB);
			cp->separation = mp->separation;

			cp->baumgarte = 0.0f;
			cp->biasCoefficient = 0.0f;

			// todo perhaps re-use effective mass across substeps

			float rtA = s2Cross(rA, tangent);
			float rtB = s2Cross(rB, tangent);
			float kTangent = mA + mB + iA * rtA * rtA + iB * rtB * rtB;
			cp->tangentMass = kTangent > 0.0f ? 1.0f / kTangent : 0.0f;

			float rnA = s2Cross(rA, normal);
			float rnB = s2Cross(rB, normal);
			float kNormal = mA + mB + iA * rnA * rnA + iB * rnB * rnB;
			cp->normalMass = kNormal > 0.0f ? 1.0f / kNormal : 0.0f;
		}
	}
}

static void s2SolveContactPositions_XPBD(s2World* world, s2ContactConstraint* constraints, int constraintCount, float h)
{
	s2Body* bodies = world->bodies;
	float inv_h = h > 0.0f ? 1.0f / h : 0.0f;

	// compliance because contacts are too energetic otherwise
	float baseCompliance = 0.00001f * inv_h* inv_h;
	// but the rush sample has too much overlap ...
	//float baseCompliance = 0.0f;

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

		float compliance = (mA == 0.0f || mB == 0.0f) ? 0.25f * baseCompliance : baseCompliance;

		s2Vec2 cA = bodyA->position;
		float aA = bodyA->angle;
		s2Vec2 cB = bodyB->position;
		float aB = bodyB->angle;

		s2Vec2 normal = constraint->normal;
		s2Vec2 tangent = s2CrossVS(normal, 1.0f);

		// non-penetration constraints
		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			s2Rot qA = s2MakeRot(aA);
			s2Rot qB = s2MakeRot(aB);

			s2Vec2 rA = s2RotateVector(qA, cp->localAnchorA);
			s2Vec2 rB = s2RotateVector(qB, cp->localAnchorB);

			// normal separation
			s2Vec2 dp = s2Sub(s2Add(cB, rB), s2Add(cA, rA));
			float C = s2Dot(dp, normal) + cp->separation;
			if (C > 0)
			{
				cp->normalImpulse = 0.0f;
				continue;
			}

			// this clamping is not in the paper, but it is used in other solvers
			//float C_clamped = S2_MAX(-s2_maxBaumgarteVelocity * h, C);

			float rnA = s2Cross(rA, normal);
			float rnB = s2Cross(rB, normal);

			// w1 and w2 in paper
			float kA = mA + iA * rnA * rnA;
			float kB = mB + iB * rnB * rnB;

			//float lambda = -C_clamped / (kA + kB + compliance);
			float lambda = -C / (kA + kB + compliance);
			cp->normalImpulse = lambda;

			s2Vec2 P = s2MulSV(lambda, normal);

			cA = s2MulSub(cA, mA, P);
			aA -= iA * s2Cross(rA, P);

			cB = s2MulAdd(cB, mB, P);
			aB += iB * s2Cross(rB, P);
		}

		// static friction constraints
		float friction = constraint->friction;

		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			s2Rot qA = s2MakeRot(aA);
			s2Rot qB = s2MakeRot(aB);

			s2Vec2 rA = s2RotateVector(qA, cp->localAnchorA);
			s2Vec2 rB = s2RotateVector(qB, cp->localAnchorB);

			// tangent separation
			s2Vec2 dp = s2Sub(s2Add(cB, rB), s2Add(cA, rA));
			float C = s2Dot(dp, tangent);

			float rtA = s2Cross(rA, tangent);
			float rtB = s2Cross(rB, tangent);

			// w1 and w2 in paper
			float kA = mA + iA * rtA * rtA;
			float kB = mB + iB * rtB * rtB;

			float lambda = -C / (kA + kB);

			float maxLambda = friction * cp->normalImpulse;
			if (lambda < -maxLambda || maxLambda < lambda)
			{
				cp->tangentImpulse = 0.0f;
				continue;
			}

			cp->tangentImpulse = lambda;

			s2Vec2 P = s2MulSV(lambda, tangent);

			cA = s2MulSub(cA, mA, P);
			aA -= iA * s2Cross(rA, P);

			cB = s2MulAdd(cB, mB, P);
			aB += iB * s2Cross(rB, P);
		}

		bodyA->position = cA;
		bodyA->angle = aA;
		bodyB->position = cB;
		bodyB->angle = aB;
	}
}

static void s2SolveContactVelocities_XPBD(s2World* world, s2ContactConstraint* constraints, int constraintCount, float h)
{
	s2Body* bodies = world->bodies;
	float threshold = 2.0f * s2Length(world->gravity) * h;
	float inv_h = h > 0.0f ? 1.0f / h : 0.0f;

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

		s2Vec2 cA = bodyA->position;
		float aA = bodyA->angle;
		s2Vec2 cB = bodyB->position;
		float aB = bodyB->angle;

		s2Rot qA = s2MakeRot(aA);
		s2Rot qB = s2MakeRot(aB);

		s2Vec2 vA = bodyA->linearVelocity;
		float wA = bodyA->angularVelocity;
		s2Vec2 vB = bodyB->linearVelocity;
		float wB = bodyB->angularVelocity;

		s2Vec2 vA0 = bodyA->linearVelocity0;
		float wA0 = bodyA->angularVelocity0;
		s2Vec2 vB0 = bodyB->linearVelocity0;
		float wB0 = bodyB->angularVelocity0;

		s2Vec2 normal = constraint->normal;
		s2Vec2 tangent = s2CrossVS(normal, 1.0f);
		float friction = constraint->friction;

		// relax non-penetration
		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;
			if (cp->normalImpulse == 0.0f)
			{
				continue;
			}

			s2Vec2 rA = s2RotateVector(qA, cp->localAnchorA);
			s2Vec2 rB = s2RotateVector(qB, cp->localAnchorB);

			// Relative velocity at contact
			s2Vec2 vrB0 = s2Add(vB0, s2CrossSV(wB0, rB));
			s2Vec2 vrA0 = s2Add(vA0, s2CrossSV(wA0, rA));
			s2Vec2 dv0 = s2Sub(vrB0, vrA0);

			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, rB));
			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, rA));
			s2Vec2 dv = s2Sub(vrB, vrA);

			float rnA = s2Cross(rA, normal);
			float rnB = s2Cross(rB, normal);

			// w1 and w2 in paper
			float kA = mA + iA * rnA * rnA;
			float kB = mB + iB * rnB * rnB;

			float vn0 = s2Dot(dv0, normal);
			float vn = s2Dot(dv, normal);

			float Cdot = vn;
			float lambda = -Cdot / (kA + kB);

			s2Vec2 P = s2MulSV(lambda, normal);
			vA = s2MulSub(vA, mA, P);
			wA -= iA * s2Cross(rA, P);
			vB = s2MulAdd(vB, mB, P);
			wB += iB * s2Cross(rB, P);
		}

		// kinetic friction
		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			s2Vec2 rA = s2RotateVector(qA, cp->localAnchorA);
			s2Vec2 rB = s2RotateVector(qB, cp->localAnchorB);

			// Relative velocity at contact
			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, rB));
			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, rA));
			s2Vec2 dv = s2Sub(vrB, vrA);

			// Compute tangent force
			float vt = s2Dot(dv, tangent);
			if (vt == 0.0f)
			{
				continue;
			}

			// eq 31
			cp->tangentImpulse = friction * cp->normalImpulse;
			float huf = cp->tangentImpulse * inv_h;
			float abs_vt = S2_ABS(vt);
			float Cdot = (vt / abs_vt) * S2_MIN(huf, abs_vt);

			float rtA = s2Cross(rA, tangent);
			float rtB = s2Cross(rB, tangent);

			// w1 and w2 in paper
			float kA = mA + iA * rtA * rtA;
			float kB = mB + iB * rtB * rtB;

			float lambda = -Cdot / (kA + kB);

			s2Vec2 P = s2MulSV(lambda, tangent);
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

// Detailed Rigid Body Simulation with Extended Position Based Dynamics, 2020
// Matthias Müller, Miles Macklin, Nuttapong Chentanez, Stefan Jeschke, Tae-Yong Kim
void s2Solve_XPBD(s2World* world, s2StepContext* context)
{
	int substepCount = context->iterations;
	if (substepCount == 0)
	{
		return;
	}

	if (context->dt == 0.0f)
	{
		return;
	}

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

	s2PrepareContacts_XPBD(world, constraints, constraintCount);

	for (int i = 0; i < jointCapacity; ++i)
	{
		s2Joint* joint = joints + i;
		if (s2IsFree(&joint->object))
		{
			continue;
		}
		s2PrepareJoint_XPBD(joint, context);
	}

	float h = context->dt / substepCount;
	float inv_h = 1.0f / h;
	s2Body* bodies = world->bodies;
	int bodyCapacity = world->bodyPool.capacity;
	s2Vec2 gravity = world->gravity;

	for (int substep = 0; substep < substepCount; ++substep)
	{
		// Integrate velocities and positions
		for (int i = 0; i < bodyCapacity; ++i)
		{
			s2Body* body = bodies + i;
			if (s2IsFree(&body->object))
			{
				continue;
			}

			if (body->type == s2_staticBody)
			{
				continue;
			}

			float invMass = body->invMass;
			float invI = body->invI;

			s2Vec2 v = body->linearVelocity;
			float w = body->angularVelocity;

			// integrate velocities
			v = s2Add(v, s2MulSV(h * invMass, s2MulAdd(body->force, body->mass * body->gravityScale, gravity)));
			w = w + h * invI * body->torque;

			// damping
			v = s2MulSV(1.0f / (1.0f + h * body->linearDamping), v);
			w *= 1.0f / (1.0f + h * body->angularDamping);

			body->linearVelocity = v;
			body->angularVelocity = w;

			s2Vec2 c = body->position;
			float a = body->angle;

			// store previous position
			body->position0 = c;
			body->angle0 = a;

			// integrate positions
			// this is unique to XPBD, no other solvers update position immediately
			body->position = s2MulAdd(c, h, v);
			body->angle += h * w;
		}

		for (int i = 0; i < jointCapacity; ++i)
		{
			s2Joint* joint = joints + i;
			if (s2IsFree(&joint->object))
			{
				continue;
			}

			s2SolveJoint_XPBD(joint, context, inv_h);
		}

		s2SolveContactPositions_XPBD(world, constraints, constraintCount, h);

		// Project velocities
		for (int i = 0; i < bodyCapacity; ++i)
		{
			s2Body* body = bodies + i;
			if (s2IsFree(&body->object))
			{
				continue;
			}

			if (body->type != s2_dynamicBody)
			{
				continue;
			}

			body->linearVelocity0 = body->linearVelocity;
			body->angularVelocity0 = body->angularVelocity;

			body->linearVelocity = s2MulSV(inv_h, s2Sub(body->position, body->position0));
			body->angularVelocity = inv_h * (body->angle - body->angle0);
		}

		// Relax contact velocities
		s2SolveContactVelocities_XPBD(world, constraints, constraintCount, h);
	}

	// warm starting is not used, this is just for reporting
	for (int i = 0; i < constraintCount; ++i)
	{
		s2ContactConstraint* constraint = constraints + i;
		s2Contact* contact = constraint->contact;
		s2Manifold* manifold = &contact->manifold;

		for (int j = 0; j < constraint->pointCount; ++j)
		{
			manifold->points[j].normalImpulse = constraint->points[j].normalImpulse * inv_h;
			manifold->points[j].tangentImpulse = constraint->points[j].tangentImpulse * inv_h;
		}
	}

	s2FreeStackItem(world->stackAllocator, constraints);
}
