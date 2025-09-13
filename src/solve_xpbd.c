// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "body.h"
#include "contact.h"
#include "core.h"
#include "joint.h"
#include "solvers.h"
#include "stack_allocator.h"
#include "world.h"

#include "solver2d/aabb.h"

#include <assert.h>
#include <float.h>
#include <stdbool.h>

static void s2PrepareContacts_XPBD(s2World* world, s2ContactConstraint* constraints, int constraintCount, float h)
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

		s2Rot qA = bodyA->rot;
		s2Rot qB = bodyB->rot;

		s2Vec2 normal = constraint->normal;
		s2Vec2 tangent = s2RightPerp(normal);

		for (int j = 0; j < pointCount; ++j)
		{
			const s2ManifoldPoint* mp = manifold->points + j;
			s2ContactConstraintPoint* cp = constraint->points + j;

			cp->normalImpulse = mp->normalImpulse * h;
			cp->tangentImpulse = mp->tangentImpulse * h;

			cp->localAnchorA = s2Sub(mp->localAnchorA, bodyA->localCenter);
			cp->localAnchorB = s2Sub(mp->localAnchorB, bodyB->localCenter);
			s2Vec2 rA = s2RotateVector(qA, cp->localAnchorA);
			s2Vec2 rB = s2RotateVector(qB, cp->localAnchorB);

			cp->rA0 = rA;
			cp->rB0 = rB;
			cp->separation = mp->separation;
			cp->adjustedSeparation = mp->separation - s2Dot(s2Sub(rB, rA), normal);

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

static void s2WarmStartContactPositions_XPBD(s2World* world, s2ContactConstraint* constraints, int constraintCount, float h)
{
	s2Body* bodies = world->bodies;
	float inv_h = h > 0.0f ? 1.0f / h : 0.0f;
	float warmStartCoefficient = 0.99;

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

		s2Vec2 dcA = bodyA->deltaPosition;
		s2Rot qA = bodyA->rot;
		s2Vec2 dcB = bodyB->deltaPosition;
		s2Rot qB = bodyB->rot;

		s2Vec2 normal = constraint->normal;
		s2Vec2 tangent = s2CrossVS(normal, 1.0f);

		// non-penetration constraints
		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			s2Vec2 rA = s2RotateVector(qA, cp->localAnchorA);
			s2Vec2 rB = s2RotateVector(qB, cp->localAnchorB);
			s2Vec2 drA = s2Sub(rA, cp->rA0);
			s2Vec2 drB = s2Sub(rB, cp->rB0);

			float rnA = s2Cross(rA, normal);
			float rnB = s2Cross(rB, normal);

			// w1 and w2 in paper
			float kA = mA + iA * rnA * rnA;
			float kB = mB + iB * rnB * rnB;

			// lambda has units of Mass * Length

			cp->normalImpulse *= warmStartCoefficient;
			s2Vec2 P = s2MulSV(cp->normalImpulse, normal);

			dcA = s2MulSub(dcA, mA, P);
			qA = s2IntegrateRot(qA, -iA * s2Cross(rA, P));

			dcB = s2MulAdd(dcB, mB, P);
			qB = s2IntegrateRot(qB, iB * s2Cross(rB, P));
		}

		// static friction constraints
		float friction = constraint->friction;

		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			s2Vec2 rA = s2RotateVector(qA, cp->localAnchorA);
			s2Vec2 rB = s2RotateVector(qB, cp->localAnchorB);
			s2Vec2 drA = s2Sub(rA, cp->rA0);
			s2Vec2 drB = s2Sub(rB, cp->rB0);

			// tangent separation
			s2Vec2 dp = s2Add(s2Sub(dcB, dcA), s2Sub(drB, drA));
			float C = s2Dot(dp, tangent);

			float rtA = s2Cross(rA, tangent);
			float rtB = s2Cross(rB, tangent);

			// w1 and w2 in paper
			float kA = mA + iA * rtA * rtA;
			float kB = mB + iB * rtB * rtB;

			cp->tangentImpulse *= warmStartCoefficient;
			s2Vec2 P = s2MulSV(cp->tangentImpulse, tangent);

			dcA = s2MulSub(dcA, mA, P);
			qA = s2IntegrateRot(qA, -iA * s2Cross(rA, P));

			dcB = s2MulAdd(dcB, mB, P);
			qB = s2IntegrateRot(qB, iB * s2Cross(rB, P));
		}

		bodyA->deltaPosition = dcA;
		bodyA->rot = qA;
		bodyB->deltaPosition = dcB;
		bodyB->rot = qB;
	}
}

static void s2SolveContactPositions_XPBD(s2World* world, s2ContactConstraint* constraints, int constraintCount, float h)
{
	s2Body* bodies = world->bodies;
	float inv_h = h > 0.0f ? 1.0f / h : 0.0f;

	// compliance because contacts are too energetic otherwise
	float compliance = 0.0001;
	float damping = 500.f;
	float beta = damping * h * h;
	float alpha = compliance * inv_h * inv_h;
	float gamma = alpha * beta * inv_h;

	// but the rush sample has too much overlap ...
	// float baseCompliance = 0.0f;

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

		s2Vec2 dcA = bodyA->deltaPosition;
		s2Rot qA = bodyA->rot;
		s2Vec2 dcB = bodyB->deltaPosition;
		s2Rot qB = bodyB->rot;

		s2Vec2 normal = constraint->normal;
		s2Vec2 tangent = s2CrossVS(normal, 1.0f);

		// non-penetration constraints
		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			s2Vec2 rA = s2RotateVector(qA, cp->localAnchorA);
			s2Vec2 rB = s2RotateVector(qB, cp->localAnchorB);
			s2Vec2 drA = s2Sub(rA, cp->rA0);
			s2Vec2 drB = s2Sub(rB, cp->rB0);

			// change in separation
			s2Vec2 ds = s2Add(s2Sub(dcB, dcA), s2Sub(drB, drA));
			float C = (1.0f + gamma) * s2Dot(ds, normal) + cp->separation + alpha * cp->normalImpulse;

			// this clamping is not in the paper, but it is used in other solvers
			//C = S2_MAX(-s2_maxBaumgarteVelocity * h, C);

			float rnA = s2Cross(rA, normal);
			float rnB = s2Cross(rB, normal);

			// w1 and w2 in paper
			float kA = mA + iA * rnA * rnA;
			float kB = mB + iB * rnB * rnB;

			// lambda has units of Mass * Length

			float lambda = S2_MAX(0.f, cp->normalImpulse - C / ((1.f + gamma) * (kA + kB) + alpha));
			float deltaLambda = lambda - cp->normalImpulse;
			cp->normalImpulse = lambda;

			s2Vec2 P = s2MulSV(deltaLambda, normal);

			dcA = s2MulSub(dcA, mA, P);
			qA = s2IntegrateRot(qA, -iA * s2Cross(rA, P));

			dcB = s2MulAdd(dcB, mB, P);
			qB = s2IntegrateRot(qB, iB * s2Cross(rB, P));
		}

		// static friction constraints
		float friction = constraint->friction;

		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			s2Vec2 rA = s2RotateVector(qA, cp->localAnchorA);
			s2Vec2 rB = s2RotateVector(qB, cp->localAnchorB);
			s2Vec2 drA = s2Sub(rA, cp->rA0);
			s2Vec2 drB = s2Sub(rB, cp->rB0);

			// tangent separation
			s2Vec2 dp = s2Add(s2Sub(dcB, dcA), s2Sub(drB, drA));
			float C = s2Dot(dp, tangent);

			float rtA = s2Cross(rA, tangent);
			float rtB = s2Cross(rB, tangent);

			// w1 and w2 in paper
			float kA = mA + iA * rtA * rtA;
			float kB = mB + iB * rtB * rtB;

			float lambda = cp->tangentImpulse - C / (kA + kB);
			float maxLambda = friction * cp->normalImpulse;

#if 0
			if (lambda < -maxLambda || maxLambda < lambda)
			{
				lambda = 0.0f;
			}
#else
			// this seems to behave better, but does not follow the paper
			lambda = S2_CLAMP(lambda, -maxLambda, maxLambda);
#endif

			float deltaLambda = lambda - cp->tangentImpulse;
			cp->tangentImpulse = lambda;

			s2Vec2 P = s2MulSV(deltaLambda, tangent);

			dcA = s2MulSub(dcA, mA, P);
			qA = s2IntegrateRot(qA, -iA * s2Cross(rA, P));

			dcB = s2MulAdd(dcB, mB, P);
			qB = s2IntegrateRot(qB, iB * s2Cross(rB, P));
		}

		bodyA->deltaPosition = dcA;
		bodyA->rot = qA;
		bodyB->deltaPosition = dcB;
		bodyB->rot = qB;
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

		s2Rot qA = bodyA->rot;
		s2Rot qB = bodyB->rot;

		s2Vec2 vA = bodyA->linearVelocity;
		float wA = bodyA->angularVelocity;
		s2Vec2 vB = bodyB->linearVelocity;
		float wB = bodyB->angularVelocity;

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
			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, rB));
			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, rA));
			s2Vec2 dv = s2Sub(vrB, vrA);

			float rnA = s2Cross(rA, normal);
			float rnB = s2Cross(rB, normal);

			// w1 and w2 in paper
			float kA = mA + iA * rnA * rnA;
			float kB = mB + iB * rnB * rnB;

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

			float rtA = s2Cross(rA, tangent);
			float rtB = s2Cross(rB, tangent);

			// w1 and w2 in paper
			float kA = mA + iA * rtA * rtA;
			float kB = mB + iB * rtB * rtB;

			// eq 31
			float maxFrictionImpulse = friction * cp->normalImpulse;

			// Length / Time (this is wrong in the paper, fixed here)
			float huf = (maxFrictionImpulse * inv_h) * (kA + kB);

			// Length / Time
			float abs_vt = S2_ABS(vt);

			float Cdot = (vt / abs_vt) * S2_MIN(huf, abs_vt);
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

	float h = context->dt / substepCount;
	float inv_h = 1.0f / h;

	// Loops
	// body: 1 + 2 * substepCount
	// constraint: 2 + 2 * substepCount

	// constraint loop
	s2PrepareContacts_XPBD(world, constraints, constraintCount, h);

	for (int i = 0; i < jointCapacity; ++i)
	{
		s2Joint* joint = joints + i;
		if (s2IsFree(&joint->object))
		{
			continue;
		}
		s2PrepareJoint_XPBD(joint, context);
	}

	s2Body* bodies = world->bodies;
	int bodyCapacity = world->bodyPool.capacity;
	s2Vec2 gravity = world->gravity;

	// body 2 * substepCount
	// constraint 2 * substepCount
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

			// store previous rotation
			body->rot0 = body->rot;

			// integrate positions
			// this is unique to XPBD, no other solvers update position immediately
			body->deltaPosition0 = body->deltaPosition;
			body->deltaPosition = s2MulAdd(body->deltaPosition, h, v);
			body->rot = s2IntegrateRot(body->rot, h * w);
		}

		s2WarmStartContactPositions_XPBD(world, constraints, constraintCount, h);

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

			body->linearVelocity = s2MulSV(inv_h, s2Sub(body->deltaPosition, body->deltaPosition0));

			if (s2Length(body->linearVelocity) > 10.0f)
			{
				body->linearVelocity.x += 0.0f;
			}

			body->angularVelocity = s2ComputeAngularVelocity(body->rot0, body->rot, inv_h);
		}

		// Relax contact velocities
		s2SolveContactVelocities_XPBD(world, constraints, constraintCount, h);
	}

	// Finalize body position
	// body loop
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

		body->position = s2Add(body->position, body->deltaPosition);
		body->deltaPosition = s2Vec2_zero;
	}

	// constraint loop
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
