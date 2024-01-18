// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "body.h"
#include "contact.h"
#include "core.h"
#include "solvers.h"
#include "world.h"

void s2IntegrateVelocities(s2World* world, float h)
{
	s2Body* bodies = world->bodies;
	int bodyCapacity = world->bodyPool.capacity;
	s2Vec2 gravity = world->gravity;

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

		float invMass = body->invMass;
		float invI = body->invI;

		s2Vec2 v = body->linearVelocity;
		float w = body->angularVelocity;

		v = s2Add(v, s2MulSV(h * invMass, s2MulAdd(body->force, body->mass * body->gravityScale, gravity)));
		w = w + h * invI * body->torque;

		// Damping
		v = s2MulSV(1.0f / (1.0f + h * body->linearDamping), v);
		w *= 1.0f / (1.0f + h * body->angularDamping);

		body->linearVelocity = v;
		body->angularVelocity = w;
	}
}

void s2IntegratePositions(s2World* world, float h)
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

		body->position = s2MulAdd(body->position, h, body->linearVelocity);
		body->transform.q = s2IntegrateRot(body->transform.q, h * body->angularVelocity);
	}
}

// Soft contact math
// float d = 2.0f * zeta * omega / kNormal;
// float k = omega * omega / kNormal;
// cp->gamma = 1.0f / (h * (d + h * k));
// cp->gamma = 1.0f / (h * (2.0f * zeta * omega / kNormal + h * omega * omega / kNormal));
// cp->gamma = kNormal / (h * omega * (2.0f * zeta + h * omega));
// cp->bias = h * k * cp->gamma * mp->separation;
// cp->bias = k / (d + h * k) * mp->separation;
// cp->bias =
//	(omega * omega / kNormal) / (2 * zeta * omega / kNormal + h * omega * omega / kNormal) * mp->separation;
// cp->gamma = 0.0f;
// cp->bias = (0.2f / h) * mp->separation;
// This can be expanded
// cp->normalMass = 1.0f / (kNormal + cp->gamma);
// meff = 1.0f / kNormal * 1.0f / (1.0f + 1.0f / (h * omega * (2 * zeta + h * omega)))
// float impulse = -cp->normalMass * (vn + bias + cp->gamma * cp->normalImpulse);
// = -meff * mscale * (vn + bias) - imp_scale * impulse

void s2PrepareContacts_Soft(s2World* world, s2ContactConstraint* constraints, int constraintCount, s2StepContext* context, float h, float hertz)
{
	s2Body* bodies = world->bodies;
	bool warmStart = context->warmStart;

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

		// Stiffer for dynamic vs static
		float contactHertz = (mA == 0.0f || mB == 0.0f) ? 2.0f * hertz : hertz;

		s2Vec2 cA = bodyA->position;
		s2Vec2 cB = bodyB->position;
		s2Rot qA = s2MakeRot(bodyA->angle);
		s2Rot qB = s2MakeRot(bodyB->angle);

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

			// static anchors
			cp->rAs = rA;
			cp->rBs = rB;

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

			// soft contact
			// should use the substep not the full time step
			const float zeta = 10.0f;
			float omega = 2.0f * s2_pi * contactHertz;
			float c = h * omega * (2.0f * zeta + h * omega);
			cp->biasCoefficient = omega / (2.0f * zeta + h * omega);
			cp->impulseCoefficient = 1.0f / (1.0f + c);
			cp->massCoefficient = c * cp->impulseCoefficient;
		}
	}
}

void s2WarmStartContacts(s2World* world, s2ContactConstraint* constraints, int constraintCount)
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

		s2Rot qA = s2MakeRot(bodyA->angle);
		s2Rot qB = s2MakeRot(bodyB->angle);

		s2Vec2 normal = constraint->normal;
		s2Vec2 tangent = s2RightPerp(normal);

		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			// Current anchors (to support TGS)
			s2Vec2 rA = s2RotateVector(qA, cp->localAnchorA);
			s2Vec2 rB = s2RotateVector(qB, cp->localAnchorB);

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

void s2SolveContact_NGS(s2World* world, s2ContactConstraint* constraints, int constraintCount, float fraction)
{
	s2Body* bodies = world->bodies;
	float slop = s2_linearSlop;

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

		s2Vec2 normal = constraint->normal;

		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			s2Rot qA = s2MakeRot(aA);
			s2Rot qB = s2MakeRot(aB);

			s2Vec2 rA = s2RotateVector(qA, cp->localAnchorA);
			s2Vec2 rB = s2RotateVector(qB, cp->localAnchorB);

			// Current separation
			s2Vec2 d = s2Sub(s2Add(cB, rB), s2Add(cA, rA));
			float separation = s2Dot(d, normal) + cp->separation;

			// Prevent large corrections. Need to maintain a small overlap to avoid overshoot.
			// This improves stacking stability significantly.
			float C = S2_CLAMP(separation + slop, -s2_maxLinearCorrection, 0.0f);

			// Compute the effective mass.
			float rnA = s2Cross(rA, normal);
			float rnB = s2Cross(rB, normal);
			float K = mA + mB + iA * rnA * rnA + iB * rnB * rnB;

			// Compute normal impulse
			float impulse = K > 0.0f ? -fraction * s2_baumgarte * C / K : 0.0f;

			s2Vec2 P = s2MulSV(impulse, normal);

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

void s2StoreContactImpulses(s2ContactConstraint* constraints, int constraintCount)
{
	for (int i = 0; i < constraintCount; ++i)
	{
		s2ContactConstraint* constraint = constraints + i;
		s2Contact* contact = constraint->contact;
		s2Manifold* manifold = &contact->manifold;

		for (int j = 0; j < constraint->pointCount; ++j)
		{
			manifold->points[j].normalImpulse = constraint->points[j].normalImpulse;
			manifold->points[j].tangentImpulse = constraint->points[j].tangentImpulse;
		}
	}
}
