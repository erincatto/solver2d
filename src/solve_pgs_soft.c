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

#define maxBaumgarteVelocity 3.0f

#if 0
static void s2IntegrateDeltaTransform(s2World* world, float h)
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

static void s2UpdatePositions(s2World* world)
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
#endif

static void s2InitializeSoftConstraints(s2World* world, float h, s2ContactConstraint* constraints, int constraintCount)
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
			cp->gamma = kNormal / (h * omega * (2.0f * zeta + h * omega));

			cp->separation = mp->separation;

			// cp->bias = h * k * cp->gamma * mp->separation;
			// cp->bias = k / (d + h * k) * mp->separation;
			// cp->bias =
			//	(omega * omega / kNormal) / (2 * zeta * omega / kNormal + h * omega * omega / kNormal) * mp->separation;
			cp->biasCoefficient = omega / (2.0f * zeta + h * omega);
			// cp->gamma = 0.0f;
			// cp->bias = (0.2f / h) * mp->separation;

			// TODO_ERIN this can be expanded
			cp->normalMass = kNormal > 0.0f ? 1.0f / kNormal : 0.0f;
			//cp->normalMass = 1.0f / (kNormal + cp->gamma);

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

#if 0
static void s2InitializeStickyConstraints(s2World* world, s2ContactConstraint* constraints, int constraintCount)
{
	s2Body* bodies = world->bodies;

	for (int i = 0; i < constraintCount; ++i)
	{
		s2ContactConstraint* constraint = constraints + i;

		s2Contact* contact = constraint->contact;
		s2Manifold* manifold = &contact->manifold;
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

			cp->rA = s2Sub(mp->point, cA);
			cp->rB = s2Sub(mp->point, cB);
			cp->localAnchorA = s2InvRotateVector(qA, cp->rA);
			cp->localAnchorB = s2InvRotateVector(qB, cp->rB);
			cp->separation = mp->separation;

			cp->baumgarte = 0.8f;

			float rtA = s2Cross(cp->rA, tangent);
			float rtB = s2Cross(cp->rB, tangent);
			float kTangent = mA + mB + iA * rtA * rtA + iB * rtB * rtB;
			cp->tangentMass = kTangent > 0.0f ? 1.0f / kTangent : 0.0f;

			float rnA = s2Cross(cp->rA, normal);
			float rnB = s2Cross(cp->rB, normal);
			float kNormal = mA + mB + iA * rnA * rnA + iB * rnB * rnB;
			cp->normalMass = kNormal > 0.0f ? 1.0f / kNormal : 0.0f;
		}

		bool frictionConfirmed = false;
		if (manifold->frictionPersisted)
		{
			int confirmCount = 0;
			for (int j = 0; j < pointCount; ++j)
			{
				const s2ManifoldPoint* mp = manifold->points + j;
				s2ContactConstraintPoint* cp = constraint->points + j;

				s2Vec2 normalA = s2RotateVector(qA, mp->localNormalA);
				s2Vec2 normalB = s2RotateVector(qB, mp->localNormalB);

				float nn = s2Dot(normalA, normalB);
				if (nn < 0.98f)
				{
					// Relative rotation has invalidated cached friction anchors
					break;
				}

				s2Vec2 anchorA = s2RotateVector(qA, mp->localAnchorA);
				s2Vec2 anchorB = s2RotateVector(qB, mp->localAnchorB);
				s2Vec2 offset = s2Add(s2Sub(cB, cA), s2Sub(anchorB, anchorA));
				float normalSeparation = s2Dot(offset, normalA);
				if (S2_ABS(normalSeparation) > 2.0f * s2_linearSlop)
				{
					// Normal separation has invalidated cached friction anchors
					break;
				}

				cp->rAf = anchorA;
				cp->rBf = anchorB;
				cp->tangentSeparation = s2Dot(offset, tangent);

				float rtA = s2Cross(anchorA, tangent);
				float rtB = s2Cross(anchorB, tangent);
				float kTangent = mA + mB + iA * rtA * rtA + iB * rtB * rtB;
				cp->tangentMass = kTangent > 0.0f ? 1.0f / kTangent : 0.0f;

				confirmCount += 1;
			}

			if (confirmCount == pointCount)
			{
				frictionConfirmed = true;
			}
		}

		if (frictionConfirmed == false)
		{
			for (int j = 0; j < pointCount; ++j)
			{
				s2ManifoldPoint* mp = manifold->points + j;
				s2ContactConstraintPoint* cp = constraint->points + j;

				mp->localNormalA = s2InvRotateVector(qA, normal);
				mp->localNormalB = s2InvRotateVector(qB, normal);
				mp->localAnchorA = s2InvRotateVector(qA, cp->rA);
				mp->localAnchorB = s2InvRotateVector(qB, cp->rB);

				cp->rAf = cp->rA;
				cp->rBf = cp->rB;
				cp->tangentSeparation = 0.0f;

				float rtA = s2Cross(cp->rAf, tangent);
				float rtB = s2Cross(cp->rBf, tangent);
				float kTangent = mA + mB + iA * rtA * rtA + iB * rtB * rtB;
				cp->tangentMass = kTangent > 0.0f ? 1.0f / kTangent : 0.0f;
			}
		}

		manifold->frictionPersisted = true;
	}
}
#endif

static void s2SolveVelocityConstraintsSoft(s2World* world, s2ContactConstraint* constraints, int constraintCount, float inv_dt, bool useBias)
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

			// Relative velocity at contact
			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, cp->rB));
			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, cp->rA));
			s2Vec2 dv = s2Sub(vrB, vrA);

			// Compute change in separation
			s2Vec2 prB = s2Add(dpB, s2CrossSV(daB, cp->rB));
			s2Vec2 prA = s2Add(dpA, s2CrossSV(daA, cp->rA));
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
				//bias = cp->biasCoefficient * s;
				massScale = cp->massCoefficient;
				impulseScale = cp->impulseCoefficient;
			}
			
			// Compute normal impulse
			float vn = s2Dot(dv, normal);
			float impulse = -cp->normalMass * massScale * (vn + bias) - impulseScale * cp->normalImpulse;
			//float impulse = -cp->normalMass * (vn + bias + cp->gamma * cp->normalImpulse);

			// Clamp the accumulated impulse
			float newImpulse = S2_MAX(cp->normalImpulse + impulse, 0.0f);
			impulse = newImpulse - cp->normalImpulse;
			cp->normalImpulse = newImpulse;

			// Apply contact impulse
			s2Vec2 P = s2MulSV(impulse, normal);
			vA = s2MulSub(vA, mA, P);
			wA -= iA * s2Cross(cp->rA, P);

			vB = s2MulAdd(vB, mB, P);
			wB += iB * s2Cross(cp->rB, P);
		}

		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			// Relative velocity at contact
			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, cp->rB));
			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, cp->rA));
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
			wA -= iA * s2Cross(cp->rA, P);

			vB = s2MulAdd(vB, mB, P);
			wB += iB * s2Cross(cp->rB, P);
		}

		bodyA->linearVelocity = vA;
		bodyA->angularVelocity = wA;
		bodyB->linearVelocity = vB;
		bodyB->angularVelocity = wB;
	}
}

#if 0
static void s2SolveVelocityConstraintsSticky(s2World* world, s2ContactConstraint* constraints, int constraintCount, float minSeparation,
											 float invh)
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
		float friction = 0.3f; //constraint->friction;

		float totalNormalImpulse = 0.0f;

		// Non-penetration constraints
		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			// Relative velocity at contact
			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, cp->rB));
			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, cp->rA));
			s2Vec2 dv = s2Sub(vrB, vrA);

			// Compute change in separation
			s2Vec2 prB = s2Add(dpB, s2CrossSV(daB, cp->rB));
			s2Vec2 prA = s2Add(dpA, s2CrossSV(daA, cp->rA));
			float ds = s2Dot(s2Sub(prB, prA), normal);
			float s = cp->separation + ds;

			float bias = 0.0f;
			if (s > 0.0f)
			{
				// Speculative
				bias = s * invh;
			
			}
			else if (minSeparation < 0.0f)
			{
				bias = S2_MAX(-maxBaumgarteVelocity, cp->baumgarte * s * invh);
			}

			// Compute normal impulse
			float vn = s2Dot(dv, normal);
			float impulse = -cp->normalMass * (vn + bias);

			// Clamp the accumulated impulse
			float newImpulse = S2_MAX(cp->normalImpulse + impulse, 0.0f);
			impulse = newImpulse - cp->normalImpulse;
			cp->normalImpulse = newImpulse;

			totalNormalImpulse += cp->normalImpulse;

			// Apply contact impulse
			s2Vec2 P = s2MulSV(impulse, normal);
			vA = s2MulSub(vA, mA, P);
			wA -= iA * s2Cross(cp->rA, P);

			vB = s2MulAdd(vB, mB, P);
			wB += iB * s2Cross(cp->rB, P);
		}

		// Sticky friction constraints
		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

			// Relative velocity at contact
			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, cp->rBf));
			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, cp->rAf));
			s2Vec2 dv = s2Sub(vrB, vrA);

			// Compute change in separation
			s2Vec2 prB = s2Add(dpB, s2CrossSV(daB, cp->rBf));
			s2Vec2 prA = s2Add(dpA, s2CrossSV(daA, cp->rAf));
			float ds = s2Dot(s2Sub(prB, prA), tangent);
			float s = cp->tangentSeparation + ds;
			float bias = 0.5f * s * invh;

			// Compute tangent impulse
			float vt = s2Dot(dv, tangent);
			float impulse = -cp->tangentMass * (vt + bias);

			// max friction uses an average of the total normal impulse because persistent friction anchors don't line up with normal anchors
			float maxFriction = 0.5f * friction * totalNormalImpulse;

			// Clamp the accumulated impulse
			float newImpulse = cp->tangentImpulse + impulse;
			if (newImpulse < -maxFriction)
			{
				newImpulse = -maxFriction;
				constraint->contact->manifold.frictionPersisted = false;
			}
			else if (newImpulse > maxFriction)
			{
				newImpulse = maxFriction;
				constraint->contact->manifold.frictionPersisted = false;
			}

			impulse = newImpulse - cp->tangentImpulse;
			cp->tangentImpulse = newImpulse;

			// Apply contact impulse
			s2Vec2 P = s2MulSV(impulse, tangent);

			vA = s2MulSub(vA, mA, P);
			wA -= iA * s2Cross(cp->rA, P);

			vB = s2MulAdd(vB, mB, P);
			wB += iB * s2Cross(cp->rB, P);
		}

		bodyA->linearVelocity = vA;
		bodyA->angularVelocity = wA;
		bodyB->linearVelocity = vB;
		bodyB->angularVelocity = wB;
	}
}
#endif

void s2Solve_PGS_Soft(s2World* world, s2StepContext* context)
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

	int velocityIterations = context->velocityIterations;
	int positionIterations = context->positionIterations;
	float h = context->dt;

	s2IntegrateVelocities(world, h);

	s2InitializeSoftConstraints(world, h, constraints, constraintCount);

	for (int i = 0; i < jointCapacity; ++i)
	{
		s2Joint* joint = joints + i;
		if (s2IsFree(&joint->object))
		{
			continue;
		}
		s2PrepareJoint(joint, context);
	}

	for (int iter = 0; iter < velocityIterations; ++iter)
	{
		for (int i = 0; i < jointCapacity; ++i)
		{
			s2Joint* joint = joints + i;
			if (s2IsFree(&joint->object))
			{
				continue;
			}

			s2SolveJointVelocity(joint, context);
		}

		s2SolveVelocityConstraintsSoft(world, constraints, constraintCount, context->inv_dt, true);
	}
	
	s2IntegratePositions(world, h);

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

		s2SolveVelocityConstraintsSoft(world, constraints, constraintCount, context->inv_dt, false);
	}

	s2StoreContactImpulses(constraints, constraintCount);

	s2FreeStackItem(world->stackAllocator, constraints);
}

#if 0
// Sticky
void s2SolveGraphStickyTGS(s2World* world, const s2StepContext* context)
{
	s2Contact* contacts = world->contacts;
	int contactCapacity = world->contactPool.capacity;

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

	s2IntegrateVelocities(world, context->dt);

	s2InitializeStickyConstraints(world, constraints, constraintCount);

	int substepCount = context->velocityIterations;
	float h = context->dt / substepCount;
	float invh = substepCount / context->dt;

	for (int substep = 0; substep < substepCount; ++substep)
	{
		// One constraint iteration
		s2SolveVelocityConstraintsSticky(world, constraints, constraintCount, -s2_huge, invh);
		s2IntegrateDeltaTransform(world, h);
	}

	s2UpdatePositions(world);

	int positionIterations = context->positionIterations;
	for (int iter = 0; iter < positionIterations; ++iter)
	{
		// relax constraints
		s2SolveVelocityConstraintsSticky(world, constraints, constraintCount, 0.0f, 0.0f);
	}

	s2FinalizeSolve(world);

	s2FreeStackItem(world->stackAllocator, constraints);
}
#endif
