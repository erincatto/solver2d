// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "graph.h"

#include "allocate.h"
#include "array.h"
#include "body.h"
#include "contact.h"
#include "core.h"
#include "shape.h"
#include "solver_data.h"
#include "stack_allocator.h"
#include "world.h"

#include "solver2d/aabb.h"

#include <stdbool.h>
#include <stdlib.h>

#define maxBaumgarteVelocity 3.0f

static void s2IntegrateVelocities(s2World* world, float h)
{
	s2Body* bodies = world->bodies;
	int bodyCapacity = world->bodyPool.capacity;
	s2Vec2 gravity = world->gravity;

	// Integrate velocities and apply damping. Initialize the body state.
	for (int i = 0; i < bodyCapacity; ++i)
	{
		s2Body* body = bodies + i;
		if (s2ObjectIsFree(&body->object))
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

		// Integrate velocities
		v = s2Add(v, s2MulSV(h * invMass, s2MulAdd(body->force, body->gravityScale * body->mass, gravity)));
		w = w + h * invI * body->torque;

		// Apply damping.
		// ODE: dv/dt + c * v = 0
		// Solution: v(t) = v0 * exp(-c * t)
		// Time step: v(t + dt) = v0 * exp(-c * (t + dt)) = v0 * exp(-c * t) * exp(-c * dt) = v * exp(-c * dt)
		// v2 = exp(-c * dt) * v1
		// Pade approximation:
		// v2 = v1 * 1 / (1 + c * dt)
		v = s2MulSV(1.0f / (1.0f + h * body->linearDamping), v);
		w *= 1.0f / (1.0f + h * body->angularDamping);

		body->linearVelocity = v;
		body->angularVelocity = w;

		body->deltaAngle = 0.0f;
		body->deltaPosition = s2Vec2_zero;
	}
}

#if 0 // no need?
static void s2IntegrateVelocitiesSoft(s2World* world, float h)
{
	s2Body* bodies = world->bodies;
	int bodyCapacity = world->bodyPool.capacity;
	s2Vec2 gravity = world->gravity;

	// Integrate velocities and apply damping. Initialize the body state.
	for (int i = 0; i < bodyCapacity; ++i)
	{
		s2Body* body = bodies + i;
		if (s2ObjectValid(&body->object) == false)
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

		// Integrate velocities
		v = s2Add(v, s2MulSV(h * invMass, s2MulAdd(body->force, body->gravityScale * body->mass, gravity)));
		w = w + h * invI * body->torque;

		// Apply damping.
		// ODE: dv/dt + c * v = 0
		// Solution: v(t) = v0 * exp(-c * t)
		// Time step: v(t + dt) = v0 * exp(-c * (t + dt)) = v0 * exp(-c * t) * exp(-c * dt) = v * exp(-c * dt)
		// v2 = exp(-c * dt) * v1
		// Pade approximation:
		// v2 = v1 * 1 / (1 + c * dt)
		v = s2MulSV(1.0f / (1.0f + h * body->linearDamping), v);
		w *= 1.0f / (1.0f + h * body->angularDamping);

		body->linearVelocity = v;
		body->angularVelocity = w;

		body->deltaAngle = 0.0f;
		body->deltaPosition = s2Vec2_zero;
	}
}
#endif

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

typedef struct s2ConstraintPoint
{
	s2Vec2 rA, rB;
	s2Vec2 rAf, rBf;
	s2Vec2 localAnchorA, localAnchorB;
	float tangentSeparation;
	float separation;
	float normalImpulse;
	float tangentImpulse;
	float normalMass;
	float tangentMass;
	float gamma;
	float massCoefficient;
	float biasCoefficient;
	float impulseCoefficient;
	float baumgarte;
	bool frictionValid;
} s2ConstraintPoint;

typedef struct s2Constraint
{
	s2Contact* contact;
	int indexA;
	int indexB;
	s2ConstraintPoint points[2];
	s2Vec2 normal;
	float friction;
	int pointCount;
} s2Constraint;

static void s2InitializeSoftConstraints(s2World* world, float h, s2Constraint* constraints, int constraintCount)
{
	s2Body* bodies = world->bodies;

	for (int i = 0; i < constraintCount; ++i)
	{
		s2Constraint* constraint = constraints + i;

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
			s2ConstraintPoint* cp = constraint->points + j;

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

static void s2InitializePGSConstraints(s2World* world, s2Constraint* constraints, int constraintCount)
{
	s2Body* bodies = world->bodies;

	for (int i = 0; i < constraintCount; ++i)
	{
		s2Constraint* constraint = constraints + i;

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
			s2ConstraintPoint* cp = constraint->points + j;

			cp->normalImpulse = mp->normalImpulse;
			cp->tangentImpulse = mp->tangentImpulse;

			cp->rA = s2Sub(mp->point, cA);
			cp->rB = s2Sub(mp->point, cB);
			cp->localAnchorA = s2InvRotateVector(qA, cp->rA);
			cp->localAnchorB = s2InvRotateVector(qB, cp->rB);
			cp->separation = mp->separation;

			cp->baumgarte = 0.0f;
			cp->biasCoefficient = mp->separation > 0.0f ? 1.0f : 0.0f;

			float rtA = s2Cross(cp->rA, tangent);
			float rtB = s2Cross(cp->rB, tangent);
			float kTangent = mA + mB + iA * rtA * rtA + iB * rtB * rtB;
			cp->tangentMass = kTangent > 0.0f ? 1.0f / kTangent : 0.0f;

			float rnA = s2Cross(cp->rA, normal);
			float rnB = s2Cross(cp->rB, normal);
			float kNormal = mA + mB + iA * rnA * rnA + iB * rnB * rnB;
			cp->normalMass = kNormal > 0.0f ? 1.0f / kNormal : 0.0f;
		}
	}
}

static void s2InitializeStickyConstraints(s2World* world, s2Constraint* constraints, int constraintCount)
{
	s2Body* bodies = world->bodies;

	for (int i = 0; i < constraintCount; ++i)
	{
		s2Constraint* constraint = constraints + i;

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
			s2ConstraintPoint* cp = constraint->points + j;

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
				s2ConstraintPoint* cp = constraint->points + j;

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
				s2ConstraintPoint* cp = constraint->points + j;

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

static void s2WarmStart(s2World* world, s2Constraint* constraints, int constraintCount)
{
	s2Body* bodies = world->bodies;

	for (int i = 0; i < constraintCount; ++i)
	{
		s2Constraint* constraint = constraints + i;

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
			s2ConstraintPoint* cp = constraint->points + j;

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

static void s2WarmStartAll(s2World* world, s2Constraint* constraints, int constraintCount)
{
	s2Body* bodies = world->bodies;

	for (int i = 0; i < constraintCount; ++i)
	{
		s2Constraint* constraint = constraints + i;

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
			s2ConstraintPoint* cp = constraint->points + j;

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

static void s2SolveVelocityConstraints(s2World* world, s2Constraint* constraints, int constraintCount, float inv_dt)
{
	s2Body* bodies = world->bodies;

	for (int i = 0; i < constraintCount; ++i)
	{
		s2Constraint* constraint = constraints + i;

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
			s2ConstraintPoint* cp = constraint->points + j;

			// Relative velocity at contact
			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, cp->rB));
			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, cp->rA));
			s2Vec2 dv = s2Sub(vrB, vrA);

			// Compute normal impulse
			float vn = s2Dot(dv, normal);
			float impulse = -cp->normalMass * (vn + cp->biasCoefficient * cp->separation * inv_dt);

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
			s2ConstraintPoint* cp = constraint->points + j;

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

static void s2SolveVelocityConstraintsSorted(s2World* world, s2Constraint* constraints, int constraintCount, float inv_dt)
{
	s2Body* bodies = world->bodies;

	for (int i = 0; i < constraintCount; ++i)
	{
		s2Constraint* constraint = constraints + i;

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
			s2ConstraintPoint* cp = constraint->points + j;

			// Relative velocity at contact
			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, cp->rB));
			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, cp->rA));
			s2Vec2 dv = s2Sub(vrB, vrA);

			// Compute normal impulse
			float vn = s2Dot(dv, normal);
			float impulse = -cp->normalMass * (vn + cp->biasCoefficient * cp->separation * inv_dt);

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
			s2ConstraintPoint* cp = constraint->points + j;

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

static void s2SolveVelocityConstraintsSoft(s2World* world, s2Constraint* constraints, int constraintCount, float inv_dt, bool useBias)
{
	s2Body* bodies = world->bodies;

	for (int i = 0; i < constraintCount; ++i)
	{
		s2Constraint* constraint = constraints + i;

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
			s2ConstraintPoint* cp = constraint->points + j;

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
			s2ConstraintPoint* cp = constraint->points + j;

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

static void s2SolveVelocityConstraintsSticky(s2World* world, s2Constraint* constraints, int constraintCount, float minSeparation,
											 float invh)
{
	s2Body* bodies = world->bodies;

	for (int i = 0; i < constraintCount; ++i)
	{
		s2Constraint* constraint = constraints + i;

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
			s2ConstraintPoint* cp = constraint->points + j;

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
			s2ConstraintPoint* cp = constraint->points + j;

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

static void s2StoreImpulses(s2Constraint* constraints, int constraintCount)
{
	for (int i = 0; i < constraintCount; ++i)
	{
		s2Constraint* constraint = constraints + i;
		s2Contact* contact = constraint->contact;

		s2Manifold* manifold = &contact->manifold;

		for (int j = 0; j < constraint->pointCount; ++j)
		{
			manifold->points[j].normalImpulse = constraint->points[j].normalImpulse;
			manifold->points[j].tangentImpulse = constraint->points[j].tangentImpulse;
		}
	}
}

static void s2IntegratePositions(s2World* world, float h)
{
	s2Body* bodies = world->bodies;
	int bodyCapacity = world->bodyPool.capacity;

	// Integrate velocities and apply damping. Initialize the body state.
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

		s2Vec2 c = body->position;
		float a = body->angle;
		s2Vec2 v = body->linearVelocity;
		float w = body->angularVelocity;

		// Clamp large velocities
		s2Vec2 translation = s2MulSV(h, v);
		if (s2Dot(translation, translation) > s2_maxTranslation * s2_maxTranslation)
		{
			float ratio = s2_maxTranslation / s2Length(translation);
			v = s2MulSV(ratio, v);
		}

		float rotation = h * w;
		if (rotation * rotation > s2_maxRotation * s2_maxRotation)
		{
			float ratio = s2_maxRotation / S2_ABS(rotation);
			w *= ratio;
		}

		// Integrate
		c = s2MulAdd(c, h, v);
		a += h * w;

		body->position = c;
		body->angle = a;
		body->linearVelocity = v;
		body->angularVelocity = w;
	}
}

static void s2SolvePositionConstraints(s2World* world, s2Constraint* constraints, int constraintCount)
{
	s2Body* bodies = world->bodies;
	float slop = s2_linearSlop;

	for (int i = 0; i < constraintCount; ++i)
	{
		s2Constraint* constraint = constraints + i;

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
			s2ConstraintPoint* cp = constraint->points + j;

			s2Rot qA = s2MakeRot(aA);
			s2Rot qB = s2MakeRot(aB);

			s2Vec2 rA = s2RotateVector(qA, cp->localAnchorA);
			s2Vec2 rB = s2RotateVector(qB, cp->localAnchorB);

			// Current separation
			s2Vec2 d = s2Sub(s2Add(cB, rB), s2Add(cA, rA));
			float separation = s2Dot(d, normal) + cp->separation;

			// Prevent large corrections. Need to maintain a small overlap to avoid overshoot.
			// This improves stacking stability significantly.
			float C = S2_CLAMP(s2_baumgarte * (separation + slop), -s2_maxLinearCorrection, 0.0f);

			// Compute the effective mass.
			float rnA = s2Cross(rA, normal);
			float rnB = s2Cross(rB, normal);
			float K = mA + mB + iA * rnA * rnA + iB * rnB * rnB;

			// Compute normal impulse
			float impulse = K > 0.0f ? -C / K : 0.0f;

			s2Vec2 P = s2MulSV(impulse, normal);

			cA = s2MulSub(cA, mA, P);
			aA -= iA * s2Cross(cp->rA, P);

			cB = s2MulAdd(cB, mB, P);
			aB += iB * s2Cross(cp->rB, P);
		}

		bodyA->position = cA;
		bodyA->angle = aA;
		bodyB->position = cB;
		bodyB->angle = aB;
	}
}

static void s2SolvePositionConstraintsSorted(s2World* world, s2Constraint* constraints, int constraintCount)
{
	s2Body* bodies = world->bodies;

	for (int i = 0; i < constraintCount; ++i)
	{
		s2Constraint* constraint = constraints + i;

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
		float slop = s2_linearSlop;

		for (int j = 0; j < pointCount; ++j)
		{
			s2ConstraintPoint* cp = constraint->points + j;

			s2Rot qA = s2MakeRot(aA);
			s2Rot qB = s2MakeRot(aB);

			s2Vec2 rA = s2RotateVector(qA, cp->localAnchorA);
			s2Vec2 rB = s2RotateVector(qB, cp->localAnchorB);

			// Current separation
			s2Vec2 d = s2Sub(s2Add(cB, rB), s2Add(cA, rA));
			float separation = s2Dot(d, normal) + cp->separation;

			// Prevent large corrections. Need to maintain a small overlap to avoid overshoot.
			// This improves stacking stability significantly.
			float C = S2_CLAMP(s2_baumgarte * (separation + slop), -s2_maxLinearCorrection, 0.0f);

			// Compute the effective mass.
			float rnA = s2Cross(rA, normal);
			float rnB = s2Cross(rB, normal);
			float K = mA + mB + iA * rnA * rnA + iB * rnB * rnB;

			// Compute normal impulse
			float impulse = K > 0.0f ? -C / K : 0.0f;

			s2Vec2 P = s2MulSV(impulse, normal);

			cA = s2MulSub(cA, mA, P);
			aA -= iA * s2Cross(cp->rA, P);

			cB = s2MulAdd(cB, mB, P);
			aB += iB * s2Cross(cp->rB, P);
		}

		bodyA->position = cA;
		bodyA->angle = aA;
		bodyB->position = cB;
		bodyB->angle = aB;
	}
}

// Update body transform, mark broadphase AABB, build awake contact bits
static void s2FinalizeSolve(s2World* world)
{
	s2Body* bodies = world->bodies;
	int bodyCapacity = world->bodyPool.capacity;

	const s2Vec2 aabbMargin = {s2_aabbMargin, s2_aabbMargin};

	// Integrate velocities and apply damping. Initialize the body state.
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

		body->transform.q = s2MakeRot(body->angle);
		body->transform.p = s2Sub(body->position, s2RotateVector(body->transform.q, body->localCenter));

		body->force = s2Vec2_zero;
		body->torque = 0.0f;

		// Update shapes AABBs
		int shapeIndex = body->shapeList;
		while (shapeIndex != S2_NULL_INDEX)
		{
			s2Shape* shape = world->shapes + shapeIndex;

			shape->aabb = s2Shape_ComputeAABB(shape, body->transform);

			if (s2AABB_Contains(shape->fatAABB, shape->aabb) == false)
			{
				shape->fatAABB.lowerBound = s2Sub(shape->aabb.lowerBound, aabbMargin);
				shape->fatAABB.upperBound = s2Add(shape->aabb.upperBound, aabbMargin);
				shape->enlargedAABB = true;
			}

			shapeIndex = shape->nextShapeIndex;
		}
	}
}

void s2SolveGraphPGS(s2World* world, const s2StepContext* stepContext)
{
	s2Contact* contacts = world->contacts;
	int contactCapacity = world->contactPool.capacity;

	s2Constraint* constraints = s2AllocateStackItem(world->stackAllocator, contactCapacity * sizeof(s2Constraint), "constraint");
	
	int constraintCount = 0;
	for (int i = 0; i < contactCapacity; ++i)
	{
		s2Contact* contact = contacts + i;
		if (s2ObjectIsFree(&contact->object))
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

	int velocityIterations = stepContext->velocityIterations;
	int positionIterations = stepContext->positionIterations;
	float h = stepContext->dt;
	float inv_h = stepContext->inv_dt;

	s2IntegrateVelocities(world, h);
	s2InitializePGSConstraints(world, constraints, constraintCount);

	s2WarmStartAll(world, constraints, constraintCount);

	for (int iter = 0; iter < velocityIterations; ++iter)
	{
		s2SolveVelocityConstraints(world, constraints, constraintCount, inv_h);
	}

	s2StoreImpulses(constraints, constraintCount);

	s2IntegratePositions(world, h);

	for (int iter = 0; iter < positionIterations; ++iter)
	{
		s2SolvePositionConstraints(world, constraints, constraintCount);
	}

	s2FinalizeSolve(world);
	
	s2FreeStackItem(world->stackAllocator, constraints);
}

#if 0
void s2SolveGraphSoftPGS(s2World* world, const s2StepContext* stepContext)
{
	s2Graph* graph = &world->graph;
	s2GraphColor* colors = graph->colors;

	int constraintCount = 0;
	for (int i = 0; i < s2_graphColorCount; ++i)
	{
		constraintCount += s2Array(colors[i].contactArray).count;
	}

	s2Constraint* constraints = s2AllocateStackItem(world->stackAllocator, constraintCount * sizeof(s2Constraint), "constraint");
	int base = 0;

	for (int i = 0; i < s2_graphColorCount; ++i)
	{
		colors[i].constraints = constraints + base;
		base += s2Array(colors[i].contactArray).count;
	}

	S2_ASSERT(base == constraintCount);

	int velocityIterations = stepContext->velocityIterations;
	int positionIterations = stepContext->positionIterations;
	float h = stepContext->dt;

	s2IntegrateVelocities(world, h);

	for (int i = 0; i < s2_graphColorCount; ++i)
	{
		s2InitializeSoftConstraints(world, colors + i, h, true);
	}

	for (int iter = 0; iter < velocityIterations; ++iter)
	{
		for (int i = 0; i < s2_graphColorCount; ++i)
		{
			s2SolveVelocityConstraintsSoft(world, colors + i, stepContext->inv_dt, true);
		}
	}
	
	s2IntegratePositions(world, h);

	for (int iter = 0; iter < positionIterations; ++iter)
	{
		for (int i = 0; i < s2_graphColorCount; ++i)
		{
			s2SolveVelocityConstraintsSoft(world, colors + i, stepContext->inv_dt, false);
		}
	}

	s2StoreImpulses(constraints, constraintCount);

	s2FinalizeSolve(world);

	s2FreeStackItem(world->stackAllocator, constraints);
}

void s2SolveGraphSoftTGS(s2World* world, const s2StepContext* stepContext)
{
	s2Graph* graph = &world->graph;
	s2GraphColor* colors = graph->colors;

	int constraintCount = 0;
	for (int i = 0; i < s2_graphColorCount; ++i)
	{
		constraintCount += s2Array(colors[i].contactArray).count;
	}

	s2Constraint* constraints = s2AllocateStackItem(world->stackAllocator, constraintCount * sizeof(s2Constraint), "constraint");
	int base = 0;

	for (int i = 0; i < s2_graphColorCount; ++i)
	{
		colors[i].constraints = constraints + base;
		base += s2Array(colors[i].contactArray).count;
	}

	S2_ASSERT(base == constraintCount);

	// Full step apply gravity
	s2IntegrateVelocities(world, stepContext->dt);

	for (int i = 0; i < s2_graphColorCount; ++i)
	{
		bool warmStart = true;
		s2InitializeSoftConstraints(world, colors + i, stepContext->dt, warmStart);
	}
	
	int substepCount = stepContext->velocityIterations;
	float h = stepContext->dt / substepCount;
	float inv_h = 1.0f / h;

	for (int substep = 0; substep < substepCount; ++substep)
	{
		// One constraint iteration
		for (int i = 0; i < s2_graphColorCount; ++i)
		{
			bool useBias = true;
			s2SolveVelocityConstraintsSoft(world, colors + i, inv_h, useBias);
		}

		s2IntegrateDeltaTransform(world, h);
	}

	s2UpdatePositions(world);

	int positionIterations = stepContext->positionIterations;
	for (int iter = 0; iter < positionIterations; ++iter)
	{
		bool useBias = false;
		s2SolveVelocityConstraintsSoft(world, 0.0f, useBias);
	}

	s2StoreImpulses(constraints, constraintCount);

	s2FinalizeSolve(world);

	s2FreeStackItem(world->stackAllocator, constraints);
}
#endif

// Sticky
void s2SolveGraphStickyTGS(s2World* world, const s2StepContext* stepContext)
{
	s2Contact* contacts = world->contacts;
	int contactCapacity = world->contactPool.capacity;

	s2Constraint* constraints = s2AllocateStackItem(world->stackAllocator, contactCapacity * sizeof(s2Constraint), "constraint");

	int constraintCount = 0;
	for (int i = 0; i < contactCapacity; ++i)
	{
		s2Contact* contact = contacts + i;
		if (s2ObjectIsFree(&contact->object))
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

	s2IntegrateVelocities(world, stepContext->dt);

	s2InitializeStickyConstraints(world, constraints, constraintCount);

	int substepCount = stepContext->velocityIterations;
	float h = stepContext->dt / substepCount;
	float invh = substepCount / stepContext->dt;

	for (int substep = 0; substep < substepCount; ++substep)
	{
		// One constraint iteration
		s2SolveVelocityConstraintsSticky(world, constraints, constraintCount, -s2_huge, invh);
		s2IntegrateDeltaTransform(world, h);
	}

	s2UpdatePositions(world);

	int positionIterations = stepContext->positionIterations;
	for (int iter = 0; iter < positionIterations; ++iter)
	{
		// relax constraints
		s2SolveVelocityConstraintsSticky(world, constraints, constraintCount, 0.0f, 0.0f);
	}

	s2FinalizeSolve(world);

	s2FreeStackItem(world->stackAllocator, constraints);
}
