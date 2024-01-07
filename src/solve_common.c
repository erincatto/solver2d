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

	// Integrate velocities and apply damping. Initialize the body state.
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

		// Integrate velocities
		v = s2Add(v, s2MulSV(h * invMass, s2MulAdd(body->force, body->mass, gravity)));
		w = w + h * invI * body->torque;

		body->linearVelocity = v;
		body->angularVelocity = w;

		body->deltaAngle = 0.0f;
		body->deltaPosition = s2Vec2_zero;
	}
}

void s2IntegratePositions(s2World* world, float h)
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

		s2Vec2 normal = constraint->normal;
		s2Vec2 tangent = s2RightPerp(normal);

		for (int j = 0; j < pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;

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
