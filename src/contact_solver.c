// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "contact_solver.h"

#include "array.h"
#include "body.h"
#include "contact.h"
#include "core.h"
#include "stack_allocator.h"
#include "world.h"

// Solver debugging is normally disabled because the block solver sometimes has to deal with a poorly conditioned
// effective mass matrix.
#define S2_DEBUG_SOLVER 0

typedef struct s2VelocityConstraintPoint
{
	s2Vec2 rA;
	s2Vec2 rB;
	float normalImpulse;
	float tangentImpulse;
	float normalMass;
	float tangentMass;
	float velocityBias;
	float relativeVelocity;
} s2VelocityConstraintPoint;

typedef struct s2ContactVelocityConstraint
{
	s2Contact* contact;
	s2VelocityConstraintPoint points[2];
	s2Vec2 normal;
	s2Mat22 normalMass;
	s2Mat22 K;
	float friction;
	float restitution;
	float tangentSpeed;
	int32_t pointCount;
} s2ContactVelocityConstraint;

typedef struct s2ContactPositionConstraint
{
	s2Contact* contact;
	s2Vec2 localAnchorsA[2];
	s2Vec2 localAnchorsB[2];
	float separations[2];
	float lambdas[2];
	s2Vec2 normal;
	int32_t pointCount;
} s2ContactPositionConstraint;

s2ContactSolver* s2CreateContactSolver(s2ContactSolverDef* def)
{
	s2StackAllocator* alloc = def->world->stackAllocator;

	s2ContactSolver* solver = s2AllocateStackItem(alloc, sizeof(s2ContactSolver), "contact solver");
	solver->context = def->context;
	solver->contactList = def->contactList;
	solver->contactCount = def->contactCount;

	// These are allocated conservatively because some island contacts may not have contact points
	solver->positionConstraints =
		s2AllocateStackItem(alloc, solver->contactCount * sizeof(s2ContactPositionConstraint), "position constraints");
	solver->velocityConstraints =
		s2AllocateStackItem(alloc, solver->contactCount * sizeof(s2ContactVelocityConstraint), "velocity constraints");

	solver->world = def->world;
	solver->constraintCount = 0;
	return solver;
}

void s2ContactSolver_Initialize(s2ContactSolver* solver)
{
	s2World* world = solver->world;
	s2Contact* contacts = world->contacts;
	const s2StepContext* context = solver->context;
	s2Body* bodies = world->bodies;

	// Initialize position independent portions of the constraints.
	int32_t constraintCount = 0;
	int32_t contactIndex = solver->contactList;
	while (contactIndex != S2_NULL_INDEX)
	{
		s2Contact* contact = contacts + contactIndex;

		const s2Manifold* manifold = &contact->manifold;
		int32_t pointCount = manifold->pointCount;

		if (pointCount == 0)
		{
			continue;
		}

		int32_t indexA = contact->edges[0].bodyIndex;
		int32_t indexB = contact->edges[1].bodyIndex;
		s2Body* bodyA = bodies + indexA;
		s2Body* bodyB = bodies + indexB;

		s2ContactVelocityConstraint* vc = solver->velocityConstraints + constraintCount;
		vc->contact = contact;
		vc->normal = manifold->normal;
		vc->friction = contact->friction;
		vc->restitution = contact->restitution;
		vc->pointCount = pointCount;
		vc->K = s2Mat22_zero;
		vc->normalMass = s2Mat22_zero;

		s2ContactPositionConstraint* pc = solver->positionConstraints + constraintCount;
		pc->contact = contact;
		pc->normal = manifold->normal;
		pc->pointCount = pointCount;

		float mA = bodyA->invMass;
		float iA = bodyA->invI;
		float mB = bodyB->invMass;
		float iB = bodyB->invI;

		s2Rot qA = bodyA->transform.q;
		s2Vec2 cA = bodyA->position;
		s2Rot qB = bodyB->transform.q;
		s2Vec2 cB = bodyB->position;

		// TODO_ERIN testing
		// qA = s2MakeRot(bodyA->angle);
		// qB = s2MakeRot(bodyB->angle);

		s2Vec2 vA = bodyA->linearVelocity;
		float wA = bodyA->angularVelocity;
		s2Vec2 vB = bodyB->linearVelocity;
		float wB = bodyB->angularVelocity;

		for (int32_t j = 0; j < pointCount; ++j)
		{
			const s2ManifoldPoint* cp = manifold->points + j;
			s2VelocityConstraintPoint* vcp = vc->points + j;

			vcp->normalImpulse = cp->normalImpulse;
			vcp->tangentImpulse = cp->tangentImpulse;

			vcp->rA = s2Sub(cp->point, cA);
			vcp->rB = s2Sub(cp->point, cB);

			float rnA = s2Cross(vcp->rA, vc->normal);
			float rnB = s2Cross(vcp->rB, vc->normal);

			float kNormal = mA + mB + iA * rnA * rnA + iB * rnB * rnB;

			vcp->normalMass = kNormal > 0.0f ? 1.0f / kNormal : 0.0f;

			s2Vec2 tangent = s2CrossVS(vc->normal, 1.0f);

			float rtA = s2Cross(vcp->rA, tangent);
			float rtB = s2Cross(vcp->rB, tangent);

			float kTangent = mA + mB + iA * rtA * rtA + iB * rtB * rtB;

			vcp->tangentMass = kTangent > 0.0f ? 1.0f / kTangent : 0.0f;

			// Velocity bias for speculative collision
			vcp->velocityBias = -S2_MAX(0.0f, cp->separation * context->inv_dt);

			// Relative velocity
			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, vcp->rB));
			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, vcp->rA));
			vcp->relativeVelocity = s2Dot(vc->normal, s2Sub(vrB, vrA));

			pc->localAnchorsA[j] = s2InvRotateVector(qA, vcp->rA);
			pc->localAnchorsB[j] = s2InvRotateVector(qB, vcp->rB);
			pc->separations[j] = cp->separation;
			pc->lambdas[j] = 0.0f;
		}

		// If we have two points, then prepare the block solver.
		if (vc->pointCount == 2)
		{
			s2VelocityConstraintPoint* vcp1 = vc->points + 0;
			s2VelocityConstraintPoint* vcp2 = vc->points + 1;

			float rn1A = s2Cross(vcp1->rA, vc->normal);
			float rn1B = s2Cross(vcp1->rB, vc->normal);
			float rn2A = s2Cross(vcp2->rA, vc->normal);
			float rn2B = s2Cross(vcp2->rB, vc->normal);

			float k11 = mA + mB + iA * rn1A * rn1A + iB * rn1B * rn1B;
			float k22 = mA + mB + iA * rn2A * rn2A + iB * rn2B * rn2B;
			float k12 = mA + mB + iA * rn1A * rn2A + iB * rn1B * rn2B;

			// Ensure a reasonable condition number.
			const float k_maxConditionNumber = 1000.0f;
			if (k11 * k11 < k_maxConditionNumber * (k11 * k22 - k12 * k12))
			{
				// K is safe to invert.
				vc->K.cx = (s2Vec2){k11, k12};
				vc->K.cy = (s2Vec2){k12, k22};
				vc->normalMass = s2GetInverse22(vc->K);
			}
			else
			{
				// The constraints are redundant, just use one.
				// TODO_ERIN use deepest?
				vc->pointCount = 1;
			}
		}

		constraintCount += 1;
	}

	solver->constraintCount = constraintCount;

	for (int32_t i = 0; i < constraintCount; ++i)
	{
		s2ContactVelocityConstraint* vc = solver->velocityConstraints + i;

		const s2Contact* contact = vc->contact;

		int32_t indexA = contact->edges[0].bodyIndex;
		int32_t indexB = contact->edges[1].bodyIndex;
		s2Body* bodyA = bodies + indexA;
		s2Body* bodyB = bodies + indexB;
		float mA = bodyA->invMass;
		float iA = bodyA->invI;
		float mB = bodyB->invMass;
		float iB = bodyB->invI;
		int32_t pointCount = vc->pointCount;

		s2Vec2 vA = bodyA->linearVelocity;
		float wA = bodyA->angularVelocity;
		s2Vec2 vB = bodyB->linearVelocity;
		float wB = bodyB->angularVelocity;

		s2Vec2 normal = vc->normal;
		s2Vec2 tangent = s2CrossVS(normal, 1.0f);

		for (int32_t j = 0; j < pointCount; ++j)
		{
			s2VelocityConstraintPoint* vcp = vc->points + j;
			s2Vec2 P = s2Add(s2MulSV(vcp->normalImpulse, normal), s2MulSV(vcp->tangentImpulse, tangent));
			wA -= iA * s2Cross(vcp->rA, P);
			vA = s2MulAdd(vA, -mA, P);
			wB += iB * s2Cross(vcp->rB, P);
			vB = s2MulAdd(vB, mB, P);
		}

		bodyA->linearVelocity = vA;
		bodyA->angularVelocity = wA;
		bodyB->linearVelocity = vB;
		bodyB->angularVelocity = wB;
	}
}

void s2ContactSolver_SolveVelocityConstraints(s2ContactSolver* solver)
{
	int32_t count = solver->constraintCount;

	s2World* world = solver->world;
	s2Body* bodies = world->bodies;

	for (int32_t i = 0; i < count; ++i)
	{
		s2ContactVelocityConstraint* vc = solver->velocityConstraints + i;

		const s2Contact* contact = vc->contact;

		int32_t indexA = contact->edges[0].bodyIndex;
		int32_t indexB = contact->edges[1].bodyIndex;
		s2Body* bodyA = bodies + indexA;
		s2Body* bodyB = bodies + indexB;

		float mA = bodyA->invMass;
		float iA = bodyA->invI;
		float mB = bodyB->invMass;
		float iB = bodyB->invI;
		int32_t pointCount = vc->pointCount;

		s2Vec2 vA = bodyA->linearVelocity;
		float wA = bodyA->angularVelocity;
		s2Vec2 vB = bodyB->linearVelocity;
		float wB = bodyB->angularVelocity;

		s2Vec2 normal = vc->normal;
		s2Vec2 tangent = s2CrossVS(normal, 1.0f);
		float friction = vc->friction;

		S2_ASSERT(pointCount == 1 || pointCount == 2);

		// Solve tangent constraints first because non-penetration is more important
		// than friction.
		for (int32_t j = 0; j < pointCount; ++j)
		{
			s2VelocityConstraintPoint* vcp = vc->points + j;

			// Relative velocity at contact
			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, vcp->rB));
			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, vcp->rA));
			s2Vec2 dv = s2Sub(vrB, vrA);

			// Compute tangent force
			float vt = s2Dot(dv, tangent) - vc->tangentSpeed;
			float lambda = vcp->tangentMass * (-vt);

			// Clamp the accumulated force
			float maxFriction = friction * vcp->normalImpulse;
			float newImpulse = S2_CLAMP(vcp->tangentImpulse + lambda, -maxFriction, maxFriction);
			lambda = newImpulse - vcp->tangentImpulse;
			vcp->tangentImpulse = newImpulse;

			// Apply contact impulse
			s2Vec2 P = s2MulSV(lambda, tangent);

			vA = s2MulSub(vA, mA, P);
			wA -= iA * s2Cross(vcp->rA, P);

			vB = s2MulAdd(vB, mB, P);
			wB += iB * s2Cross(vcp->rB, P);
		}

		// Solve normal constraints
		if (pointCount == 1)
		{
			for (int32_t j = 0; j < pointCount; ++j)
			{
				s2VelocityConstraintPoint* vcp = vc->points + j;

				// Relative velocity at contact
				s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, vcp->rB));
				s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, vcp->rA));
				s2Vec2 dv = s2Sub(vrB, vrA);

				// Compute normal impulse
				float vn = s2Dot(dv, normal);
				float lambda = -vcp->normalMass * (vn - vcp->velocityBias);

				// Clamp the accumulated impulse
				float newImpulse = S2_MAX(vcp->normalImpulse + lambda, 0.0f);
				lambda = newImpulse - vcp->normalImpulse;
				vcp->normalImpulse = newImpulse;

				// Apply contact impulse
				s2Vec2 P = s2MulSV(lambda, normal);
				vA = s2MulSub(vA, mA, P);
				wA -= iA * s2Cross(vcp->rA, P);

				vB = s2MulAdd(vB, mB, P);
				wB += iB * s2Cross(vcp->rB, P);
			}
		}
		else
		{
			// Block solver developed in collaboration with Dirk Gregorius (back in 01/07 on Box2D_Lite).
			// Build the mini LCP for this contact patch
			//
			// vn = A * x + b, vn >= 0, x >= 0 and vn_i * x_i = 0 with i = 1..2
			//
			// A = J * W * JT and J = ( -n, -r1 x n, n, r2 x n )
			// b = vn0 - velocityBias
			//
			// The system is solved using the "Total enumeration method" (s. Murty). The complementary constraint vn_i *
			// x_i implies that we must have in any solution either vn_i = 0 or x_i = 0. So for the 2D contact problem
			// the cases vn1 = 0 and vn2 = 0, x1 = 0 and x2 = 0, x1 = 0 and vn2 = 0, x2 = 0 and vn1 = 0 need to be
			// tested. The first valid solution that satisfies the problem is chosen.
			//
			// In order to account of the accumulated impulse 'a' (because of the iterative nature of the solver which
			// only requires that the accumulated impulse is clamped and not the incremental impulse) we change the
			// impulse variable (x_i).
			//
			// Substitute:
			//
			// x = a + d
			//
			// a := old total impulse
			// x := new total impulse
			// d := incremental impulse
			//
			// For the current iteration we extend the formula for the incremental impulse
			// to compute the new total impulse:
			//
			// vn = A * d + b
			//    = A * (x - a) + b
			//    = A * x + b - A * a
			//    = A * x + b'
			// b' = b - A * a;

			s2VelocityConstraintPoint* cp1 = vc->points + 0;
			s2VelocityConstraintPoint* cp2 = vc->points + 1;

			s2Vec2 a = {cp1->normalImpulse, cp2->normalImpulse};
			S2_ASSERT(a.x >= 0.0f && a.y >= 0.0f);

			// Relative velocity at contact
			s2Vec2 vrA, vrB;
			vrA = s2Add(vA, s2CrossSV(wA, cp1->rA));
			vrB = s2Add(vB, s2CrossSV(wB, cp1->rB));
			s2Vec2 dv1 = s2Sub(vrB, vrA);
			vrA = s2Add(vA, s2CrossSV(wA, cp2->rA));
			vrB = s2Add(vB, s2CrossSV(wB, cp2->rB));
			s2Vec2 dv2 = s2Sub(vrB, vrA);

			// Compute normal velocity
			float vn1 = s2Dot(dv1, normal);
			float vn2 = s2Dot(dv2, normal);

			s2Vec2 b = {vn1 - cp1->velocityBias, vn2 - cp2->velocityBias};

			// Compute b'
			b = s2Sub(b, s2MulMV(vc->K, a));

			const float k_errorTol = 1e-3f;
			S2_MAYBE_UNUSED(k_errorTol);

			for (;;)
			{
				//
				// Case 1: vn = 0
				//
				// 0 = A * x + b'
				//
				// Solve for x:
				//
				// x = - inv(A) * b'
				//
				s2Vec2 x = s2Neg(s2MulMV(vc->normalMass, b));

				if (x.x >= 0.0f && x.y >= 0.0f)
				{
					// Get the incremental impulse
					s2Vec2 d = s2Sub(x, a);

					// Apply incremental impulse
					s2Vec2 P1 = s2MulSV(d.x, normal);
					s2Vec2 P2 = s2MulSV(d.y, normal);
					vA = s2MulSub(vA, mA, s2Add(P1, P2));
					wA -= iA * (s2Cross(cp1->rA, P1) + s2Cross(cp2->rA, P2));

					vB = s2MulAdd(vB, mB, s2Add(P1, P2));
					wB += iB * (s2Cross(cp1->rB, P1) + s2Cross(cp2->rB, P2));

					// Accumulate
					cp1->normalImpulse = x.x;
					cp2->normalImpulse = x.y;

#if S2_DEBUG_SOLVER == 1
					// Postconditions
					dv1 = vB + s2Cross(wB, cp1->rB) - vA - s2Cross(wA, cp1->rA);
					dv2 = vB + s2Cross(wB, cp2->rB) - vA - s2Cross(wA, cp2->rA);

					// Compute normal velocity
					vn1 = s2Dot(dv1, normal);
					vn2 = s2Dot(dv2, normal);

					S2_ASSERT(s2Abs(vn1 - cp1->velocityBias) < k_errorTol);
					S2_ASSERT(s2Abs(vn2 - cp2->velocityBias) < k_errorTol);
#endif
					break;
				}

				//
				// Case 2: vn1 = 0 and x2 = 0
				//
				//   0 = a11 * x1 + a12 * 0 + b1'
				// vn2 = a21 * x1 + a22 * 0 + s2'
				//
				x.x = -cp1->normalMass * b.x;
				x.y = 0.0f;
				vn1 = 0.0f;
				vn2 = vc->K.cx.y * x.x + b.y;
				if (x.x >= 0.0f && vn2 >= 0.0f)
				{
					// Get the incremental impulse
					s2Vec2 d = s2Sub(x, a);

					// Apply incremental impulse
					s2Vec2 P1 = s2MulSV(d.x, normal);
					s2Vec2 P2 = s2MulSV(d.y, normal);

					vA = s2MulSub(vA, mA, s2Add(P1, P2));
					wA -= iA * (s2Cross(cp1->rA, P1) + s2Cross(cp2->rA, P2));

					vB = s2MulAdd(vB, mB, s2Add(P1, P2));
					wB += iB * (s2Cross(cp1->rB, P1) + s2Cross(cp2->rB, P2));

					// Accumulate
					cp1->normalImpulse = x.x;
					cp2->normalImpulse = x.y;

#if S2_DEBUG_SOLVER == 1
					// Postconditions
					dv1 = vB + s2Cross(wB, cp1->rB) - vA - s2Cross(wA, cp1->rA);

					// Compute normal velocity
					vn1 = s2Dot(dv1, normal);

					S2_ASSERT(s2Abs(vn1 - cp1->velocityBias) < k_errorTol);
#endif
					break;
				}

				//
				// Case 3: vn2 = 0 and x1 = 0
				//
				// vn1 = a11 * 0 + a12 * x2 + b1'
				//   0 = a21 * 0 + a22 * x2 + s2'
				//
				x.x = 0.0f;
				x.y = -cp2->normalMass * b.y;
				vn1 = vc->K.cy.x * x.y + b.x;
				vn2 = 0.0f;

				if (x.y >= 0.0f && vn1 >= 0.0f)
				{
					// Resubstitute for the incremental impulse
					s2Vec2 d = s2Sub(x, a);

					// Apply incremental impulse
					s2Vec2 P1 = s2MulSV(d.x, normal);
					s2Vec2 P2 = s2MulSV(d.y, normal);

					vA = s2MulSub(vA, mA, s2Add(P1, P2));
					wA -= iA * (s2Cross(cp1->rA, P1) + s2Cross(cp2->rA, P2));

					vB = s2MulAdd(vB, mB, s2Add(P1, P2));
					wB += iB * (s2Cross(cp1->rB, P1) + s2Cross(cp2->rB, P2));

					// Accumulate
					cp1->normalImpulse = x.x;
					cp2->normalImpulse = x.y;

#if S2_DEBUG_SOLVER == 1
					// Postconditions
					dv2 = vB + s2Cross(wB, cp2->rB) - vA - s2Cross(wA, cp2->rA);

					// Compute normal velocity
					vn2 = s2Dot(dv2, normal);

					S2_ASSERT(s2Abs(vn2 - cp2->velocityBias) < k_errorTol);
#endif
					break;
				}

				//
				// Case 4: x1 = 0 and x2 = 0
				//
				// vn1 = b1
				// vn2 = s2;
				x.x = 0.0f;
				x.y = 0.0f;
				vn1 = b.x;
				vn2 = b.y;

				if (vn1 >= 0.0f && vn2 >= 0.0f)
				{
					// Resubstitute for the incremental impulse
					s2Vec2 d = s2Sub(x, a);

					// Apply incremental impulse
					s2Vec2 P1 = s2MulSV(d.x, normal);
					s2Vec2 P2 = s2MulSV(d.y, normal);

					vA = s2MulSub(vA, mA, s2Add(P1, P2));
					wA -= iA * (s2Cross(cp1->rA, P1) + s2Cross(cp2->rA, P2));

					vB = s2MulAdd(vB, mB, s2Add(P1, P2));
					wB += iB * (s2Cross(cp1->rB, P1) + s2Cross(cp2->rB, P2));

					// Accumulate
					cp1->normalImpulse = x.x;
					cp2->normalImpulse = x.y;

					break;
				}

				// No solution, give up. This is hit sometimes, but it doesn't seem to matter.
				break;
			}
		}

		bodyA->linearVelocity = vA;
		bodyA->angularVelocity = wA;
		bodyB->linearVelocity = vB;
		bodyB->angularVelocity = wB;
	}
}

void s2ContactSolver_ApplyRestitution(s2ContactSolver* solver)
{
	int32_t count = solver->constraintCount;
	float threshold = solver->context->restitutionThreshold;

	s2World* world = solver->world;
	s2Body* bodies = world->bodies;

	for (int32_t i = 0; i < count; ++i)
	{
		s2ContactVelocityConstraint* vc = solver->velocityConstraints + i;
		const s2Contact* contact = vc->contact;

		int32_t indexA = contact->edges[0].bodyIndex;
		int32_t indexB = contact->edges[1].bodyIndex;
		s2Body* bodyA = bodies + indexA;
		s2Body* bodyB = bodies + indexB;

		if (vc->restitution == 0.0f)
		{
			continue;
		}

		float mA = bodyA->invMass;
		float iA = bodyA->invI;
		float mB = bodyB->invMass;
		float iB = bodyB->invI;
		int32_t pointCount = vc->pointCount;

		s2Vec2 vA = bodyA->linearVelocity;
		float wA = bodyA->angularVelocity;
		s2Vec2 vB = bodyB->linearVelocity;
		float wB = bodyB->angularVelocity;

		s2Vec2 normal = vc->normal;

		for (int32_t j = 0; j < pointCount; ++j)
		{
			s2VelocityConstraintPoint* vcp = vc->points + j;

			// if the normal impulse is zero then there was no collision
			if (vcp->relativeVelocity > -threshold || vcp->normalImpulse == 0.0f)
			{
				continue;
			}

			// Relative velocity at contact
			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, vcp->rB));
			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, vcp->rA));
			s2Vec2 dv = s2Sub(vrB, vrA);

			// Compute normal impulse
			float vn = s2Dot(dv, normal);
			float lambda = -vcp->normalMass * (vn + vc->restitution * vcp->relativeVelocity);

			// Apply contact impulse
			s2Vec2 P = s2MulSV(lambda, normal);
			vA = s2MulSub(vA, mA, P);
			wA -= iA * s2Cross(vcp->rA, P);

			vB = s2MulAdd(vB, mB, P);
			wB += iB * s2Cross(vcp->rB, P);
		}

		bodyA->linearVelocity = vA;
		bodyA->angularVelocity = wA;
		bodyB->linearVelocity = vB;
		bodyB->angularVelocity = wB;
	}
}

void s2ContactSolver_StoreImpulses(s2ContactSolver* solver)
{
	int32_t count = solver->constraintCount;

	for (int32_t i = 0; i < count; ++i)
	{
		s2ContactVelocityConstraint* vc = solver->velocityConstraints + i;
		s2Contact* contact = vc->contact;

		s2Manifold* manifold = &contact->manifold;

		for (int32_t j = 0; j < vc->pointCount; ++j)
		{
			manifold->points[j].normalImpulse = vc->points[j].normalImpulse;
			manifold->points[j].tangentImpulse = vc->points[j].tangentImpulse;
		}
	}
}

bool s2ContactSolver_SolvePositionConstraintsBlock(s2ContactSolver* solver)
{
	float minSeparation = 0.0f;
	int32_t count = solver->constraintCount;
	float slop = s2_linearSlop;

	s2World* world = solver->world;
	s2Body* bodies = world->bodies;

	for (int32_t i = 0; i < count; ++i)
	{
		s2ContactPositionConstraint* pc = solver->positionConstraints + i;
		const s2Contact* contact = pc->contact;

		int32_t indexA = contact->edges[0].bodyIndex;
		int32_t indexB = contact->edges[1].bodyIndex;
		s2Body* bodyA = bodies + indexA;
		s2Body* bodyB = bodies + indexB;

		float mA = bodyA->invMass;
		float iA = bodyA->invI;
		float mB = bodyB->invMass;
		float iB = bodyB->invI;

		int32_t pointCount = pc->pointCount;

		s2Vec2 cA = bodyA->position;
		float aA = bodyA->angle;
		s2Vec2 cB = bodyB->position;
		float aB = bodyB->angle;

		s2Vec2 normal = pc->normal;

		if (pointCount == 2)
		{
			s2Rot qA = s2MakeRot(aA);
			s2Rot qB = s2MakeRot(aB);

			s2Vec2 rA1 = s2RotateVector(qA, pc->localAnchorsA[0]);
			s2Vec2 rB1 = s2RotateVector(qB, pc->localAnchorsB[0]);
			s2Vec2 rA2 = s2RotateVector(qA, pc->localAnchorsA[1]);
			s2Vec2 rS2 = s2RotateVector(qB, pc->localAnchorsB[1]);

			// Current separation
			s2Vec2 d1 = s2Sub(s2Add(cB, rB1), s2Add(cA, rA1));
			float separation1 = s2Dot(d1, normal) + pc->separations[0];

			s2Vec2 d2 = s2Sub(s2Add(cB, rS2), s2Add(cA, rA2));
			float separation2 = s2Dot(d2, normal) + pc->separations[1];

			// Track max constraint error.
			minSeparation = S2_MIN(minSeparation, separation1);
			minSeparation = S2_MIN(minSeparation, separation2);

			float C1 = S2_CLAMP(s2_baumgarte * (separation1 + slop), -s2_maxLinearCorrection, 0.0f);
			float C2 = S2_CLAMP(s2_baumgarte * (separation2 + slop), -s2_maxLinearCorrection, 0.0f);

			s2Vec2 b = {C1, C2};

			float rn1A = s2Cross(rA1, normal);
			float rn1B = s2Cross(rB1, normal);
			float rn2A = s2Cross(rA2, normal);
			float rn2B = s2Cross(rS2, normal);

			float k11 = mA + mB + iA * rn1A * rn1A + iB * rn1B * rn1B;
			float k22 = mA + mB + iA * rn2A * rn2A + iB * rn2B * rn2B;
			float k12 = mA + mB + iA * rn1A * rn2A + iB * rn1B * rn2B;

			s2Mat22 K, invK;

			// Ensure a reasonable condition number.
			const float k_maxConditionNumber = 10000.0f;
			if (k11 * k11 < k_maxConditionNumber * (k11 * k22 - k12 * k12))
			{
				// K is safe to invert.
				K.cx = (s2Vec2){k11, k12};
				K.cy = (s2Vec2){k12, k22};
				invK = s2GetInverse22(K);
			}
			else
			{
				// The constraints are redundant, however one may be deeper than the other.
				// This can happen when a capsule is deeply embedded in a box.
				goto manifold_degenerate;
			}

			const float k_errorTol = 1e-3f;
			S2_MAYBE_UNUSED(k_errorTol);

			for (;;)
			{
				//
				// Case 1: vn = 0
				//
				// 0 = A * x + b'
				//
				// Solve for x:
				//
				// x = - inv(A) * b'
				//
				s2Vec2 x = s2Neg(s2MulMV(invK, b));

				if (x.x >= 0.0f && x.y >= 0.0f)
				{
					// Get the incremental impulse
					s2Vec2 d = x;

					// Apply incremental impulse
					s2Vec2 P1 = s2MulSV(d.x, normal);
					s2Vec2 P2 = s2MulSV(d.y, normal);

					cA = s2MulSub(cA, mA, s2Add(P1, P2));
					aA -= iA * (s2Cross(rA1, P1) + s2Cross(rA2, P2));

					cB = s2MulAdd(cB, mB, s2Add(P1, P2));
					aB += iB * (s2Cross(rB1, P1) + s2Cross(rS2, P2));
					break;
				}

				//
				// Case 2: vn1 = 0 and x2 = 0
				//
				//   0 = a11 * x1 + a12 * 0 + b1'
				// vn2 = a21 * x1 + a22 * 0 + s2'
				//
				x.x = -b.x / k11;
				x.y = 0.0f;
				float vn2 = K.cx.y * x.x + b.y;
				if (x.x >= 0.0f && vn2 >= 0.0f)
				{
					// Get the incremental impulse
					s2Vec2 d = x;

					// Apply incremental impulse
					s2Vec2 P1 = s2MulSV(d.x, normal);
					s2Vec2 P2 = s2MulSV(d.y, normal);

					cA = s2MulSub(cA, mA, s2Add(P1, P2));
					aA -= iA * (s2Cross(rA1, P1) + s2Cross(rA2, P2));

					cB = s2MulAdd(cB, mB, s2Add(P1, P2));
					aB += iB * (s2Cross(rB1, P1) + s2Cross(rS2, P2));
					break;
				}

				//
				// Case 3: vn2 = 0 and x1 = 0
				//
				// vn1 = a11 * 0 + a12 * x2 + b1'
				//   0 = a21 * 0 + a22 * x2 + s2'
				//
				x.x = 0.0f;
				x.y = -b.y / k22;
				float vn1 = K.cy.x * x.y + b.x;
				if (x.y >= 0.0f && vn1 >= 0.0f)
				{
					// Resubstitute for the incremental impulse
					s2Vec2 d = x;

					// Apply incremental impulse
					s2Vec2 P1 = s2MulSV(d.x, normal);
					s2Vec2 P2 = s2MulSV(d.y, normal);

					cA = s2MulSub(cA, mA, s2Add(P1, P2));
					aA -= iA * (s2Cross(rA1, P1) + s2Cross(rA2, P2));

					cB = s2MulAdd(cB, mB, s2Add(P1, P2));
					aB += iB * (s2Cross(rB1, P1) + s2Cross(rS2, P2));
					break;
				}
				break;
			}
		}
		else
		{
		manifold_degenerate:
			for (int32_t j = 0; j < pointCount; ++j)
			{
				s2Rot qA = s2MakeRot(aA);
				s2Rot qB = s2MakeRot(aB);

				s2Vec2 rA = s2RotateVector(qA, pc->localAnchorsA[j]);
				s2Vec2 rB = s2RotateVector(qB, pc->localAnchorsB[j]);

				// Current separation
				s2Vec2 d = s2Sub(s2Add(cB, rB), s2Add(cA, rA));
				float separation = s2Dot(d, normal) + pc->separations[j];

				// Track max constraint error.
				minSeparation = S2_MIN(minSeparation, separation);

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
				aA -= iA * s2Cross(rA, P);

				cB = s2MulAdd(cB, mB, P);
				aB += iB * s2Cross(rB, P);
			}
		}

		bodyA->position = cA;
		bodyA->angle = aA;
		bodyB->position = cB;
		bodyB->angle = aB;
	}

	// We can't expect minSpeparation >= -s2_linearSlop because we don't
	// push the separation above -s2_linearSlop.
	return minSeparation >= -3.0f * s2_linearSlop;
}
