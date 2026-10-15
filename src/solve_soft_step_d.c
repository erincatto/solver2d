// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "body.h"
#include "contact.h"
#include "core.h"
#include "joint.h"
#include "joint_factor.h"
#include "solvers.h"
#include "stack_allocator.h"
#include "world.h"

#include <stdbool.h>

// Soft Step with the point rows of revolute joints solved together by a direct factorization. An isolated
// joint matches Soft Step. Joints sharing a body are solved exactly instead of by one Gauss-Seidel pass.
// Soft Step's per joint update dlambda = -m K^-1 (Cdot + b) - i lambda becomes, over all joints,
// (D + G) dlambda = -(J v + b) - G lambda with D = J M^-1 J^T and compliance G = K / c per joint.

// Jacobian anchors follow the bodies each substep. With 0 they stay at the start of the step, as Soft Step
// does for contacts. Both refactor before every direct solve for now.
#define S2_SOFT_STEP_D_REFACTOR 1

// Solve every revolute joint directly, ignoring the per joint opt in
#define S2_SOFT_STEP_D_ALL_DIRECT 1

// Relax with the soft compliance kept on the left. Reuses the soft matrix but is inexact on long chains.
#define S2_SOFT_STEP_D_RELAX_COMPLIANCE 0

// Added to the unit diagonal of the scaled relax block, which redundant joints make singular.
// Well above float roundoff so impulses that move nothing stay bounded.
#define S2_SOFT_STEP_D_RELAX_REGULARIZATION 1.0e-5f

static bool s2UseDirect(s2World* world, s2Joint* joint)
{
	if (joint->type != s2_revoluteJoint)
	{
		return false;
	}

#if S2_SOFT_STEP_D_ALL_DIRECT == 0
	if (joint->revoluteJoint.enableDirect == false)
	{
		return false;
	}
#endif

	s2Body* bodyA = world->bodies + joint->edges[0].bodyIndex;
	s2Body* bodyB = world->bodies + joint->edges[1].bodyIndex;
	return bodyA->type == s2_dynamicBody || bodyB->type == s2_dynamicBody;
}

static void s2UpdateDirectAnchors(s2World* world, s2JointFactor* factor)
{
	for (int i = 0; i < factor->jointCount; ++i)
	{
		s2RevoluteJoint* joint = &world->joints[factor->jointIndices[i]].revoluteJoint;
		s2PointRows* rows = factor->rows + i;
		rows->rA = s2RotateVector(world->bodies[rows->indexA].rot, joint->localAnchorA);
		rows->rB = s2RotateVector(world->bodies[rows->indexB].rot, joint->localAnchorB);
	}
}

static void s2FactorDirect(s2World* world, s2JointFactor* factor, bool soft, float regularization)
{
#if S2_SOFT_STEP_D_REFACTOR == 1
	s2UpdateDirectAnchors(world, factor);
#endif

	for (int i = 0; i < factor->jointCount; ++i)
	{
		s2RevoluteJoint* joint = &world->joints[factor->jointIndices[i]].revoluteJoint;
		factor->compliance[i] = soft ? joint->impulseCoefficient / joint->massCoefficient : 0.0f;
	}

	s2FactorJoints(world, factor, regularization);
}

static void s2WarmStartDirect(s2World* world, s2JointFactor* factor)
{
	for (int i = 0; i < factor->jointCount; ++i)
	{
		s2RevoluteJoint* joint = &world->joints[factor->jointIndices[i]].revoluteJoint;
		const s2PointRows* rows = factor->rows + i;
		s2Body* bodyA = world->bodies + rows->indexA;
		s2Body* bodyB = world->bodies + rows->indexB;

		float mA = joint->invMassA, mB = joint->invMassB;
		float iA = joint->invIA, iB = joint->invIB;

		float axialImpulse = joint->motorImpulse + joint->lowerImpulse - joint->upperImpulse;
		s2Vec2 P = joint->impulse;

		bodyA->linearVelocity = s2MulSub(bodyA->linearVelocity, mA, P);
		bodyA->angularVelocity -= iA * (s2Cross(rows->rA, P) + axialImpulse);
		bodyB->linearVelocity = s2MulAdd(bodyB->linearVelocity, mB, P);
		bodyB->angularVelocity += iB * (s2Cross(rows->rB, P) + axialImpulse);
	}
}

// Point rows of all direct joints in one solve against the current factor
static void s2SolveDirect(s2World* world, s2JointFactor* factor, bool useBias)
{
	s2Body* bodies = world->bodies;
	int jointCount = factor->jointCount;

	for (int i = 0; i < jointCount; ++i)
	{
		s2RevoluteJoint* joint = &world->joints[factor->jointIndices[i]].revoluteJoint;
		const s2PointRows* rows = factor->rows + i;
		s2Body* bodyA = bodies + rows->indexA;
		s2Body* bodyB = bodies + rows->indexB;

		s2Vec2 vrA = s2Add(bodyA->linearVelocity, s2CrossSV(bodyA->angularVelocity, rows->rA));
		s2Vec2 vrB = s2Add(bodyB->linearVelocity, s2CrossSV(bodyB->angularVelocity, rows->rB));
		s2Vec2 rhs = s2Sub(vrA, vrB);

		if (useBias)
		{
			// Position error from current anchors even when the Jacobian anchors are fixed
			s2Vec2 rA = s2RotateVector(bodyA->rot, joint->localAnchorA);
			s2Vec2 rB = s2RotateVector(bodyB->rot, joint->localAnchorB);
			s2Vec2 dcA = bodyA->deltaPosition;
			s2Vec2 dcB = bodyB->deltaPosition;
			s2Vec2 separation = s2Add(s2Add(s2Sub(dcB, dcA), s2Sub(rB, rA)), joint->centerDiff0);
			rhs = s2MulSub(rhs, joint->biasCoefficient, separation);

			// Soft Step's impulse leak
			float K[2][2];
			s2ComputeJointBlock(world, rows, rows, K);
			float compliance = factor->compliance[i];
			s2Vec2 lambda = joint->impulse;
			rhs.x -= compliance * (K[0][0] * lambda.x + K[0][1] * lambda.y);
			rhs.y -= compliance * (K[1][0] * lambda.x + K[1][1] * lambda.y);
		}

		factor->rhs[i] = rhs;
	}

	s2SolveJointFactor(factor);

	for (int i = 0; i < jointCount; ++i)
	{
		s2RevoluteJoint* joint = &world->joints[factor->jointIndices[i]].revoluteJoint;
		const s2PointRows* rows = factor->rows + i;
		s2Body* bodyA = bodies + rows->indexA;
		s2Body* bodyB = bodies + rows->indexB;

		s2Vec2 impulse = factor->rhs[i];
		joint->impulse = s2Add(joint->impulse, impulse);

		bodyA->linearVelocity = s2MulSub(bodyA->linearVelocity, joint->invMassA, impulse);
		bodyA->angularVelocity -= joint->invIA * s2Cross(rows->rA, impulse);
		bodyB->linearVelocity = s2MulAdd(bodyB->linearVelocity, joint->invMassB, impulse);
		bodyB->angularVelocity += joint->invIB * s2Cross(rows->rB, impulse);
	}
}

// Joints outside the direct solve get Soft Step's solve. Direct joints get their motor and limits here.
static void s2SolveJointsPGS(s2World* world, s2StepContext* context, const s2JointFactor* factor, float h, float inv_h,
							 bool useBias)
{
	s2Joint* joints = world->joints;
	int jointCapacity = world->jointPool.capacity;

	for (int i = 0; i < jointCapacity; ++i)
	{
		s2Joint* joint = joints + i;
		if (s2IsFree(&joint->object))
		{
			continue;
		}

		if (factor->jointMap[i] == S2_NULL_INDEX)
		{
			s2SolveJoint_Soft(joint, context, h, inv_h, useBias);
		}
		else
		{
			s2SolveRevoluteAxial_Soft(joint, context, h, inv_h, useBias);
		}
	}
}

void s2Solve_SoftStepD(s2World* world, s2StepContext* context)
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
	float jointHertz = S2_MIN(s2_jointHertz, 0.25f * inv_h);

	s2PrepareContacts_Soft(world, constraints, constraintCount, context, h, contactHertz);

	int directCount = 0;
	for (int i = 0; i < jointCapacity; ++i)
	{
		s2Joint* joint = joints + i;
		if (s2IsFree(&joint->object))
		{
			continue;
		}

		bool warmStart = true;
		s2PrepareJoint_Soft(joint, context, h, jointHertz, warmStart);

		if (s2UseDirect(world, joint))
		{
			directCount += 1;
		}
	}

	s2JointFactor factor;
	s2CreateJointFactor(world, &factor, directCount);

	directCount = 0;
	for (int i = 0; i < jointCapacity; ++i)
	{
		s2Joint* joint = joints + i;
		if (s2IsFree(&joint->object) || s2UseDirect(world, joint) == false)
		{
			continue;
		}

		factor.jointIndices[directCount] = i;
		factor.rows[directCount].indexA = joint->edges[0].bodyIndex;
		factor.rows[directCount].indexB = joint->edges[1].bodyIndex;
		directCount += 1;
	}

	s2UpdateDirectAnchors(world, &factor);
	s2OrderJointFactor(world, &factor);

	for (int substep = 0; substep < substepCount; ++substep)
	{
		s2IntegrateVelocities(world, h);

		if (context->warmStart)
		{
#if S2_SOFT_STEP_D_REFACTOR == 1
			s2UpdateDirectAnchors(world, &factor);
#endif

			for (int i = 0; i < jointCapacity; ++i)
			{
				s2Joint* joint = joints + i;
				if (s2IsFree(&joint->object) || factor.jointMap[i] != S2_NULL_INDEX)
				{
					continue;
				}

				s2WarmStartJoint(joint, context);
			}

			s2WarmStartDirect(world, &factor);
			s2WarmStartContacts_Fixed(world, constraints, constraintCount);
		}

		// Direct joints go last so each pass ends with their point rows exact
		bool useBias = true;
		s2SolveJointsPGS(world, context, &factor, h, inv_h, useBias);
		s2SolveContacts_TGS_Fixed(world, constraints, constraintCount, inv_h, useBias);
		s2FactorDirect(world, &factor, true, 0.0f);
		s2SolveDirect(world, &factor, useBias);

		s2IntegratePositions(world, h);

		if (context->extraIterations > 0)
		{
			useBias = false;
			s2SolveJointsPGS(world, context, &factor, h, inv_h, useBias);
			s2SolveContacts_TGS_Fixed(world, constraints, constraintCount, inv_h, useBias);

#if S2_SOFT_STEP_D_RELAX_COMPLIANCE == 1
			s2FactorDirect(world, &factor, true, 0.0f);
#else
			s2FactorDirect(world, &factor, false, S2_SOFT_STEP_D_RELAX_REGULARIZATION);
#endif
			s2SolveDirect(world, &factor, useBias);
		}
	}

	s2FinalizePositions(world);

	s2StoreContactImpulses(constraints, constraintCount);

	s2DestroyJointFactor(world, &factor);
	s2FreeStackItem(world->stackAllocator, constraints);
}
