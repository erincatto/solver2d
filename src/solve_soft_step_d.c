// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "body.h"
#include "contact.h"
#include "core.h"
#include "joint.h"
#include "joint_cr.h"
#include "joint_direct.h"
#include "joint_factor.h"
#include "solvers.h"
#include "stack_allocator.h"
#include "world.h"

#include <stdbool.h>

// Soft Step with the equality rows of revolute, weld, and prismatic joints solved together by a direct
// factorization. An isolated revolute joint matches Soft Step. Joints sharing a body are solved exactly instead of
// by one Gauss-Seidel pass. See joint_direct.h for the coupled system.

// Jacobian anchors follow the bodies each substep. With 0 they stay at the start of the step, as Soft Step
// does for contacts. Both refactor before every direct solve for now.
#define S2_SOFT_STEP_D_REFACTOR 1

// Solve every revolute, weld, and prismatic joint directly, ignoring the per joint opt in
#define S2_SOFT_STEP_D_ALL_DIRECT 1

// Relax with the soft compliance kept on the left. Reuses the soft matrix but is inexact on long chains.
#define S2_SOFT_STEP_D_RELAX_COMPLIANCE 0

// Added to the unit diagonal of the scaled relax block, which redundant joints make singular.
// Well above float roundoff so impulses that move nothing stay bounded.
#define S2_SOFT_STEP_D_RELAX_REGULARIZATION 1.0e-5f

static float s2GetJointHertz(s2World* world, float inv_h)
{
	if (world->jointHertzCap > 0.0f)
	{
		return S2_MIN(world->jointHertz, world->jointHertzCap * inv_h);
	}

	return world->jointHertz;
}

static void s2FactorDirect(s2World* world, s2JointFactor* factor, bool soft, float regularization)
{
#if S2_SOFT_STEP_D_REFACTOR == 1
	s2UpdateDirectRows(world, factor);
#endif

	for (int i = 0; i < factor->jointCount; ++i)
	{
		s2DirectJoint direct = s2GetDirectJoint(world->joints + factor->jointIndices[i]);
		factor->compliance[i] = soft ? direct.impulseCoefficient / direct.massCoefficient : 0.0f;
	}

	if (world->jointCRIterations > 0)
	{
		factor->regularization = regularization;
		return;
	}

	s2FactorJoints(world, factor, regularization);
}

// Equality rows of all direct joints in one solve against the current factor
static void s2SolveDirect(s2World* world, s2JointFactor* factor, bool useBias)
{
	int jointCount = factor->jointCount;

	s2BuildDirectRhs(world, factor, useBias, NULL, 1.0f);

	if (world->jointCRIterations > 0)
	{
		s2SolveJointCR(world, factor, world->jointCRIterations);
	}
	else
	{
		s2JointVector* measuredRhs = NULL;
		if (world->measureJointCR)
		{
			measuredRhs = s2AllocateStackItem(world->stackAllocator, S2_MAX(1, jointCount) * sizeof(s2JointVector), "cr rhs");
			for (int i = 0; i < jointCount; ++i)
			{
				measuredRhs[i] = factor->rhs[i];
			}
		}

		s2SolveJointFactor(factor);

		if (measuredRhs != NULL)
		{
			s2MeasureJointCR(world, factor, measuredRhs, useBias ? &world->softCRStats : &world->relaxCRStats);
			s2FreeStackItem(world->stackAllocator, measuredRhs);
		}
	}

	s2ApplyDirectSolution(world, factor, 1.0f, true);
}

// Joints outside the direct solve get Soft Step's solve. Direct revolute and prismatic joints get their motor,
// spring, and limits here.
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
		else if (joint->type == s2_revoluteJoint)
		{
			s2SolveRevoluteAxial_Soft(joint, context, h, inv_h, useBias);
		}
		else if (joint->type == s2_prismaticJoint)
		{
			s2SolvePrismaticAxial_Soft(joint, context, inv_h, useBias);
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

	world->softCRStats = (s2JointCRStats){0};
	world->relaxCRStats = (s2JointCRStats){0};

	float contactHertz = S2_MIN(s2_contactHertz, 0.25f * inv_h);
	float jointHertz = s2GetJointHertz(world, inv_h);

	// Prismatic limits stay Gauss-Seidel, so they keep Soft Step's cap
	float limitHertz = S2_MIN(world->jointHertz, 0.25f * inv_h);

	s2PrepareContacts_Soft(world, constraints, constraintCount, context, h, contactHertz);

	for (int i = 0; i < jointCapacity; ++i)
	{
		s2Joint* joint = joints + i;
		if (s2IsFree(&joint->object))
		{
			continue;
		}

		bool warmStart = true;
		s2PrepareJoint_Soft(joint, context, h, jointHertz, warmStart);

		if (joint->type == s2_prismaticJoint)
		{
			s2PreparePrismaticLimit_Soft(joint, h, limitHertz);
		}
	}

	s2JointFactor factor;
	s2CreateDirectFactor(world, &factor, S2_SOFT_STEP_D_ALL_DIRECT == 1);

	for (int substep = 0; substep < substepCount; ++substep)
	{
		s2IntegrateVelocities(world, h);

		if (context->warmStart)
		{
#if S2_SOFT_STEP_D_REFACTOR == 1
			s2UpdateDirectRows(world, &factor);
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

			s2WarmStartDirect(world, &factor, 1.0f);
			s2WarmStartContacts_Fixed(world, constraints, constraintCount);
		}

		// Direct joints go last so each pass ends with their equality rows exact
		bool useBias = true;
		s2SolveJointsPGS(world, context, &factor, h, inv_h, useBias);
		s2SolveContacts_TGS_Fixed(world, constraints, constraintCount, inv_h, useBias);
		s2FactorDirect(world, &factor, true, 0.0f);
		s2SolveDirect(world, &factor, useBias);

		s2IntegratePositions(world, h);

		if (context->extraIterations > 0 && world->enableDirectRelax == false)
		{
			// Soft Step's relax, every joint by Gauss-Seidel. One factorization per substep instead of two.
			useBias = false;
			for (int i = 0; i < jointCapacity; ++i)
			{
				s2Joint* joint = joints + i;
				if (s2IsFree(&joint->object))
				{
					continue;
				}

				s2SolveJoint_Soft(joint, context, h, inv_h, useBias);
			}

			s2SolveContacts_TGS_Fixed(world, constraints, constraintCount, inv_h, useBias);
		}
		else if (context->extraIterations > 0)
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
