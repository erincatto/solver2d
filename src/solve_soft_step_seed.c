// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "body.h"
#include "contact.h"
#include "core.h"
#include "joint.h"
#include "joint_direct.h"
#include "joint_factor.h"
#include "solvers.h"
#include "stack_allocator.h"
#include "world.h"

#include <stdbool.h>

// Soft Step with an exact joint solve feeding the warm start. The substeps are Soft Step's, Gauss-Seidel joints
// included. Only the stored impulses of the revolute, weld, and prismatic equality rows change: the exact
// solution of the biased joint system the next substep starts from. The velocities of that solve are discarded.

static void s2WarmStartContactsScaled(s2World* world, s2ContactConstraint* constraints, int constraintCount, float scale)
{
	s2Body* bodies = world->bodies;

	for (int i = 0; i < constraintCount; ++i)
	{
		s2ContactConstraint* constraint = constraints + i;
		s2Body* bodyA = bodies + constraint->indexA;
		s2Body* bodyB = bodies + constraint->indexB;

		s2Vec2 normal = constraint->normal;
		s2Vec2 tangent = s2RightPerp(normal);

		for (int j = 0; j < constraint->pointCount; ++j)
		{
			s2ContactConstraintPoint* cp = constraint->points + j;
			s2Vec2 P = s2MulSV(scale, s2Add(s2MulSV(cp->normalImpulse, normal), s2MulSV(cp->tangentImpulse, tangent)));
			bodyA->angularVelocity -= bodyA->invI * s2Cross(cp->rA0, P);
			bodyA->linearVelocity = s2MulAdd(bodyA->linearVelocity, -bodyA->invMass, P);
			bodyB->angularVelocity += bodyB->invI * s2Cross(cp->rB0, P);
			bodyB->linearVelocity = s2MulAdd(bodyB->linearVelocity, bodyB->invMass, P);
		}
	}
}

// Solves the biased joint system from the warm started state over an interval of scale substeps, with the given
// softness or each joint's own. Body velocities are restored afterward.
static void s2SeedJoints(s2World* world, s2StepContext* context, s2JointFactor* factor, s2ContactConstraint* constraints,
						 int constraintCount, float interval, float scale, const s2Softness* softness, s2Vec2* linear,
						 float* angular)
{
	s2Body* bodies = world->bodies;
	int bodyCapacity = world->bodyPool.capacity;
	s2Joint* joints = world->joints;
	int jointCapacity = world->jointPool.capacity;

	for (int i = 0; i < bodyCapacity; ++i)
	{
		linear[i] = bodies[i].linearVelocity;
		angular[i] = bodies[i].angularVelocity;
	}

	// The substep loop already integrated a single substep
	if (scale > 1.0f)
	{
		s2IntegrateVelocities(world, interval);
	}

	for (int i = 0; i < jointCapacity; ++i)
	{
		s2Joint* joint = joints + i;
		if (s2IsFree(&joint->object) || factor->jointMap[i] != S2_NULL_INDEX)
		{
			continue;
		}

		s2WarmStartJoint(joint, context);
	}

	s2UpdateDirectRows(world, factor);
	s2WarmStartDirect(world, factor, scale);
	s2WarmStartContactsScaled(world, constraints, constraintCount, scale);

	for (int i = 0; i < factor->jointCount; ++i)
	{
		s2DirectJoint direct = s2GetDirectJoint(joints + factor->jointIndices[i]);
		factor->compliance[i] = softness != NULL ? softness->impulseScale / softness->massScale
												 : direct.impulseCoefficient / direct.massCoefficient;
	}

	s2FactorJoints(world, factor, 0.0f);
	s2BuildDirectRhs(world, factor, true, softness, scale);
	s2SolveJointFactor(factor);
	s2ApplyDirectSolution(world, factor, scale, false);

	for (int i = 0; i < bodyCapacity; ++i)
	{
		bodies[i].linearVelocity = linear[i];
		bodies[i].angularVelocity = angular[i];
	}
}

void s2Solve_SoftStepSeed(s2World* world, s2StepContext* context)
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

	// Defaults to Soft Step's cap. The seed may allow stiffer joints.
	float jointHertz = world->jointHertz;
	if (world->jointHertzCap > 0.0f)
	{
		jointHertz = S2_MIN(jointHertz, world->jointHertzCap * inv_h);
	}

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
	}

	s2JointFactor factor;
	s2CreateDirectFactor(world, &factor, true);

	int bodyCapacity = world->bodyPool.capacity;
	s2Vec2* linear = s2AllocateStackItem(world->stackAllocator, bodyCapacity * sizeof(s2Vec2), "seed linear");
	float* angular = s2AllocateStackItem(world->stackAllocator, bodyCapacity * sizeof(float), "seed angular");

	s2SeedMode mode = world->jointSeedMode;

	if (context->warmStart && mode == s2_seedWholeStep)
	{
		// Softness over the whole step, so the impulses are a step's worth until divided among the substeps
		float dt = h * substepCount;
		s2Softness softness = s2MakeSoft(jointHertz, 10.0f, dt);
		s2SeedJoints(world, context, &factor, constraints, constraintCount, dt, (float)substepCount, &softness, linear,
					 angular);
	}

	for (int substep = 0; substep < substepCount; ++substep)
	{
		// Integrate gravity and forces
		s2IntegrateVelocities(world, h);

		// Apply warm starting
		if (context->warmStart)
		{
			if (mode == s2_seedEverySubstep || (mode == s2_seedFirstSubstep && substep == 0))
			{
				s2SeedJoints(world, context, &factor, constraints, constraintCount, h, 1.0f, NULL, linear, angular);
			}

			for (int i = 0; i < jointCapacity; ++i)
			{
				s2Joint* joint = joints + i;
				if (s2IsFree(&joint->object))
				{
					continue;
				}

				s2WarmStartJoint(joint, context);
			}

			s2WarmStartContacts_Fixed(world, constraints, constraintCount);
		}

		// Solve velocities using position bias
		bool useBias = true;
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

		// Integrate positions using biased velocities
		s2IntegratePositions(world, h);

		// Relax biased velocities and impulses
		if (context->extraIterations > 0)
		{
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
	}

	s2FinalizePositions(world);

	s2StoreContactImpulses(constraints, constraintCount);

	s2FreeStackItem(world->stackAllocator, angular);
	s2FreeStackItem(world->stackAllocator, linear);
	s2DestroyJointFactor(world, &factor);
	s2FreeStackItem(world->stackAllocator, constraints);
}
