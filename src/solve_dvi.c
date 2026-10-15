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

#include <float.h>
#include <math.h>
#include <stdbool.h>

// DVI solver from Kamino in Newton physics. One solve per step followed by semi-implicit Euler.
// Joint point rows are solved directly with a sparse Cholesky factorization of their Delassus block.
// Contacts, joint limits, and motors use projected Gauss-Seidel. The two alternate, so each block
// sees the latest impulses of the other.
// Tsounis, Grandia, Bacher: On Solving the Dynamics of Constrained Rigid Multi-Body Systems with Kinematic Loops,
// arXiv:2504.19771

// Kamino defaults. The Baumgarte factors are the fraction of position error removed per step.
#define S2_DVI_JOINT_BAUMGARTE 0.01f
#define S2_DVI_LIMIT_BAUMGARTE 0.01f
#define S2_DVI_CONTACT_BAUMGARTE 0.01f
#define S2_DVI_CONTACT_DEAD_ZONE 1.0e-6f
#define S2_DVI_TANGENT_WARM_START_SCALE 0.97f

// Kamino drops contact points that are not yet touching
#define S2_DVI_CULL_SPECULATIVE 1

// Added to the unit diagonal of the scaled joint block. Joint loops make the block singular.
#define S2_DVI_JOINT_REGULARIZATION 7.0e-7f

typedef struct s2DVIContactPoint
{
	s2Vec2 rA, rB;
	float normalImpulse;
	float tangentImpulse;
	float normalMass;
	float tangentMass;

	// Position recovery velocity, kept out of the friction load
	float stabilization;

	// Newton impact target from the normal velocity at the start of the step
	float restitution;

	bool active;
} s2DVIContactPoint;

typedef struct s2DVIContact
{
	s2Contact* contact;
	int indexA;
	int indexB;
	s2Vec2 normal;
	float friction;
	int pointCount;
	s2DVIContactPoint points[2];
} s2DVIContact;

typedef struct s2DVIJoint
{
	s2Joint* joint;
	int indexA;
	int indexB;
	s2Vec2 rA, rB;
	s2Vec2 bias;
	float angularBias;
	float axialMass;
	float limitBias;

	// Two point rows, plus the angle row for welds
	int rowCount;

	// 1 at the lower limit, -1 at the upper limit, 0 within limits
	float limitSide;
} s2DVIJoint;

static void s2PrepareContacts_DVI(s2World* world, s2DVIContact* constraints, int constraintCount, float inv_h,
								  bool warmStart)
{
	s2Body* bodies = world->bodies;

	for (int i = 0; i < constraintCount; ++i)
	{
		s2DVIContact* constraint = constraints + i;

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

		s2Vec2 vA = bodyA->linearVelocity;
		float wA = bodyA->angularVelocity;
		s2Vec2 vB = bodyB->linearVelocity;
		float wB = bodyB->angularVelocity;

		s2Vec2 normal = constraint->normal;
		s2Vec2 tangent = s2RightPerp(normal);
		float restitution = contact->restitution;

		for (int j = 0; j < pointCount; ++j)
		{
			const s2ManifoldPoint* mp = manifold->points + j;
			s2DVIContactPoint* cp = constraint->points + j;

			s2Vec2 rA = s2RotateVector(bodyA->rot, s2Sub(mp->localAnchorA, bodyA->localCenter));
			s2Vec2 rB = s2RotateVector(bodyB->rot, s2Sub(mp->localAnchorB, bodyB->localCenter));
			cp->rA = rA;
			cp->rB = rB;

			float rnA = s2Cross(rA, normal);
			float rnB = s2Cross(rB, normal);
			float kNormal = mA + mB + iA * rnA * rnA + iB * rnB * rnB;
			cp->normalMass = kNormal > 0.0f ? 1.0f / kNormal : 0.0f;

			float rtA = s2Cross(rA, tangent);
			float rtB = s2Cross(rB, tangent);
			float kTangent = mA + mB + iA * rtA * rtA + iB * rtB * rtB;
			cp->tangentMass = kTangent > 0.0f ? 1.0f / kTangent : 0.0f;

			float separation = mp->separation;

#if S2_DVI_CULL_SPECULATIVE == 1
			cp->active = separation <= 0.0f;
#else
			cp->active = true;
#endif

			if (cp->active == false)
			{
				cp->normalImpulse = 0.0f;
				cp->tangentImpulse = 0.0f;
				cp->stabilization = 0.0f;
				cp->restitution = 0.0f;
				continue;
			}

			// The dead zone filters noise on touching contacts and stays continuous in separation
			float depth = S2_MAX(0.0f, S2_ABS(separation) - S2_DVI_CONTACT_DEAD_ZONE);
			float xi = inv_h * (separation < 0.0f ? -depth : depth);

			// One sided, approach a speculative gap fully but push out of overlap slowly
			cp->stabilization = S2_DVI_CONTACT_BAUMGARTE * S2_MIN(0.0f, xi) + S2_MAX(0.0f, xi);

			// Elastic impacts are destabilized by the extra push
			if (restitution >= 1.0f)
			{
				cp->stabilization = 0.0f;
			}

			s2Vec2 vrA = s2Add(vA, s2CrossSV(wA, rA));
			s2Vec2 vrB = s2Add(vB, s2CrossSV(wB, rB));
			cp->restitution = restitution * s2Dot(s2Sub(vrB, vrA), normal);

			if (warmStart)
			{
				cp->normalImpulse = mp->normalImpulse;
				cp->tangentImpulse = S2_DVI_TANGENT_WARM_START_SCALE * mp->tangentImpulse;
			}
			else
			{
				cp->normalImpulse = 0.0f;
				cp->tangentImpulse = 0.0f;
			}
		}
	}
}

static void s2PrepareJoint_DVI(s2World* world, s2DVIJoint* dj, float inv_h, bool warmStart)
{
	s2Joint* base = dj->joint;

	s2Body* bodyA = world->bodies + dj->indexA;
	s2Body* bodyB = world->bodies + dj->indexB;

	s2Rot qA = bodyA->rot;
	s2Rot qB = bodyB->rot;

	dj->rA = s2RotateVector(qA, s2Sub(base->localOriginAnchorA, bodyA->localCenter));
	dj->rB = s2RotateVector(qB, s2Sub(base->localOriginAnchorB, bodyB->localCenter));

	s2Vec2 separation = s2Sub(s2Add(bodyB->position, dj->rB), s2Add(bodyA->position, dj->rA));
	dj->bias = s2MulSV(S2_DVI_JOINT_BAUMGARTE * inv_h, separation);

	float k = bodyA->invI + bodyB->invI;
	dj->axialMass = k > 0.0f ? 1.0f / k : 0.0f;

	dj->limitSide = 0.0f;
	dj->limitBias = 0.0f;

	if (base->type == s2_weldJoint)
	{
		s2WeldJoint* joint = &base->weldJoint;
		float angle = s2RelativeAngle(qB, qA) - joint->referenceAngle;
		dj->angularBias = S2_DVI_JOINT_BAUMGARTE * inv_h * angle;
		dj->rowCount = 3;

		if (warmStart == false)
		{
			joint->linearImpulse = s2Vec2_zero;
			joint->angularImpulse = 0.0f;
		}

		return;
	}

	S2_ASSERT(base->type == s2_revoluteJoint);
	s2RevoluteJoint* joint = &base->revoluteJoint;
	dj->angularBias = 0.0f;
	dj->rowCount = 2;

	// Kamino only creates a limit row once the limit is violated
	if (joint->enableLimit)
	{
		float angle = s2RelativeAngle(qB, qA) - joint->referenceAngle;
		if (angle < joint->lowerAngle)
		{
			dj->limitSide = 1.0f;
			dj->limitBias = S2_DVI_LIMIT_BAUMGARTE * inv_h * (angle - joint->lowerAngle);
		}
		else if (angle > joint->upperAngle)
		{
			dj->limitSide = -1.0f;
			dj->limitBias = S2_DVI_LIMIT_BAUMGARTE * inv_h * (joint->upperAngle - angle);
		}
	}

	if (warmStart == false)
	{
		joint->impulse = s2Vec2_zero;
		joint->motorImpulse = 0.0f;
		joint->lowerImpulse = 0.0f;
		joint->upperImpulse = 0.0f;
	}

	if (dj->limitSide != 1.0f)
	{
		joint->lowerImpulse = 0.0f;
	}

	if (dj->limitSide != -1.0f)
	{
		joint->upperImpulse = 0.0f;
	}

	if (joint->enableMotor == false)
	{
		joint->motorImpulse = 0.0f;
	}
}

static void s2ApplyImpulse(s2Body* bodyA, s2Body* bodyB, s2Vec2 rA, s2Vec2 rB, s2Vec2 P)
{
	bodyA->linearVelocity = s2MulSub(bodyA->linearVelocity, bodyA->invMass, P);
	bodyA->angularVelocity -= bodyA->invI * s2Cross(rA, P);
	bodyB->linearVelocity = s2MulAdd(bodyB->linearVelocity, bodyB->invMass, P);
	bodyB->angularVelocity += bodyB->invI * s2Cross(rB, P);
}

static void s2ApplyAngularImpulse(s2Body* bodyA, s2Body* bodyB, float impulse)
{
	bodyA->angularVelocity -= bodyA->invI * impulse;
	bodyB->angularVelocity += bodyB->invI * impulse;
}

static void s2WarmStart_DVI(s2World* world, s2DVIContact* contacts, int contactCount, s2DVIJoint* joints, int jointCount)
{
	s2Body* bodies = world->bodies;

	for (int i = 0; i < contactCount; ++i)
	{
		s2DVIContact* constraint = contacts + i;
		s2Body* bodyA = bodies + constraint->indexA;
		s2Body* bodyB = bodies + constraint->indexB;
		s2Vec2 normal = constraint->normal;
		s2Vec2 tangent = s2RightPerp(normal);

		for (int j = 0; j < constraint->pointCount; ++j)
		{
			s2DVIContactPoint* cp = constraint->points + j;
			s2Vec2 P = s2Add(s2MulSV(cp->normalImpulse, normal), s2MulSV(cp->tangentImpulse, tangent));
			s2ApplyImpulse(bodyA, bodyB, cp->rA, cp->rB, P);
		}
	}

	for (int i = 0; i < jointCount; ++i)
	{
		s2DVIJoint* dj = joints + i;
		s2Body* bodyA = bodies + dj->indexA;
		s2Body* bodyB = bodies + dj->indexB;

		if (dj->joint->type == s2_weldJoint)
		{
			s2WeldJoint* joint = &dj->joint->weldJoint;
			s2ApplyImpulse(bodyA, bodyB, dj->rA, dj->rB, joint->linearImpulse);
			s2ApplyAngularImpulse(bodyA, bodyB, joint->angularImpulse);
			continue;
		}

		s2RevoluteJoint* joint = &dj->joint->revoluteJoint;
		s2ApplyImpulse(bodyA, bodyB, dj->rA, dj->rB, joint->impulse);
		s2ApplyAngularImpulse(bodyA, bodyB, joint->motorImpulse + joint->lowerImpulse - joint->upperImpulse);
	}
}

// Solve the joint block given the current inequality impulses. Matrix free this is the impulse that zeroes the
// current biased joint velocity, which matches Kamino's D_bb * lambda_b = -(v_f + D_bu * lambda_u).
static void s2SolveJointBlock(s2World* world, s2DVIJoint* joints, s2JointFactor* factor)
{
	s2Body* bodies = world->bodies;
	int jointCount = factor->jointCount;

	for (int i = 0; i < jointCount; ++i)
	{
		const s2DVIJoint* dj = joints + i;
		s2Body* bodyA = bodies + dj->indexA;
		s2Body* bodyB = bodies + dj->indexB;

		s2Vec2 vrA = s2Add(bodyA->linearVelocity, s2CrossSV(bodyA->angularVelocity, dj->rA));
		s2Vec2 vrB = s2Add(bodyB->linearVelocity, s2CrossSV(bodyB->angularVelocity, dj->rB));
		s2Vec2 Cdot = s2Add(s2Sub(vrB, vrA), dj->bias);

		factor->rhs[i].linear = s2Neg(Cdot);
		factor->rhs[i].angular = -(bodyB->angularVelocity - bodyA->angularVelocity + dj->angularBias);
	}

	s2SolveJointFactor(factor);

	for (int p = 0; p < jointCount; ++p)
	{
		s2DVIJoint* dj = joints + factor->order[p];
		s2JointVector impulse = factor->rhs[factor->order[p]];
		s2Body* bodyA = bodies + dj->indexA;
		s2Body* bodyB = bodies + dj->indexB;

		s2ApplyImpulse(bodyA, bodyB, dj->rA, dj->rB, impulse.linear);

		if (dj->joint->type == s2_weldJoint)
		{
			s2WeldJoint* joint = &dj->joint->weldJoint;
			joint->linearImpulse = s2Add(joint->linearImpulse, impulse.linear);
			joint->angularImpulse += impulse.angular;
			s2ApplyAngularImpulse(bodyA, bodyB, impulse.angular);
		}
		else
		{
			s2RevoluteJoint* joint = &dj->joint->revoluteJoint;
			joint->impulse = s2Add(joint->impulse, impulse.linear);
		}
	}
}

// Kamino's bounded rows are the motor and its limit rows only exist while violated
static void s2SolveJointInequalities(s2World* world, s2DVIJoint* joints, int jointCount, float h)
{
	s2Body* bodies = world->bodies;

	for (int i = 0; i < jointCount; ++i)
	{
		s2DVIJoint* dj = joints + i;
		if (dj->joint->type != s2_revoluteJoint)
		{
			continue;
		}

		s2RevoluteJoint* joint = &dj->joint->revoluteJoint;
		s2Body* bodyA = bodies + dj->indexA;
		s2Body* bodyB = bodies + dj->indexB;

		if (joint->enableMotor)
		{
			float Cdot = bodyB->angularVelocity - bodyA->angularVelocity - joint->motorSpeed;
			float maxImpulse = h * joint->maxMotorTorque;
			float oldImpulse = joint->motorImpulse;
			joint->motorImpulse = S2_CLAMP(oldImpulse - dj->axialMass * Cdot, -maxImpulse, maxImpulse);
			s2ApplyAngularImpulse(bodyA, bodyB, joint->motorImpulse - oldImpulse);
		}

		if (dj->limitSide == 1.0f)
		{
			float Cdot = bodyB->angularVelocity - bodyA->angularVelocity + dj->limitBias;
			float oldImpulse = joint->lowerImpulse;
			joint->lowerImpulse = S2_MAX(0.0f, oldImpulse - dj->axialMass * Cdot);
			s2ApplyAngularImpulse(bodyA, bodyB, joint->lowerImpulse - oldImpulse);
		}
		else if (dj->limitSide == -1.0f)
		{
			float Cdot = bodyA->angularVelocity - bodyB->angularVelocity + dj->limitBias;
			float oldImpulse = joint->upperImpulse;
			joint->upperImpulse = S2_MAX(0.0f, oldImpulse - dj->axialMass * Cdot);
			s2ApplyAngularImpulse(bodyA, bodyB, oldImpulse - joint->upperImpulse);
		}
	}
}

static void s2SolveContactNormals(s2World* world, s2DVIContact* contacts, int contactCount)
{
	s2Body* bodies = world->bodies;

	for (int i = 0; i < contactCount; ++i)
	{
		s2DVIContact* constraint = contacts + i;
		s2Body* bodyA = bodies + constraint->indexA;
		s2Body* bodyB = bodies + constraint->indexB;
		s2Vec2 normal = constraint->normal;

		for (int j = 0; j < constraint->pointCount; ++j)
		{
			s2DVIContactPoint* cp = constraint->points + j;
			if (cp->active == false)
			{
				continue;
			}

			s2Vec2 vrA = s2Add(bodyA->linearVelocity, s2CrossSV(bodyA->angularVelocity, cp->rA));
			s2Vec2 vrB = s2Add(bodyB->linearVelocity, s2CrossSV(bodyB->angularVelocity, cp->rB));
			float vn = s2Dot(s2Sub(vrB, vrA), normal) + cp->stabilization + cp->restitution;

			float oldImpulse = cp->normalImpulse;
			cp->normalImpulse = S2_MAX(0.0f, oldImpulse - cp->normalMass * vn);
			s2ApplyImpulse(bodyA, bodyB, cp->rA, cp->rB, s2MulSV(cp->normalImpulse - oldImpulse, normal));
		}
	}
}

static void s2SolveContactTangent(s2World* world, s2DVIContact* constraint)
{
	s2Body* bodyA = world->bodies + constraint->indexA;
	s2Body* bodyB = world->bodies + constraint->indexB;
	s2Vec2 tangent = s2RightPerp(constraint->normal);
	float friction = constraint->friction;

	for (int j = 0; j < constraint->pointCount; ++j)
	{
		s2DVIContactPoint* cp = constraint->points + j;
		if (cp->active == false)
		{
			continue;
		}

		// Friction is bounded by the support load, not by the extra push out of overlap
		float load = S2_CLAMP(cp->normalImpulse + cp->normalMass * cp->stabilization, 0.0f, cp->normalImpulse);
		float maxFriction = friction * load;

		s2Vec2 vrA = s2Add(bodyA->linearVelocity, s2CrossSV(bodyA->angularVelocity, cp->rA));
		s2Vec2 vrB = s2Add(bodyB->linearVelocity, s2CrossSV(bodyB->angularVelocity, cp->rB));
		float vt = s2Dot(s2Sub(vrB, vrA), tangent);

		float oldImpulse = cp->tangentImpulse;
		cp->tangentImpulse = S2_CLAMP(oldImpulse - cp->tangentMass * vt, -maxFriction, maxFriction);
		s2ApplyImpulse(bodyA, bodyB, cp->rA, cp->rB, s2MulSV(cp->tangentImpulse - oldImpulse, tangent));
	}
}

// One projected Gauss-Seidel sweep over the inequalities. All normal rows go before any friction row.
static void s2SolveInequalities(s2World* world, s2StepContext* context, s2DVIJoint* joints, int jointCount,
								s2DVIContact* contacts, int contactCount, bool solveFriction, bool reverseFriction)
{
	s2Joint* worldJoints = world->joints;
	int jointCapacity = world->jointPool.capacity;
	for (int i = 0; i < jointCapacity; ++i)
	{
		s2Joint* joint = worldJoints + i;
		if (s2IsFree(&joint->object) == false && joint->type == s2_mouseJoint)
		{
			s2SolveJoint(joint, context, context->dt);
		}
	}

	s2SolveJointInequalities(world, joints, jointCount, context->dt);
	s2SolveContactNormals(world, contacts, contactCount);

	if (solveFriction == false)
	{
		return;
	}

	// Alternating the order reduces load bias in redundant sticking patches
	if (reverseFriction)
	{
		for (int i = contactCount - 1; i >= 0; --i)
		{
			s2SolveContactTangent(world, contacts + i);
		}
	}
	else
	{
		for (int i = 0; i < contactCount; ++i)
		{
			s2SolveContactTangent(world, contacts + i);
		}
	}
}

void s2Solve_DVI(s2World* world, s2StepContext* context)
{
	s2Body* bodies = world->bodies;
	s2Contact* worldContacts = world->contacts;
	int contactCapacity = world->contactPool.capacity;
	s2Joint* worldJoints = world->joints;
	int jointCapacity = world->jointPool.capacity;
	s2StackAllocator* stack = world->stackAllocator;

	float h = context->dt;
	float inv_h = context->inv_dt;
	bool warmStart = context->warmStart;
	int alternations = S2_MAX(1, context->iterations);
	int sweeps = S2_MAX(1, context->extraIterations);

	s2DVIContact* contacts = s2AllocateStackItem(stack, contactCapacity * sizeof(s2DVIContact), "dvi contacts");
	int contactCount = 0;
	for (int i = 0; i < contactCapacity; ++i)
	{
		s2Contact* contact = worldContacts + i;
		if (s2IsFree(&contact->object) || contact->manifold.pointCount == 0)
		{
			continue;
		}

		contacts[contactCount].contact = contact;
		contact->manifold.constraintIndex = contactCount;
		contactCount += 1;
	}

	s2DVIJoint* joints = s2AllocateStackItem(stack, jointCapacity * sizeof(s2DVIJoint), "dvi joints");
	int jointCount = 0;
	bool hasInequalities = contactCount > 0;
	for (int i = 0; i < jointCapacity; ++i)
	{
		s2Joint* joint = worldJoints + i;
		if (s2IsFree(&joint->object))
		{
			continue;
		}

		if (joint->type == s2_mouseJoint)
		{
			s2PrepareJoint(joint, context, warmStart);
			hasInequalities = true;
			continue;
		}

		int indexA = joint->edges[0].bodyIndex;
		int indexB = joint->edges[1].bodyIndex;
		if (bodies[indexA].type != s2_dynamicBody && bodies[indexB].type != s2_dynamicBody)
		{
			continue;
		}

		s2DVIJoint* dj = joints + jointCount;
		dj->joint = joint;
		dj->indexA = indexA;
		dj->indexB = indexB;
		s2PrepareJoint_DVI(world, dj, inv_h, warmStart);
		if (joint->type == s2_revoluteJoint)
		{
			hasInequalities = hasInequalities || joint->revoluteJoint.enableMotor || dj->limitSide != 0.0f;
		}

		jointCount += 1;
	}

	// Restitution needs the normal velocity before gravity is applied
	s2PrepareContacts_DVI(world, contacts, contactCount, inv_h, warmStart);

	s2IntegrateVelocities(world, h);

	s2WarmStart_DVI(world, contacts, contactCount, joints, jointCount);
	for (int i = 0; i < jointCapacity; ++i)
	{
		s2Joint* joint = worldJoints + i;
		if (warmStart && s2IsFree(&joint->object) == false && joint->type == s2_mouseJoint)
		{
			s2WarmStartJoint(joint, context);
		}
	}

	s2JointFactor factor;
	s2CreateJointFactor(world, &factor, jointCount);
	for (int i = 0; i < jointCount; ++i)
	{
		const s2DVIJoint* dj = joints + i;
		factor.jointIndices[i] = dj->joint->object.index;
		factor.rows[i] = (s2JointRows){dj->indexA, dj->indexB, dj->rA, dj->rB, dj->rowCount};
	}
	s2OrderJointFactor(world, &factor);

	if (jointCount > 0)
	{
		s2FactorJoints(world, &factor, S2_DVI_JOINT_REGULARIZATION);
		s2SolveJointBlock(world, joints, &factor);

		if (hasInequalities)
		{
			for (int iteration = 0; iteration < alternations; ++iteration)
			{
				for (int sweep = 0; sweep < sweeps; ++sweep)
				{
					s2SolveInequalities(world, context, joints, jointCount, contacts, contactCount, true, sweep % 2 == 1);
				}

				s2SolveJointBlock(world, joints, &factor);
			}
		}
	}
	else
	{
		// Without joint rows Kamino runs one long schedule. Friction joins late so the normal rows first
		// establish the support load.
		int sweepCount = (sweeps + 1) * alternations;
		int frictionSweepCount = sweeps * alternations / 2;
		for (int sweep = 0; sweep < sweepCount; ++sweep)
		{
			bool solveFriction = sweep >= sweepCount - frictionSweepCount;
			s2SolveInequalities(world, context, joints, jointCount, contacts, contactCount, solveFriction, sweep % 2 == 1);
		}
	}

	s2IntegratePositions(world, h);
	s2FinalizePositions(world);

	for (int i = 0; i < contactCount; ++i)
	{
		s2DVIContact* constraint = contacts + i;
		s2Manifold* manifold = &constraint->contact->manifold;
		for (int j = 0; j < constraint->pointCount; ++j)
		{
			manifold->points[j].normalImpulse = constraint->points[j].normalImpulse;
			manifold->points[j].tangentImpulse = constraint->points[j].tangentImpulse;
		}
	}

	s2DestroyJointFactor(world, &factor);
	s2FreeStackItem(stack, joints);
	s2FreeStackItem(stack, contacts);
}
