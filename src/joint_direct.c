// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "joint_direct.h"

#include "body.h"
#include "core.h"
#include "world.h"

s2DirectJoint s2GetDirectJoint(s2Joint* base)
{
	s2DirectJoint direct = {0};
	s2GetJointRowImpulses(base, direct.impulses);

	switch (base->type)
	{
		case s2_revoluteJoint:
		{
			s2RevoluteJoint* joint = &base->revoluteJoint;
			direct.biasCoefficient = joint->biasCoefficient;
			direct.massCoefficient = joint->massCoefficient;
			direct.impulseCoefficient = joint->impulseCoefficient;
			direct.axialImpulse = joint->motorImpulse + joint->lowerImpulse - joint->upperImpulse;
		}
		break;

		case s2_weldJoint:
		{
			s2WeldJoint* joint = &base->weldJoint;
			direct.biasCoefficient = joint->biasCoefficient;
			direct.massCoefficient = joint->massCoefficient;
			direct.impulseCoefficient = joint->impulseCoefficient;
		}
		break;

		case s2_prismaticJoint:
		{
			s2PrismaticJoint* joint = &base->prismaticJoint;
			direct.biasCoefficient = joint->softness.biasRate;
			direct.massCoefficient = joint->softness.massScale;
			direct.impulseCoefficient = joint->softness.impulseScale;
			direct.axialImpulse = joint->springImpulse + joint->lowerImpulse - joint->upperImpulse;
		}
		break;

		default:
			S2_ASSERT(false);
	}

	return direct;
}

static s2Vec2 s2GetSeparation(const s2Body* bodyA, const s2Body* bodyB, s2Vec2 localAnchorA, s2Vec2 localAnchorB,
							  s2Vec2 centerDiff0)
{
	s2Vec2 rA = s2RotateVector(bodyA->rot, localAnchorA);
	s2Vec2 rB = s2RotateVector(bodyB->rot, localAnchorB);
	return s2Add(s2Add(s2Sub(bodyB->deltaPosition, bodyA->deltaPosition), s2Sub(rB, rA)), centerDiff0);
}

void s2GetDirectError(s2World* world, s2Joint* base, float C[3])
{
	s2Body* bodyA = world->bodies + base->edges[0].bodyIndex;
	s2Body* bodyB = world->bodies + base->edges[1].bodyIndex;

	switch (base->type)
	{
		case s2_revoluteJoint:
		{
			s2RevoluteJoint* joint = &base->revoluteJoint;
			s2Vec2 separation = s2GetSeparation(bodyA, bodyB, joint->localAnchorA, joint->localAnchorB, joint->centerDiff0);
			C[0] = separation.x;
			C[1] = separation.y;
		}
		break;

		case s2_weldJoint:
		{
			s2WeldJoint* joint = &base->weldJoint;
			s2Vec2 separation = s2GetSeparation(bodyA, bodyB, joint->localAnchorA, joint->localAnchorB, joint->centerDiff0);
			C[0] = separation.x;
			C[1] = separation.y;
			C[2] = s2RelativeAngle(bodyB->rot, bodyA->rot) - joint->referenceAngle;
		}
		break;

		case s2_prismaticJoint:
		{
			s2PrismaticJoint* joint = &base->prismaticJoint;
			s2Vec2 d = s2GetSeparation(bodyA, bodyB, joint->localAnchorA, joint->localAnchorB, joint->centerDiff0);
			s2Vec2 perp = s2LeftPerp(s2RotateVector(bodyA->rot, joint->localAxisA));
			C[0] = s2Dot(perp, d);
			C[1] = s2RelativeAngle(bodyB->rot, bodyA->rot) - joint->referenceAngle;
		}
		break;

		default:
			S2_ASSERT(false);
	}
}

static bool s2UseDirect(s2World* world, s2Joint* joint, bool all)
{
	if (joint->type != s2_revoluteJoint && joint->type != s2_weldJoint && joint->type != s2_prismaticJoint)
	{
		return false;
	}

	if (all == false)
	{
		bool enableDirect = joint->type == s2_revoluteJoint ? joint->revoluteJoint.enableDirect
							: joint->type == s2_weldJoint	? joint->weldJoint.enableDirect
															: joint->prismaticJoint.enableDirect;
		if (enableDirect == false)
		{
			return false;
		}
	}

	s2Body* bodyA = world->bodies + joint->edges[0].bodyIndex;
	s2Body* bodyB = world->bodies + joint->edges[1].bodyIndex;
	return bodyA->type == s2_dynamicBody || bodyB->type == s2_dynamicBody;
}

void s2CreateDirectFactor(s2World* world, s2JointFactor* factor, bool all)
{
	s2Joint* joints = world->joints;
	int jointCapacity = world->jointPool.capacity;

	int directCount = 0;
	for (int i = 0; i < jointCapacity; ++i)
	{
		s2Joint* joint = joints + i;
		if (s2IsFree(&joint->object) == false && s2UseDirect(world, joint, all))
		{
			directCount += 1;
		}
	}

	s2CreateJointFactor(world, factor, directCount);

	directCount = 0;
	for (int i = 0; i < jointCapacity; ++i)
	{
		s2Joint* joint = joints + i;
		if (s2IsFree(&joint->object) || s2UseDirect(world, joint, all) == false)
		{
			continue;
		}

		factor->jointIndices[directCount] = i;
		directCount += 1;
	}

	s2UpdateDirectRows(world, factor);
	s2OrderJointFactor(world, factor);
}

void s2UpdateDirectRows(s2World* world, s2JointFactor* factor)
{
	for (int i = 0; i < factor->jointCount; ++i)
	{
		s2GetJointRows(world, world->joints + factor->jointIndices[i], factor->rows + i);
	}
}

void s2WarmStartDirect(s2World* world, s2JointFactor* factor, float scale)
{
	for (int i = 0; i < factor->jointCount; ++i)
	{
		s2Joint* base = world->joints + factor->jointIndices[i];
		s2DirectJoint direct = s2GetDirectJoint(base);
		const s2JointRows* rows = factor->rows + i;
		s2Body* bodyA = world->bodies + rows->indexA;
		s2Body* bodyB = world->bodies + rows->indexB;

		float lambda[3] = {0.0f, 0.0f, 0.0f};
		for (int r = 0; r < rows->rowCount; ++r)
		{
			lambda[r] = scale * *direct.impulses[r];
		}

		s2Vec2 P = s2GetPointImpulse(rows, lambda);
		float angularImpulse = s2GetAngleImpulse(rows, lambda);
		float axialImpulse = scale * direct.axialImpulse;

		if (base->type == s2_prismaticJoint)
		{
			// Limits push along the axis at the anchors of the perpendicular row
			P = s2MulAdd(P, axialImpulse, s2RightPerp(rows->axes[0]));
		}
		else
		{
			angularImpulse += axialImpulse;
		}

		bodyA->linearVelocity = s2MulSub(bodyA->linearVelocity, bodyA->invMass, P);
		bodyA->angularVelocity -= bodyA->invI * (s2Cross(rows->rA, P) + angularImpulse);
		bodyB->linearVelocity = s2MulAdd(bodyB->linearVelocity, bodyB->invMass, P);
		bodyB->angularVelocity += bodyB->invI * (s2Cross(rows->rB, P) + angularImpulse);
	}
}

void s2BuildDirectRhs(s2World* world, s2JointFactor* factor, bool useBias, const s2Softness* softness, float scale)
{
	s2Body* bodies = world->bodies;
	int jointCount = factor->jointCount;

	for (int i = 0; i < jointCount; ++i)
	{
		s2Joint* base = world->joints + factor->jointIndices[i];
		const s2JointRows* rows = factor->rows + i;
		s2Body* bodyA = bodies + rows->indexA;
		s2Body* bodyB = bodies + rows->indexB;

		s2Vec2 vrA = s2Add(bodyA->linearVelocity, s2CrossSV(bodyA->angularVelocity, rows->rA));
		s2Vec2 vrB = s2Add(bodyB->linearVelocity, s2CrossSV(bodyB->angularVelocity, rows->rB));
		s2Vec2 dv = s2Sub(vrA, vrB);
		float dw = bodyA->angularVelocity - bodyB->angularVelocity;

		s2JointVector rhs = {0};
		for (int r = 0; r < rows->rowCount; ++r)
		{
			rhs.row[r] = r < rows->pointCount ? s2Dot(rows->axes[r], dv) : dw;
		}

		if (useBias)
		{
			s2DirectJoint direct = s2GetDirectJoint(base);
			float biasCoefficient = softness != NULL ? softness->biasRate : direct.biasCoefficient;

			float C[3];
			s2GetDirectError(world, base, C);

			// Soft Step's impulse leak
			float K[3][3];
			s2ComputeJointBlock(world, rows, rows, K);
			float compliance = factor->compliance[i];

			for (int r = 0; r < rows->rowCount; ++r)
			{
				float leak = 0.0f;
				for (int c = 0; c < rows->rowCount; ++c)
				{
					leak += K[r][c] * (scale * *direct.impulses[c]);
				}

				rhs.row[r] -= biasCoefficient * C[r];
				rhs.row[r] -= compliance * leak;
			}
		}

		factor->rhs[i] = rhs;
	}
}

void s2ApplyDirectSolution(s2World* world, s2JointFactor* factor, float scale, bool applyToBodies)
{
	s2Body* bodies = world->bodies;

	for (int i = 0; i < factor->jointCount; ++i)
	{
		const s2JointRows* rows = factor->rows + i;
		const float* lambda = factor->rhs[i].row;
		float* impulses[3];
		s2GetJointRowImpulses(world->joints + factor->jointIndices[i], impulses);
		for (int r = 0; r < rows->rowCount; ++r)
		{
			*impulses[r] += lambda[r] / scale;
		}

		if (applyToBodies == false)
		{
			continue;
		}

		s2Body* bodyA = bodies + rows->indexA;
		s2Body* bodyB = bodies + rows->indexB;

		s2Vec2 P = s2GetPointImpulse(rows, lambda);
		float angularImpulse = s2GetAngleImpulse(rows, lambda);

		bodyA->linearVelocity = s2MulSub(bodyA->linearVelocity, bodyA->invMass, P);
		bodyA->angularVelocity -= bodyA->invI * (s2Cross(rows->rA, P) + angularImpulse);
		bodyB->linearVelocity = s2MulAdd(bodyB->linearVelocity, bodyB->invMass, P);
		bodyB->angularVelocity += bodyB->invI * (s2Cross(rows->rB, P) + angularImpulse);
	}
}
