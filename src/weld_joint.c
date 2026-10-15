// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "body.h"
#include "core.h"
#include "joint.h"
#include "solvers.h"
#include "world.h"

// Weld joint without a spring. The angle row is solved before the point rows, as Box2D does.

// Point constraint
// C = p2 - p1
// Cdot = v2 + cross(w2, r2) - v1 - cross(w1, r1)
// J = [-I -r1_skew I r2_skew ]

// Angle constraint
// C = angle2 - angle1 - referenceAngle
// Cdot = w2 - w1
// J = [0 0 -1 0 0 1]
// K = invI1 + invI2

void s2PrepareWeld(s2Joint* base, s2StepContext* context, bool warmStart)
{
	S2_ASSERT(base->type == s2_weldJoint);

	s2Body* bodyA = context->bodies + base->edges[0].bodyIndex;
	s2Body* bodyB = context->bodies + base->edges[1].bodyIndex;

	s2WeldJoint* joint = &base->weldJoint;
	joint->localAnchorA = s2Sub(base->localOriginAnchorA, bodyA->localCenter);
	joint->localAnchorB = s2Sub(base->localOriginAnchorB, bodyB->localCenter);
	joint->centerDiff0 = s2Sub(bodyB->position, bodyA->position);
	joint->invMassA = bodyA->invMass;
	joint->invMassB = bodyB->invMass;
	joint->invIA = bodyA->invI;
	joint->invIB = bodyB->invI;

	float k = joint->invIA + joint->invIB;
	joint->axialMass = k > 0.0f ? 1.0f / k : 0.0f;

	joint->biasCoefficient = 0.0f;
	joint->massCoefficient = 1.0f;
	joint->impulseCoefficient = 0.0f;

	if (warmStart == false)
	{
		joint->linearImpulse = s2Vec2_zero;
		joint->angularImpulse = 0.0f;
	}
}

// Same softness as the revolute joint
void s2PrepareWeld_Soft(s2Joint* base, s2StepContext* context, float h, float hertz, bool warmStart)
{
	s2PrepareWeld(base, context, warmStart);

	s2WeldJoint* joint = &base->weldJoint;
	const float zeta = 10.0f;
	float omega = 2.0f * s2_pi * hertz;
	joint->biasCoefficient = omega / (2.0f * zeta + h * omega);
	float c = h * omega * (2.0f * zeta + h * omega);
	joint->impulseCoefficient = 1.0f / (1.0f + c);
	joint->massCoefficient = c * joint->impulseCoefficient;
}

void s2PrepareWeld_XPBD(s2Joint* base, s2StepContext* context)
{
	s2PrepareWeld(base, context, false);
}

void s2WarmStartWeld(s2Joint* base, s2StepContext* context)
{
	s2WeldJoint* joint = &base->weldJoint;

	s2Body* bodyA = context->bodies + base->edges[0].bodyIndex;
	s2Body* bodyB = context->bodies + base->edges[1].bodyIndex;

	s2Vec2 rA = s2RotateVector(bodyA->rot, joint->localAnchorA);
	s2Vec2 rB = s2RotateVector(bodyB->rot, joint->localAnchorB);
	s2Vec2 P = joint->linearImpulse;

	bodyA->linearVelocity = s2MulSub(bodyA->linearVelocity, joint->invMassA, P);
	bodyA->angularVelocity -= joint->invIA * (s2Cross(rA, P) + joint->angularImpulse);
	bodyB->linearVelocity = s2MulAdd(bodyB->linearVelocity, joint->invMassB, P);
	bodyB->angularVelocity += joint->invIB * (s2Cross(rB, P) + joint->angularImpulse);
}

// Rigid with a zero bias rate. Soft with the soft coefficients.
static void s2SolveWeldVelocity(s2Joint* base, s2StepContext* context, float biasRate, float massScale, float impulseScale)
{
	s2WeldJoint* joint = &base->weldJoint;

	s2Body* bodyA = context->bodies + base->edges[0].bodyIndex;
	s2Body* bodyB = context->bodies + base->edges[1].bodyIndex;

	s2Vec2 vA = bodyA->linearVelocity;
	float wA = bodyA->angularVelocity;
	s2Vec2 vB = bodyB->linearVelocity;
	float wB = bodyB->angularVelocity;

	float mA = joint->invMassA, mB = joint->invMassB;
	float iA = joint->invIA, iB = joint->invIB;

	{
		float C = s2RelativeAngle(bodyB->rot, bodyA->rot) - joint->referenceAngle;
		float Cdot = wB - wA;
		float impulse = -massScale * joint->axialMass * (Cdot + biasRate * C) - impulseScale * joint->angularImpulse;
		joint->angularImpulse += impulse;

		wA -= iA * impulse;
		wB += iB * impulse;
	}

	{
		s2Vec2 rA = s2RotateVector(bodyA->rot, joint->localAnchorA);
		s2Vec2 rB = s2RotateVector(bodyB->rot, joint->localAnchorB);

		s2Vec2 dcA = bodyA->deltaPosition;
		s2Vec2 dcB = bodyB->deltaPosition;
		s2Vec2 C = s2Add(s2Add(s2Sub(dcB, dcA), s2Sub(rB, rA)), joint->centerDiff0);
		s2Vec2 Cdot = s2Sub(s2Add(vB, s2CrossSV(wB, rB)), s2Add(vA, s2CrossSV(wA, rA)));

		s2Mat22 K;
		K.cx.x = mA + mB + rA.y * rA.y * iA + rB.y * rB.y * iB;
		K.cy.x = -rA.y * rA.x * iA - rB.y * rB.x * iB;
		K.cx.y = K.cy.x;
		K.cy.y = mA + mB + rA.x * rA.x * iA + rB.x * rB.x * iB;
		s2Vec2 b = s2Solve22(K, s2MulAdd(Cdot, biasRate, C));

		s2Vec2 impulse;
		impulse.x = -massScale * b.x - impulseScale * joint->linearImpulse.x;
		impulse.y = -massScale * b.y - impulseScale * joint->linearImpulse.y;
		joint->linearImpulse = s2Add(joint->linearImpulse, impulse);

		vA = s2MulSub(vA, mA, impulse);
		wA -= iA * s2Cross(rA, impulse);
		vB = s2MulAdd(vB, mB, impulse);
		wB += iB * s2Cross(rB, impulse);
	}

	bodyA->linearVelocity = vA;
	bodyA->angularVelocity = wA;
	bodyB->linearVelocity = vB;
	bodyB->angularVelocity = wB;
}

void s2SolveWeld(s2Joint* base, s2StepContext* context)
{
	s2SolveWeldVelocity(base, context, 0.0f, 1.0f, 0.0f);
}

void s2SolveWeld_Soft(s2Joint* base, s2StepContext* context, bool useBias)
{
	if (useBias)
	{
		s2WeldJoint* joint = &base->weldJoint;
		s2SolveWeldVelocity(base, context, joint->biasCoefficient, joint->massCoefficient, joint->impulseCoefficient);
	}
	else
	{
		s2SolveWeldVelocity(base, context, 0.0f, 1.0f, 0.0f);
	}
}

void s2SolveWeld_Baumgarte(s2Joint* base, s2StepContext* context, float inv_h, bool useBias)
{
	s2SolveWeldVelocity(base, context, useBias ? s2_baumgarte * inv_h : 0.0f, 1.0f, 0.0f);
}

void s2SolveWeldPosition(s2Joint* base, s2StepContext* context)
{
	s2WeldJoint* joint = &base->weldJoint;

	s2Body* bodyA = context->bodies + base->edges[0].bodyIndex;
	s2Body* bodyB = context->bodies + base->edges[1].bodyIndex;

	s2Vec2 dcA = bodyA->deltaPosition;
	s2Rot qA = bodyA->rot;
	s2Vec2 dcB = bodyB->deltaPosition;
	s2Rot qB = bodyB->rot;

	float mA = joint->invMassA, mB = joint->invMassB;
	float iA = joint->invIA, iB = joint->invIB;

	{
		float C = s2RelativeAngle(qB, qA) - joint->referenceAngle;
		float impulse = -joint->axialMass * C;
		qA = s2IntegrateRot(qA, -iA * impulse);
		qB = s2IntegrateRot(qB, iB * impulse);
	}

	{
		s2Vec2 rA = s2RotateVector(qA, joint->localAnchorA);
		s2Vec2 rB = s2RotateVector(qB, joint->localAnchorB);
		s2Vec2 C = s2Add(s2Add(s2Sub(dcB, dcA), s2Sub(rB, rA)), joint->centerDiff0);

		s2Mat22 K;
		K.cx.x = mA + mB + rA.y * rA.y * iA + rB.y * rB.y * iB;
		K.cy.x = -rA.y * rA.x * iA - rB.y * rB.x * iB;
		K.cx.y = K.cy.x;
		K.cy.y = mA + mB + rA.x * rA.x * iA + rB.x * rB.x * iB;
		s2Vec2 impulse = s2Solve22(K, s2Neg(C));

		dcA = s2MulSub(dcA, mA, impulse);
		qA = s2IntegrateRot(qA, -iA * s2Cross(rA, impulse));
		dcB = s2MulAdd(dcB, mB, impulse);
		qB = s2IntegrateRot(qB, iB * s2Cross(rB, impulse));
	}

	bodyA->deltaPosition = dcA;
	bodyA->rot = qA;
	bodyB->deltaPosition = dcB;
	bodyB->rot = qB;
}

void s2SolveWeld_XPBD(s2Joint* base, s2StepContext* context)
{
	s2WeldJoint* joint = &base->weldJoint;

	s2Body* bodyA = context->bodies + base->edges[0].bodyIndex;
	s2Body* bodyB = context->bodies + base->edges[1].bodyIndex;

	float mA = joint->invMassA, mB = joint->invMassB;
	float iA = joint->invIA, iB = joint->invIB;

	if (mA == 0.0f && mB == 0.0f)
	{
		return;
	}

	s2Vec2 dcA = bodyA->deltaPosition;
	s2Rot qA = bodyA->rot;
	s2Vec2 dcB = bodyB->deltaPosition;
	s2Rot qB = bodyB->rot;

	if (iA + iB > 0.0f)
	{
		float C = s2RelativeAngle(qB, qA) - joint->referenceAngle;
		float lambda = -C / (iA + iB);
		qA = s2IntegrateRot(qA, -iA * lambda);
		qB = s2IntegrateRot(qB, iB * lambda);
	}

	{
		s2Vec2 rA = s2RotateVector(qA, joint->localAnchorA);
		s2Vec2 rB = s2RotateVector(qB, joint->localAnchorB);
		s2Vec2 separation = s2Add(s2Add(s2Sub(dcB, dcA), s2Sub(rB, rA)), joint->centerDiff0);

		float c = s2Length(separation);
		s2Vec2 n = s2Normalize(separation);

		// w1 and w2 in the XPBD paper
		float rnA = s2Cross(rA, n);
		float rnB = s2Cross(rB, n);
		float kA = mA + iA * rnA * rnA;
		float kB = mB + iB * rnB * rnB;

		float lambda = -c / (kA + kB);
		s2Vec2 p = s2MulSV(lambda, n);

		dcA = s2MulSub(dcA, mA, p);
		qA = s2IntegrateRot(qA, -iA * s2Cross(rA, p));
		dcB = s2MulAdd(dcB, mB, p);
		qB = s2IntegrateRot(qB, iB * s2Cross(rB, p));
	}

	bodyA->deltaPosition = dcA;
	bodyA->rot = qA;
	bodyB->deltaPosition = dcB;
	bodyB->rot = qB;
}
