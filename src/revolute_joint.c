// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "body.h"
#include "core.h"
#include "joint.h"
#include "solver_data.h"
#include "world.h"

#include "solver2d/debug_draw.h"

// Point-to-point constraint
// C = p2 - p1
// Cdot = v2 - v1
//      = v2 + cross(w2, r2) - v1 - cross(w1, r1)
// J = [-I -r1_skew I r2_skew ]
// Identity used:
// w k % (rx i + ry j) = w * (-ry i + rx j)

// Motor constraint
// Cdot = w2 - w1
// J = [0 0 -1 0 0 1]
// K = invI1 + invI2

void s2InitializeRevolute(s2Joint* base, s2StepContext* context)
{
	S2_ASSERT(base->type == s2_revoluteJoint);

	int32_t indexA = base->edges[0].bodyIndex;
	int32_t indexB = base->edges[1].bodyIndex;
	S2_ASSERT(0 <= indexA && indexA < context->bodyCapacity);
	S2_ASSERT(0 <= indexB && indexB < context->bodyCapacity);

	s2Body* bodyA = context->bodies + indexA;
	s2Body* bodyB = context->bodies + indexB;
	S2_ASSERT(bodyA->object.index == bodyA->object.next);
	S2_ASSERT(bodyB->object.index == bodyB->object.next);

	s2RevoluteJoint* joint = &base->revoluteJoint;
	joint->localCenterA = bodyA->localCenter;
	joint->invMassA = bodyA->invMass;
	joint->invIA = bodyA->invI;

	joint->localCenterB = bodyB->localCenter;
	joint->invMassB = bodyB->invMass;
	joint->invIB = bodyB->invI;

	float aA = bodyA->angle;
	s2Vec2 vA = bodyA->linearVelocity;
	float wA = bodyA->angularVelocity;

	float aB = bodyB->angle;
	s2Vec2 vB = bodyB->linearVelocity;
	float wB = bodyB->angularVelocity;

	s2Rot qA = s2MakeRot(aA);
	s2Rot qB = s2MakeRot(aB);

	joint->rA = s2RotateVector(qA, s2Sub(base->localAnchorA, joint->localCenterA));
	joint->rB = s2RotateVector(qB, s2Sub(base->localAnchorB, joint->localCenterB));

	// J = [-I -r1_skew I r2_skew]
	// r_skew = [-ry; rx]

	// Matlab
	// K = [ mA+r1y^2*iA+mB+r2y^2*iB,  -r1y*iA*r1x-r2y*iB*r2x]
	//     [  -r1y*iA*r1x-r2y*iB*r2x, mA+r1x^2*iA+mB+r2x^2*iB]

	float mA = joint->invMassA, mB = joint->invMassB;
	float iA = joint->invIA, iB = joint->invIB;

	joint->K.cx.x = mA + mB + joint->rA.y * joint->rA.y * iA + joint->rB.y * joint->rB.y * iB;
	joint->K.cy.x = -joint->rA.y * joint->rA.x * iA - joint->rB.y * joint->rB.x * iB;
	joint->K.cx.y = joint->K.cy.x;
	joint->K.cy.y = mA + mB + joint->rA.x * joint->rA.x * iA + joint->rB.x * joint->rB.x * iB;

	joint->axialMass = iA + iB;
	bool fixedRotation;
	if (joint->axialMass > 0.0f)
	{
		joint->axialMass = 1.0f / joint->axialMass;
		fixedRotation = false;
	}
	else
	{
		fixedRotation = true;
	}

	joint->angle = aB - aA - joint->referenceAngle;
	if (joint->enableLimit == false || fixedRotation)
	{
		joint->lowerImpulse = 0.0f;
		joint->upperImpulse = 0.0f;
	}

	if (joint->enableMotor == false || fixedRotation)
	{
		joint->motorImpulse = 0.0f;
	}

	float axialImpulse = joint->motorImpulse + joint->lowerImpulse - joint->upperImpulse;
	s2Vec2 P = {joint->impulse.x, joint->impulse.y};

	vA = s2MulSub(vA, mA, P);
	wA -= iA * (s2Cross(joint->rA, P) + axialImpulse);

	vB = s2MulAdd(vB, mB, P);
	wB += iB * (s2Cross(joint->rB, P) + axialImpulse);

	bodyA->linearVelocity = vA;
	bodyA->angularVelocity = wA;
	bodyB->linearVelocity = vB;
	bodyB->angularVelocity = wB;
}

void s2SolveRevoluteVelocity(s2Joint* base, s2StepContext* context)
{
	S2_ASSERT(base->type == s2_revoluteJoint);

	s2RevoluteJoint* joint = &base->revoluteJoint;

	s2Body* bodyA = context->bodies + base->edges[0].bodyIndex;
	s2Body* bodyB = context->bodies + base->edges[1].bodyIndex;

	s2Vec2 vA = bodyA->linearVelocity;
	float wA = bodyA->angularVelocity;
	s2Vec2 vB = bodyB->linearVelocity;
	float wB = bodyB->angularVelocity;

	float mA = joint->invMassA, mB = joint->invMassB;
	float iA = joint->invIA, iB = joint->invIB;

	bool fixedRotation = (iA + iB == 0.0f);

	// Solve motor constraint.
	if (joint->enableMotor && fixedRotation == false)
	{
		float Cdot = wB - wA - joint->motorSpeed;
		float impulse = -joint->axialMass * Cdot;
		float oldImpulse = joint->motorImpulse;
		float maxImpulse = context->dt * joint->maxMotorTorque;
		joint->motorImpulse = S2_CLAMP(joint->motorImpulse + impulse, -maxImpulse, maxImpulse);
		impulse = joint->motorImpulse - oldImpulse;

		wA -= iA * impulse;
		wB += iB * impulse;
	}

	if (joint->enableLimit && fixedRotation == false)
	{
		// Lower limit
		{
			float C = joint->angle - joint->lowerAngle;
			float Cdot = wB - wA;
			float impulse = -joint->axialMass * (Cdot + S2_MAX(C, 0.0f) * context->inv_dt);
			float oldImpulse = joint->lowerImpulse;
			joint->lowerImpulse = S2_MAX(joint->lowerImpulse + impulse, 0.0f);
			impulse = joint->lowerImpulse - oldImpulse;

			wA -= iA * impulse;
			wB += iB * impulse;
		}

		// Upper limit
		// Note: signs are flipped to keep C positive when the constraint is satisfied.
		// This also keeps the impulse positive when the limit is active.
		{
			float C = joint->upperAngle - joint->angle;
			float Cdot = wA - wB;
			float impulse = -joint->axialMass * (Cdot + S2_MAX(C, 0.0f) * context->inv_dt);
			float oldImpulse = joint->upperImpulse;
			joint->upperImpulse = S2_MAX(joint->upperImpulse + impulse, 0.0f);
			impulse = joint->upperImpulse - oldImpulse;

			wA += iA * impulse;
			wB -= iB * impulse;
		}
	}

	// Solve point-to-point constraint
	{
		s2Vec2 Cdot = s2Sub(s2Add(vB, s2CrossSV(wB, joint->rB)), s2Add(vA, s2CrossSV(wA, joint->rA)));
		s2Vec2 impulse = s2Solve22(joint->K, s2Neg(Cdot));

		joint->impulse.x += impulse.x;
		joint->impulse.y += impulse.y;

		vA = s2MulSub(vA, mA, impulse);
		wA -= iA * s2Cross(joint->rA, impulse);

		vB = s2MulAdd(vB, mB, impulse);
		wB += iB * s2Cross(joint->rB, impulse);
	}

	bodyA->linearVelocity = vA;
	bodyA->angularVelocity = wA;
	bodyB->linearVelocity = vB;
	bodyB->angularVelocity = wB;
}

bool s2SolveRevolutePosition(s2Joint* base, s2StepContext* context)
{
	S2_ASSERT(base->type == s2_revoluteJoint);

	s2RevoluteJoint* joint = &base->revoluteJoint;

	s2Body* bodyA = context->bodies + base->edges[0].bodyIndex;
	s2Body* bodyB = context->bodies + base->edges[1].bodyIndex;

	s2Vec2 cA = bodyA->position;
	float aA = bodyA->angle;
	s2Vec2 cB = bodyB->position;
	float aB = bodyB->angle;

	s2Rot qA = s2MakeRot(aA), qB = s2MakeRot(aB);

	float angularError = 0.0f;
	float positionError = 0.0f;

	bool fixedRotation = (joint->invIA + joint->invIB == 0.0f);

	// Solve angular limit constraint
	if (joint->enableLimit && fixedRotation == false)
	{
		float angle = aB - aA - joint->referenceAngle;
		float C = 0.0f;

		if (S2_ABS(joint->upperAngle - joint->lowerAngle) < 2.0f * s2_angularSlop)
		{
			// Prevent large angular corrections
			C = S2_CLAMP(angle - joint->lowerAngle, -s2_maxAngularCorrection, s2_maxAngularCorrection);
		}
		else if (angle <= joint->lowerAngle)
		{
			// Prevent large angular corrections and allow some slop.
			C = S2_CLAMP(angle - joint->lowerAngle + s2_angularSlop, -s2_maxAngularCorrection, 0.0f);
		}
		else if (angle >= joint->upperAngle)
		{
			// Prevent large angular corrections and allow some slop.
			C = S2_CLAMP(angle - joint->upperAngle - s2_angularSlop, 0.0f, s2_maxAngularCorrection);
		}

		float limitImpulse = -joint->axialMass * C;
		aA -= joint->invIA * limitImpulse;
		aB += joint->invIB * limitImpulse;
		angularError = S2_ABS(C);
	}

	// Solve point-to-point constraint.
	{
		qA = s2MakeRot(aA);
		qB = s2MakeRot(aB);
		s2Vec2 rA = s2RotateVector(qA, s2Sub(base->localAnchorA, joint->localCenterA));
		s2Vec2 rB = s2RotateVector(qB, s2Sub(base->localAnchorB, joint->localCenterB));

		s2Vec2 C = s2Sub(s2Add(cB, rB), s2Add(cA, rA));
		positionError = s2Length(C);

		float mA = joint->invMassA, mB = joint->invMassB;
		float iA = joint->invIA, iB = joint->invIB;

		s2Mat22 K;
		K.cx.x = mA + mB + iA * rA.y * rA.y + iB * rB.y * rB.y;
		K.cx.y = -iA * rA.x * rA.y - iB * rB.x * rB.y;
		K.cy.x = K.cx.y;
		K.cy.y = mA + mB + iA * rA.x * rA.x + iB * rB.x * rB.x;

		s2Vec2 impulse = s2Solve22(K, s2Neg(C));

		cA = s2MulSub(cA, mA, impulse);
		aA -= iA * s2Cross(rA, impulse);

		cB = s2MulAdd(cB, mB, impulse);
		aB += iB * s2Cross(rB, impulse);
	}

	bodyA->position = cA;
	bodyA->angle = aA;
	bodyB->position = cB;
	bodyB->angle = aB;

	return positionError <= s2_linearSlop && angularError <= s2_angularSlop;
}

void s2RevoluteJoint_EnableLimit(s2JointId jointId, bool enableLimit)
{
	s2World* world = s2GetWorldFromIndex(jointId.world);

	S2_ASSERT(0 <= jointId.index && jointId.index < world->jointPool.capacity);

	s2Joint* joint = world->joints + jointId.index;
	S2_ASSERT(joint->object.index == joint->object.next);
	S2_ASSERT(joint->object.revision == jointId.revision);
	S2_ASSERT(joint->type == s2_revoluteJoint);
	joint->revoluteJoint.enableLimit = enableLimit;
}

void s2RevoluteJoint_EnableMotor(s2JointId jointId, bool enableMotor)
{
	s2World* world = s2GetWorldFromIndex(jointId.world);

	S2_ASSERT(0 <= jointId.index && jointId.index < world->jointPool.capacity);

	s2Joint* joint = world->joints + jointId.index;
	S2_ASSERT(joint->object.index == joint->object.next);
	S2_ASSERT(joint->object.revision == jointId.revision);
	S2_ASSERT(joint->type == s2_revoluteJoint);
	joint->revoluteJoint.enableMotor = enableMotor;
}

void s2RevoluteJoint_SetMotorSpeed(s2JointId jointId, float motorSpeed)
{
	s2World* world = s2GetWorldFromIndex(jointId.world);

	S2_ASSERT(0 <= jointId.index && jointId.index < world->jointPool.capacity);

	s2Joint* joint = world->joints + jointId.index;
	S2_ASSERT(joint->object.index == joint->object.next);
	S2_ASSERT(joint->object.revision == jointId.revision);
	S2_ASSERT(joint->type == s2_revoluteJoint);
	joint->revoluteJoint.motorSpeed = motorSpeed;
}

float s2RevoluteJoint_GetMotorTorque(s2JointId jointId, float inverseTimeStep)
{
	s2World* world = s2GetWorldFromIndex(jointId.world);

	S2_ASSERT(0 <= jointId.index && jointId.index < world->jointPool.capacity);

	s2Joint* joint = world->joints + jointId.index;
	S2_ASSERT(joint->object.index == joint->object.next);
	S2_ASSERT(joint->object.revision == jointId.revision);
	S2_ASSERT(joint->type == s2_revoluteJoint);
	return inverseTimeStep * joint->revoluteJoint.motorImpulse;
}

#if 0
void s2RevoluteJoint::Dump()
{
	int32 indexA = joint->bodyA->joint->islandIndex;
	int32 indexB = joint->bodyB->joint->islandIndex;

	s2Dump("  s2RevoluteJointDef jd;\n");
	s2Dump("  jd.bodyA = bodies[%d];\n", indexA);
	s2Dump("  jd.bodyB = bodies[%d];\n", indexB);
	s2Dump("  jd.collideConnected = bool(%d);\n", joint->collideConnected);
	s2Dump("  jd.localAnchorA.Set(%.9g, %.9g);\n", joint->localAnchorA.x, joint->localAnchorA.y);
	s2Dump("  jd.localAnchorB.Set(%.9g, %.9g);\n", joint->localAnchorB.x, joint->localAnchorB.y);
	s2Dump("  jd.referenceAngle = %.9g;\n", joint->referenceAngle);
	s2Dump("  jd.enableLimit = bool(%d);\n", joint->enableLimit);
	s2Dump("  jd.lowerAngle = %.9g;\n", joint->lowerAngle);
	s2Dump("  jd.upperAngle = %.9g;\n", joint->upperAngle);
	s2Dump("  jd.enableMotor = bool(%d);\n", joint->enableMotor);
	s2Dump("  jd.motorSpeed = %.9g;\n", joint->motorSpeed);
	s2Dump("  jd.maxMotorTorque = %.9g;\n", joint->maxMotorTorque);
	s2Dump("  joints[%d] = joint->world->CreateJoint(&jd);\n", joint->index);
}
#endif

void s2DrawRevolute(s2DebugDraw* draw, s2Joint* base, s2Body* bodyA, s2Body* bodyB)
{
	S2_ASSERT(base->type == s2_revoluteJoint);

	s2RevoluteJoint* joint = &base->revoluteJoint;

	s2Transform xfA = bodyA->transform;
	s2Transform xfB = bodyB->transform;
	s2Vec2 pA = s2TransformPoint(xfA, base->localAnchorA);
	s2Vec2 pB = s2TransformPoint(xfB, base->localAnchorB);

	s2Color c1 = {0.7f, 0.7f, 0.7f, 1.0f};
	s2Color c2 = {0.3f, 0.9f, 0.3f, 1.0f};
	s2Color c3 = {0.9f, 0.3f, 0.3f, 1.0f};
	s2Color c4 = {0.3f, 0.3f, 0.9f, 1.0f};
	s2Color c5 = {0.4f, 0.4f, 0.4f, 1.0f};

	draw->DrawPoint(pA, 5.0f, c4, draw->context);
	draw->DrawPoint(pB, 5.0f, c5, draw->context);

	float aA = bodyA->angle;
	float aB = bodyB->angle;
	float angle = aB - aA - joint->referenceAngle;

	const float L = 0.5f;

	s2Vec2 r = {L * cosf(angle), L * sinf(angle)};
	draw->DrawSegment(pB, s2Add(pB, r), c1, draw->context);
	draw->DrawCircle(pB, L, c1, draw->context);

	if (joint->enableLimit)
	{
		s2Vec2 rlo = {L * cosf(joint->lowerAngle), L * sinf(joint->lowerAngle)};
		s2Vec2 rhi = {L * cosf(joint->upperAngle), L * sinf(joint->upperAngle)};

		draw->DrawSegment(pB, s2Add(pB, rlo), c2, draw->context);
		draw->DrawSegment(pB, s2Add(pB, rhi), c3, draw->context);
	}

	s2Color color = {0.5f, 0.8f, 0.8f, 1.0f};
	draw->DrawSegment(xfA.p, pA, color, draw->context);
	draw->DrawSegment(pA, pB, color, draw->context);
	draw->DrawSegment(xfB.p, pB, color, draw->context);
}
