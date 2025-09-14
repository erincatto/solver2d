// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "body.h"
#include "core.h"
#include "joint.h"
#include "solvers.h"
#include "world.h"

#include "solver2d/debug_draw.h"

#include <assert.h>

// stale mass can slightly increase the jitter on ragdolls far from the origin
#define S2_FRESH_PIVOT_MASS 1

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

void s2PrepareRevolute(s2Joint* base, s2StepContext* context, bool warmStart)
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
	joint->localAnchorA = s2Sub(base->localOriginAnchorA, bodyA->localCenter);
	joint->invMassA = bodyA->invMass;
	joint->invIA = joint->inertiaScale * bodyA->invI;

	joint->localAnchorB = s2Sub(base->localOriginAnchorB, bodyB->localCenter);
	joint->invMassB = bodyB->invMass;
	joint->invIB = joint->inertiaScale * bodyB->invI;

	joint->centerDiff0 = s2Sub(bodyB->position, bodyA->position);

	s2Rot qA = bodyA->rot;
	s2Rot qB = bodyB->rot;

	s2Vec2 rA = s2RotateVector(qA, joint->localAnchorA);
	s2Vec2 rB = s2RotateVector(qB, joint->localAnchorB);

	// J = [-I -r1_skew I r2_skew]
	// r_skew = [-ry; rx]

	// Matlab
	// K = [ mA+r1y^2*iA+mB+r2y^2*iB,  -r1y*iA*r1x-r2y*iB*r2x]
	//     [  -r1y*iA*r1x-r2y*iB*r2x, mA+r1x^2*iA+mB+r2x^2*iB]

	float mA = joint->invMassA, mB = joint->invMassB;
	float iA = joint->invIA, iB = joint->invIB;

	s2Mat22 K;
	K.cx.x = mA + mB + rA.y * rA.y * iA + rB.y * rB.y * iB;
	K.cy.x = -rA.y * rA.x * iA - rB.y * rB.x * iB;
	K.cx.y = K.cy.x;
	K.cy.y = mA + mB + rA.x * rA.x * iA + rB.x * rB.x * iB;
	joint->pivotMass = s2GetInverse22(K);

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

	if (joint->enableLimit == false || fixedRotation || warmStart == false)
	{
		joint->lowerImpulse = 0.0f;
		joint->upperImpulse = 0.0f;
	}

	if (joint->enableMotor == false || fixedRotation || warmStart == false)
	{
		joint->motorImpulse = 0.0f;
	}

	if (warmStart == false)
	{
		joint->impulse = s2Vec2_zero;
	}
}

void s2WarmStartRevolute(s2Joint* base, s2StepContext* context)
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

	s2Rot qA = bodyA->rot;
	s2Vec2 vA = bodyA->linearVelocity;
	float wA = bodyA->angularVelocity;

	s2Rot qB = bodyB->rot;
	s2Vec2 vB = bodyB->linearVelocity;
	float wB = bodyB->angularVelocity;

	s2Vec2 rA = s2RotateVector(qA, joint->localAnchorA);
	s2Vec2 rB = s2RotateVector(qB, joint->localAnchorB);

	float mA = joint->invMassA, mB = joint->invMassB;
	float iA = joint->invIA, iB = joint->invIB;

	float axialImpulse = joint->motorImpulse + joint->lowerImpulse - joint->upperImpulse;
	s2Vec2 P = {joint->impulse.x, joint->impulse.y};

	vA = s2MulSub(vA, mA, P);
	wA -= iA * (s2Cross(rA, P) + axialImpulse);

	vB = s2MulAdd(vB, mB, P);
	wB += iB * (s2Cross(rB, P) + axialImpulse);

	bodyA->linearVelocity = vA;
	bodyA->angularVelocity = wA;
	bodyB->linearVelocity = vB;
	bodyB->angularVelocity = wB;
}

void s2SolveRevolute(s2Joint* base, s2StepContext* context, float h)
{
	S2_ASSERT(base->type == s2_revoluteJoint);

	s2RevoluteJoint* joint = &base->revoluteJoint;

	s2Body* bodyA = context->bodies + base->edges[0].bodyIndex;
	s2Body* bodyB = context->bodies + base->edges[1].bodyIndex;

	s2Rot qA = bodyA->rot;
	s2Rot qB = bodyB->rot;

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
		float maxImpulse = h * joint->maxMotorTorque;
		joint->motorImpulse = S2_CLAMP(joint->motorImpulse + impulse, -maxImpulse, maxImpulse);
		impulse = joint->motorImpulse - oldImpulse;

		wA -= iA * impulse;
		wB += iB * impulse;
	}

	if (joint->enableLimit && fixedRotation == false)
	{
		float angle = s2RelativeAngle(qB, qA) - joint->referenceAngle;

		// Lower limit
		{
			float C = angle - joint->lowerAngle;
			float Cdot = wB - wA;
			float impulse = -joint->axialMass * (Cdot + S2_MAX(C, 0.0f) / h);
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
			float C = joint->upperAngle - angle;
			float Cdot = wA - wB;
			float impulse = -joint->axialMass * (Cdot + S2_MAX(C, 0.0f) / h);
			float oldImpulse = joint->upperImpulse;
			joint->upperImpulse = S2_MAX(joint->upperImpulse + impulse, 0.0f);
			impulse = joint->upperImpulse - oldImpulse;

			wA += iA * impulse;
			wB -= iB * impulse;
		}
	}

	// Solve point-to-point constraint
	{
		s2Vec2 rA = s2RotateVector(qA, joint->localAnchorA);
		s2Vec2 rB = s2RotateVector(qB, joint->localAnchorB);

#if 0
		// geometric stiffness
		// J = [-I -r1_skew I r2_skew ]
		// JT * lambda = [-lambda;
		// dJT/dx * lambda = -skew(lambda) * skew(r)
		//
		float kA = fabsf(s2Cross(joint->impulse, rA));
		float kB = fabsf(s2Cross(joint->impulse, rB));

		//float iAt = 1.0f / (1.0f / iA + h * kA);
		//float iBt = 1.0f / (1.0f / iB + h * kB);
		float iAt = iA;
		float iBt = iB;

		if (h * kA * iA > 4.0f)
		{
			iAt = 1.0f / (1.0f / iA + 0.5f * (h * kA - 4.0f / iA) );
		}

		if (h * kB * iB > 4.0f)
		{
			iBt = 1.0f / (1.0f / iB + 0.5f * (h * kB - 4.0f / iB) );
		}

		if (iAt > iA)
		{
			iAt += 0.0f;
		}

		if (iBt > iB)
		{
			iBt += 0.0f;
		}

		s2Mat22 K;
		K.cx.x = mA + mB + rA.y * rA.y * iAt + rB.y * rB.y * iBt;
		K.cy.x = -rA.y * rA.x * iAt - rB.y * rB.x * iBt;
		K.cx.y = K.cy.x;
		K.cy.y = mA + mB + rA.x * rA.x * iAt + rB.x * rB.x * iBt;
		joint->pivotMass = s2GetInverse22(K);
#else
		float iAt = iA;
		float iBt = iB;
#endif

		s2Vec2 Cdot = s2Sub(s2Add(vB, s2CrossSV(wB, rB)), s2Add(vA, s2CrossSV(wA, rA)));
		s2Vec2 impulse = s2MulMV(joint->pivotMass, s2Neg(Cdot));

		joint->impulse.x += impulse.x;
		joint->impulse.y += impulse.y;

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

void s2SolveRevolutePosition(s2Joint* base, s2StepContext* context)
{
	S2_ASSERT(base->type == s2_revoluteJoint);

	s2RevoluteJoint* joint = &base->revoluteJoint;

	s2Body* bodyA = context->bodies + base->edges[0].bodyIndex;
	s2Body* bodyB = context->bodies + base->edges[1].bodyIndex;

	s2Vec2 dcA = bodyA->deltaPosition;
	s2Rot qA = bodyA->rot;
	s2Vec2 dcB = bodyB->deltaPosition;
	s2Rot qB = bodyB->rot;

	bool fixedRotation = (joint->invIA + joint->invIB == 0.0f);

	// Solve angular limit constraint
	if (joint->enableLimit && fixedRotation == false)
	{
		float angle = s2RelativeAngle(qB, qA) - joint->referenceAngle;
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
		qA = s2IntegrateRot(qA, -joint->invIA * limitImpulse);
		qB = s2IntegrateRot(qB, joint->invIB * limitImpulse);
	}

	// Solve point-to-point constraint.
	for (int i = 0; i < 1; ++i)
	{
		s2Vec2 rA = s2RotateVector(qA, joint->localAnchorA);
		s2Vec2 rB = s2RotateVector(qB, joint->localAnchorB);

		s2Vec2 C = s2Add(s2Add(s2Sub(dcB, dcA), s2Sub(rB, rA)), joint->centerDiff0);

		//float largeTol = 100.0f * s2_linearSlop;
		//joint->largeError = C.x < -largeTol | largeTol < C.x | C.y < -largeTol | largeTol < C.y;
		//if (joint->largeError)
		//{
		//	bodyA->largeError = true;
		//	bodyB->largeError = true;
		//}

		float lengthC = s2Length(C);
		float threshold = 50.0f * s2_linearSlop;
		float inertiaScale = 1.0f;
		if (lengthC > threshold)
		{
			float s = threshold / lengthC;
			joint->impulse.x *= powf(s, 0.25f);
			joint->impulse.y *= powf(s, 0.25f);
			inertiaScale = s * s;
		}

		// smoother
		joint->inertiaScale = 0.9f * joint->inertiaScale + 0.1f * inertiaScale;

		//joint->impulse.x *= powf(joint->inertiaScale, 0.25f);
		//joint->impulse.y *= powf(joint->inertiaScale, 0.25f);

		//float tol = 0.5f * s2_linearSlop;
		//if (i > 0 & -tol < C.x & C.x < tol & -tol < C.y & C.y < tol)
		//{
		//	break;
		//}

		float mA = joint->invMassA, mB = joint->invMassB;
		float iA = joint->inertiaScale * bodyA->invI, iB = joint->inertiaScale * bodyB->invI;

#if S2_FRESH_PIVOT_MASS == 1
		s2Mat22 K;
		K.cx.x = mA + mB + iA * rA.y * rA.y + iB * rB.y * rB.y;
		K.cx.y = -iA * rA.x * rA.y - iB * rB.x * rB.y;
		K.cy.x = K.cx.y;
		K.cy.y = mA + mB + iA * rA.x * rA.x + iB * rB.x * rB.x;
		s2Vec2 impulse = s2Solve22(K, s2Neg(C));
#else
		s2Vec2 impulse = s2MulMV(joint->pivotMass, s2Neg(C));
#endif

		dcA = s2MulSub(dcA, mA, impulse);
		qA = s2IntegrateRot(qA, -iA * s2Cross(rA, impulse));

		dcB = s2MulAdd(dcB, mB, impulse);
		qB = s2IntegrateRot(qB, iB * s2Cross(rB, impulse));


		//rA = s2RotateVector(qA, joint->localAnchorA);
		//rB = s2RotateVector(qB, joint->localAnchorB);
		//s2Vec2 C2 = s2Add(s2Add(s2Sub(dcB, dcA), s2Sub(rB, rA)), joint->centerDiff0);
		//C2.x += 0.0f;
		//(void)C2;
	}

	bodyA->deltaPosition = dcA;
	bodyA->rot = qA;
	bodyB->deltaPosition = dcB;
	bodyB->rot = qB;
}

void s2PrepareRevolute_Soft(s2Joint* base, s2StepContext* context, float h, float hertz, bool warmStart)
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
	joint->localAnchorA = s2Sub(base->localOriginAnchorA, bodyA->localCenter);
	joint->invMassA = bodyA->invMass;
	joint->invIA = bodyA->invI;

	joint->localAnchorB = s2Sub(base->localOriginAnchorB, bodyB->localCenter);
	joint->invMassB = bodyB->invMass;
	joint->invIB = bodyB->invI;

	joint->centerDiff0 = s2Sub(bodyB->position, bodyA->position);

	// J = [-I -r1_skew I r2_skew]
	// r_skew = [-ry; rx]

	// Matlab
	// K = [ mA+r1y^2*iA+mB+r2y^2*iB,  -r1y*iA*r1x-r2y*iB*r2x]
	//     [  -r1y*iA*r1x-r2y*iB*r2x, mA+r1x^2*iA+mB+r2x^2*iB]

	float mA = joint->invMassA, mB = joint->invMassB;
	float iA = joint->invIA, iB = joint->invIB;

	s2Rot qA = bodyA->rot;
	s2Rot qB = bodyB->rot;

	s2Vec2 rA = s2RotateVector(qA, joint->localAnchorA);
	s2Vec2 rB = s2RotateVector(qB, joint->localAnchorB);

	s2Mat22 K;
	K.cx.x = mA + mB + rA.y * rA.y * iA + rB.y * rB.y * iB;
	K.cy.x = -rA.y * rA.x * iA - rB.y * rB.x * iB;
	K.cx.y = K.cy.x;
	K.cy.y = mA + mB + rA.x * rA.x * iA + rB.x * rB.x * iB;

	joint->pivotMass = s2GetInverse22(K);

	{
		const float zeta = 1.0f;
		float omega = 2.0f * s2_pi * hertz;
		joint->biasCoefficient = omega / (2.0f * zeta + h * omega);
		float c = h * omega * (2.0f * zeta + h * omega);
		joint->impulseCoefficient = 1.0f / (1.0f + c);
		joint->massCoefficient = c * joint->impulseCoefficient;
	}

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

	if (joint->enableLimit == false || fixedRotation || warmStart == false)
	{
		joint->lowerImpulse = 0.0f;
		joint->upperImpulse = 0.0f;
	}

	if (joint->enableMotor == false || fixedRotation || warmStart == false)
	{
		joint->motorImpulse = 0.0f;
	}

	if (warmStart == false)
	{
		joint->impulse = s2Vec2_zero;
	}
}

void s2SolveRevolute_Soft(s2Joint* base, s2StepContext* context, float h, float inv_h, bool useBias)
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
		float maxImpulse = h * joint->maxMotorTorque;
		joint->motorImpulse = S2_CLAMP(joint->motorImpulse + impulse, -maxImpulse, maxImpulse);
		impulse = joint->motorImpulse - oldImpulse;

		wA -= iA * impulse;
		wB += iB * impulse;
	}

	if (joint->enableLimit && fixedRotation == false)
	{
		float jointAngle = s2RelativeAngle(bodyB->rot, bodyA->rot) - joint->referenceAngle;

		// Lower limit
		{
			float C = jointAngle - joint->lowerAngle;
			float bias = 0.0f;
			float massScale = 1.0f;
			float impulseScale = 0.0f;
			if (C > 0.0f)
			{
				// speculation
				bias = C * inv_h;
			}
			else if (useBias)
			{
				bias = joint->biasCoefficient * C;
				massScale = joint->massCoefficient;
				impulseScale = joint->impulseCoefficient;
			}

			float Cdot = wB - wA;
			float impulse = -joint->axialMass * massScale * (Cdot + bias) - impulseScale * joint->lowerImpulse;
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
			float C = joint->upperAngle - jointAngle;

			float bias = 0.0f;
			float massScale = 1.0f;
			float impulseScale = 0.0f;
			if (C > 0.0f)
			{
				// speculation
				bias = C * inv_h;
			}
			else if (useBias)
			{
				bias = joint->biasCoefficient * C;
				massScale = joint->massCoefficient;
				impulseScale = joint->impulseCoefficient;
			}

			float Cdot = wA - wB;
			float impulse = -joint->axialMass * massScale * (Cdot + bias) - impulseScale * joint->lowerImpulse;
			float oldImpulse = joint->upperImpulse;
			joint->upperImpulse = S2_MAX(joint->upperImpulse + impulse, 0.0f);
			impulse = joint->upperImpulse - oldImpulse;

			wA += iA * impulse;
			wB -= iB * impulse;
		}
	}

	// Solve point-to-point constraint
	{
		// Update anchors for TGS solvers.
		// Anchors are wastefully recomputed for PGS solvers or relax stages.
		s2Rot qA = bodyA->rot;
		s2Rot qB = bodyB->rot;
		s2Vec2 rA = s2RotateVector(qA, joint->localAnchorA);
		s2Vec2 rB = s2RotateVector(qB, joint->localAnchorB);

		s2Vec2 Cdot = s2Sub(s2Add(vB, s2CrossSV(wB, rB)), s2Add(vA, s2CrossSV(wA, rA)));

		s2Vec2 bias = s2Vec2_zero;
		float massScale = 1.0f;
		float impulseScale = 0.0f;
		if (useBias)
		{
			s2Vec2 dcA = bodyA->deltaPosition;
			s2Vec2 dcB = bodyB->deltaPosition;

			s2Vec2 separation = s2Add(s2Add(s2Sub(dcB, dcA), s2Sub(rB, rA)), joint->centerDiff0);
			bias = s2MulSV(joint->biasCoefficient, separation);
			massScale = joint->massCoefficient;
			impulseScale = joint->impulseCoefficient;
		}

#if S2_FRESH_PIVOT_MASS == 1
		s2Mat22 K;
		K.cx.x = mA + mB + rA.y * rA.y * iA + rB.y * rB.y * iB;
		K.cy.x = -rA.y * rA.x * iA - rB.y * rB.x * iB;
		K.cx.y = K.cy.x;
		K.cy.y = mA + mB + rA.x * rA.x * iA + rB.x * rB.x * iB;
		s2Vec2 b = s2Solve22(K, s2Add(Cdot, bias));
#else
		s2Vec2 b = s2MulMV(joint->pivotMass, s2Add(Cdot, bias));
#endif

		s2Vec2 impulse;
		impulse.x = -massScale * b.x - impulseScale * joint->impulse.x;
		impulse.y = -massScale * b.y - impulseScale * joint->impulse.y;
		joint->impulse.x += impulse.x;
		joint->impulse.y += impulse.y;

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

// similar to box2d_lite
void s2SolveRevolute_Baumgarte(s2Joint* base, s2StepContext* context, float h, float inv_h, bool useBias)
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
		float maxImpulse = h * joint->maxMotorTorque;
		joint->motorImpulse = S2_CLAMP(joint->motorImpulse + impulse, -maxImpulse, maxImpulse);
		impulse = joint->motorImpulse - oldImpulse;

		wA -= iA * impulse;
		wB += iB * impulse;
	}

	if (joint->enableLimit && fixedRotation == false)
	{
		float jointAngle = s2RelativeAngle(bodyB->rot, bodyA->rot) - joint->referenceAngle;

		// Lower limit
		{
			float C = jointAngle - joint->lowerAngle;
			float bias = 0.0f;
			if (C > 0.0f)
			{
				// speculation
				bias = C * inv_h;
			}
			else if (useBias)
			{
				bias = s2_baumgarte * inv_h * C;
			}

			float Cdot = wB - wA;
			float impulse = -joint->axialMass * (Cdot + bias);
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
			float C = joint->upperAngle - jointAngle;

			float bias = 0.0f;
			if (C > 0.0f)
			{
				// speculation
				bias = C * inv_h;
			}
			else if (useBias)
			{
				bias = s2_baumgarte * inv_h * C;
			}

			float Cdot = wA - wB;
			float impulse = -joint->axialMass * (Cdot + bias);
			float oldImpulse = joint->upperImpulse;
			joint->upperImpulse = S2_MAX(joint->upperImpulse + impulse, 0.0f);
			impulse = joint->upperImpulse - oldImpulse;

			wA += iA * impulse;
			wB -= iB * impulse;
		}
	}

	// Solve point-to-point constraint
	{
		s2Rot qA = bodyA->rot;
		s2Rot qB = bodyB->rot;
		s2Vec2 dcA = bodyA->deltaPosition;
		s2Vec2 dcB = bodyB->deltaPosition;

		s2Vec2 rA = s2RotateVector(qA, joint->localAnchorA);
		s2Vec2 rB = s2RotateVector(qB, joint->localAnchorB);


		s2Vec2 Cdot = s2Sub(s2Add(vB, s2CrossSV(wB, rB)), s2Add(vA, s2CrossSV(wA, rA)));

		s2Vec2 separation = s2Add(s2Add(s2Sub(dcB, dcA), s2Sub(rB, rA)), joint->centerDiff0);
		s2Vec2 bias = s2MulSV(s2_baumgarte * inv_h, separation);

#if S2_FRESH_PIVOT_MASS == 1
		s2Mat22 K;
		K.cx.x = mA + mB + rA.y * rA.y * iA + rB.y * rB.y * iB;
		K.cy.x = -rA.y * rA.x * iA - rB.y * rB.x * iB;
		K.cx.y = K.cy.x;
		K.cy.y = mA + mB + rA.x * rA.x * iA + rB.x * rB.x * iB;
		s2Vec2 b = s2Solve22(K, s2Add(Cdot, bias));
#else
		s2Vec2 b = s2MulMV(joint->pivotMass, s2Add(Cdot, bias));
#endif

		s2Vec2 impulse = {-b.x, -b.y};
		joint->impulse.x += impulse.x;
		joint->impulse.y += impulse.y;

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

void s2PrepareRevolute_XPBD(s2Joint* base, s2StepContext* context)
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
	joint->localAnchorA = s2Sub(base->localOriginAnchorA, bodyA->localCenter);
	joint->invMassA = bodyA->invMass;
	joint->invIA = bodyA->invI;

	joint->localAnchorB = s2Sub(base->localOriginAnchorB, bodyB->localCenter);
	joint->invMassB = bodyB->invMass;
	joint->invIB = bodyB->invI;

	joint->centerDiff0 = s2Sub(bodyB->position, bodyA->position);

	joint->pivotMass = (s2Mat22){0};
	joint->axialMass = 0.0f;
	joint->impulse = s2Vec2_zero;
	joint->lowerImpulse = 0.0f;
	joint->upperImpulse = 0.0f;
	joint->motorImpulse = 0.0f;
}

void s2SolveRevolute_XPBD(s2Joint* base, s2StepContext* context, float inv_h)
{
	S2_ASSERT(base->type == s2_revoluteJoint);

	// joint grid sample blows up (more quickly) without compliance
	// float compliance = 0.00001f * inv_h * inv_h;
	float compliance = 0.0f;

	s2RevoluteJoint* joint = &base->revoluteJoint;

	s2Body* bodyA = context->bodies + base->edges[0].bodyIndex;
	s2Body* bodyB = context->bodies + base->edges[1].bodyIndex;

	s2Vec2 dcA = bodyA->deltaPosition;
	s2Rot qA = bodyA->rot;
	s2Vec2 dcB = bodyB->deltaPosition;
	s2Rot qB = bodyB->rot;

	// Solve point-to-point constraint.
	{
		s2Vec2 rA = s2RotateVector(qA, joint->localAnchorA);
		s2Vec2 rB = s2RotateVector(qB, joint->localAnchorB);

		s2Vec2 separation = s2Add(s2Add(s2Sub(dcB, dcA), s2Sub(rB, rA)), joint->centerDiff0);

		float c = s2Length(separation);
		s2Vec2 n = s2Normalize(separation);

		float mA = joint->invMassA, mB = joint->invMassB;

		if (mA == 0.0f && mB == 0.0f)
		{
			// static connection
			return;
		}

		float iA = joint->invIA, iB = joint->invIB;

		float rnA = s2Cross(rA, n);
		float rnB = s2Cross(rB, n);

		// w1 and w2 in paper
		float kA = mA + iA * rnA * rnA;
		float kB = mB + iB * rnB * rnB;

		assert(kA + kB > 0.0f);

		// float lambda = -c / (kA + kB);
		float lambda = -c / (kA + kB + compliance);

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

void s2DrawRevolute(s2DebugDraw* draw, s2Joint* base, s2Body* bodyA, s2Body* bodyB)
{
	S2_ASSERT(base->type == s2_revoluteJoint);

	s2RevoluteJoint* joint = &base->revoluteJoint;

	s2Transform xfA = S2_TRANSFORM(bodyA);
	s2Transform xfB = S2_TRANSFORM(bodyB);
	s2Vec2 pA = s2TransformPoint(xfA, base->localOriginAnchorA);
	s2Vec2 pB = s2TransformPoint(xfB, base->localOriginAnchorB);

	s2Color c1 = {0.7f, 0.7f, 0.7f, 1.0f};
	s2Color c2 = {0.3f, 0.9f, 0.3f, 1.0f};
	s2Color c3 = {0.9f, 0.3f, 0.3f, 1.0f};
	s2Color c4 = {0.3f, 0.3f, 0.9f, 1.0f};
	s2Color c5 = {0.4f, 0.4f, 0.4f, 1.0f};

	draw->DrawPoint(pA, 5.0f, c4, draw->context);
	draw->DrawPoint(pB, 5.0f, c5, draw->context);

	s2Rot qA = bodyA->rot;
	s2Rot qB = bodyB->rot;

	const float L = base->drawSize;

	s2Vec2 r = s2RotateVector(qB, (s2Vec2){L * cosf(joint->referenceAngle), L * sinf(joint->referenceAngle)});
	draw->DrawSegment(pB, s2Add(pB, r), c1, draw->context);
	draw->DrawCircle(pB, L, c1, draw->context);

	if (joint->enableLimit)
	{
		s2Vec2 rlo = s2RotateVector(qA, (s2Vec2){L * cosf(joint->lowerAngle), L * sinf(joint->lowerAngle)});
		s2Vec2 rhi = s2RotateVector(qA, (s2Vec2){L * cosf(joint->upperAngle), L * sinf(joint->upperAngle)});

		draw->DrawSegment(pB, s2Add(pB, rlo), c2, draw->context);
		draw->DrawSegment(pB, s2Add(pB, rhi), c3, draw->context);
	}

	s2Color color = {0.5f, 0.8f, 0.8f, 1.0f};
	draw->DrawSegment(xfA.p, pA, color, draw->context);
	draw->DrawSegment(pA, pB, color, draw->context);
	draw->DrawSegment(xfB.p, pB, color, draw->context);
}
