// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "body.h"
#include "core.h"
#include "joint.h"
#include "solvers.h"
#include "world.h"

// p = attached point, m = mouse point
// C = p - m
// Cdot = v
//      = v + cross(w, r)
// J = [I r_skew]
// Identity used:
// w k % (rx i + ry mouse) = w * (-ry i + rx mouse)

void s2MouseJoint_SetTarget(s2JointId jointId, s2Vec2 target)
{
	s2World* world = s2GetWorldFromIndex(jointId.world);

	S2_ASSERT(0 <= jointId.index && jointId.index < world->jointPool.capacity);

	s2Joint* base = world->joints + jointId.index;
	S2_ASSERT(base->object.index == base->object.next);
	S2_ASSERT(base->object.revision == jointId.revision);
	S2_ASSERT(base->type == s2_mouseJoint);
	base->mouseJoint.targetA = target;
}

void s2PrepareMouse(s2Joint* base, s2StepContext* context)
{
	S2_ASSERT(base->type == s2_mouseJoint);

	int32_t indexB = base->edges[1].bodyIndex;
	S2_ASSERT(0 <= indexB && indexB < context->bodyCapacity);

	s2Body* bodyB = context->bodies + indexB;
	S2_ASSERT(bodyB->object.index == bodyB->object.next);

	float mB = bodyB->invMass;
	float iB = bodyB->invI;

	s2MouseJoint* joint = &base->mouseJoint;
	joint->localAnchorB = s2Sub(base->localOriginAnchorB, bodyB->localCenter);
	joint->invMassB = mB;
	joint->invIB = iB;


	{
		float h = context->dt;
		float zeta = joint->dampingRatio;
		float omega = 2.0f * s2_pi * joint->hertz;
		joint->biasCoefficient = omega / (2.0f * zeta + h * omega);
		float c = h * omega * (2.0f * zeta + h * omega);
		joint->impulseCoefficient = 1.0f / (1.0f + c);
		joint->massCoefficient = c * joint->impulseCoefficient;
	}

	// Compute the effective mass matrix.
	s2Rot qB = bodyB->rot;
	s2Vec2 rB = s2RotateVector(qB, joint->localAnchorB);

	// K    = [(1/m1 + 1/m2) * eye(2) - skew(r1) * invI1 * skew(r1) - skew(r2) * invI2 * skew(r2)]
	//      = [1/m1+1/m2     0    ] + invI1 * [r1.y*r1.y -r1.x*r1.y] + invI2 * [r1.y*r1.y -r1.x*r1.y]
	//        [    0     1/m1+1/m2]           [-r1.x*r1.y r1.x*r1.x]           [-r1.x*r1.y r1.x*r1.x]
	s2Mat22 K;
	K.cx.x = mB + iB * rB.y * rB.y;
	K.cx.y = -iB * rB.x * rB.y;
	K.cy.x = K.cx.y;
	K.cy.y = mB + iB * rB.x * rB.x;

	joint->pivotMass = s2GetInverse22(K);

	s2Vec2 cB = bodyB->position;
	joint->centerDiff0 = s2Sub(cB, joint->targetA);

	// disable motor warm start
	//joint->motorImpulse = 0.0f;

	// Cheat with some damping
	//wB *= S2_MAX(0.0f, 1.0f - 0.02f * (60.0f * h));
}

void s2WarmStartMouse(s2Joint* base, s2StepContext* context)
{
	S2_ASSERT(base->type == s2_mouseJoint);

	int32_t indexB = base->edges[1].bodyIndex;
	S2_ASSERT(0 <= indexB && indexB < context->bodyCapacity);

	s2Body* bodyB = context->bodies + indexB;
	S2_ASSERT(bodyB->object.index == bodyB->object.next);

	s2MouseJoint* joint = &base->mouseJoint;
	s2Rot qB = bodyB->rot;
	s2Vec2 rB = s2RotateVector(qB, joint->localAnchorB);

	s2Vec2 vB = bodyB->linearVelocity;
	float wB = bodyB->angularVelocity;

	vB = s2MulAdd(vB, joint->invMassB, joint->impulse);
	wB += joint->invIB * (s2Cross(rB, joint->impulse) + joint->motorImpulse);

	bodyB->linearVelocity = vB;
	bodyB->angularVelocity = wB;
}

void s2SolveMouse(s2Joint* base, s2StepContext* context)
{
	s2MouseJoint* joint = &base->mouseJoint;
	s2Body* bodyB = context->bodies + base->edges[1].bodyIndex;

	s2Vec2 vB = bodyB->linearVelocity;
	float wB = bodyB->angularVelocity;

	float mB = joint->invMassB;
	float iB = joint->invIB;

	{
		float impulse = -bodyB->I * wB;
		float oldImpulse = joint->motorImpulse;
		float maxImpulse = 0.00002f * bodyB->I * context->inv_h;
		joint->motorImpulse = S2_CLAMP(joint->motorImpulse + impulse, -maxImpulse, maxImpulse);
		impulse = joint->motorImpulse - oldImpulse;
		wB += iB * impulse;
	}

	{
		s2Rot qB = bodyB->rot;
		s2Vec2 rB = s2RotateVector(qB, joint->localAnchorB);
		s2Vec2 Cdot = s2Add(vB, s2CrossSV(wB, rB));

		s2Vec2 dcB = bodyB->deltaPosition;

		s2Vec2 separation = s2Add(s2Add(dcB, rB), joint->centerDiff0);
		s2Vec2 bias = s2MulSV(joint->biasCoefficient, separation);

		float massScale = joint->massCoefficient;
		float impulseScale = joint->impulseCoefficient;
		
		//s2Mat22 K;
		//K.cx.x = mB + iB * rB.y * rB.y;
		//K.cx.y = -iB * rB.x * rB.y;
		//K.cy.x = K.cx.y;
		//K.cy.y = mB + iB * rB.x * rB.x;
		//s2Vec2 b = s2Solve22(K, s2Add(Cdot, bias));

		s2Vec2 b = s2MulMV(joint->pivotMass, s2Add(Cdot, bias));

		s2Vec2 impulse;
		impulse.x = -massScale * b.x - impulseScale * joint->impulse.x;
		impulse.y = -massScale * b.y - impulseScale * joint->impulse.y;
		joint->impulse.x += impulse.x;
		joint->impulse.y += impulse.y;

		vB = s2MulAdd(vB, mB, impulse);
		wB += iB * s2Cross(rB, impulse);
	}

	bodyB->linearVelocity = vB;
	bodyB->angularVelocity = wB;
}

#if 0
void s2PrepareMouse_Soft(s2Joint* base, s2StepContext* context)
{
	S2_ASSERT(base->type == s2_mouseJoint);

	int32_t indexB = base->edges[1].bodyIndex;
	S2_ASSERT(0 <= indexB && indexB < context->bodyCapacity);

	s2Body* bodyB = context->bodies + indexB;
	S2_ASSERT(bodyB->object.index == bodyB->object.next);

	s2MouseJoint* joint = &base->mouseJoint;
	joint->localCenterB = bodyB->localCenter;
	joint->invMassB = bodyB->invMass;
	joint->invIB = bodyB->invI;

	s2Vec2 cB = bodyB->position;
	s2Rot qB = bodyB->rot;
	s2Vec2 vB = bodyB->linearVelocity;
	float wB = bodyB->angularVelocity;

	float d = joint->damping;
	float k = joint->stiffness;

	// magic formulas
	// gamma has units of inverse mass.
	// beta has units of inverse time.
	float h = context->dt;
	joint->gamma = h * (d + h * k);
	if (joint->gamma != 0.0f)
	{
		joint->gamma = 1.0f / joint->gamma;
	}
	joint->beta = h * k * joint->gamma;

	// Compute the effective mass matrix.
	joint->rB = s2RotateVector(qB, s2Sub(base->localOriginAnchorB, joint->localCenterB));

	// K    = [(1/m1 + 1/m2) * eye(2) - skew(r1) * invI1 * skew(r1) - skew(r2) * invI2 * skew(r2)]
	//      = [1/m1+1/m2     0    ] + invI1 * [r1.y*r1.y -r1.x*r1.y] + invI2 * [r1.y*r1.y -r1.x*r1.y]
	//        [    0     1/m1+1/m2]           [-r1.x*r1.y r1.x*r1.x]           [-r1.x*r1.y r1.x*r1.x]
	s2Mat22 K;
	K.cx.x = joint->invMassB + joint->invIB * joint->rB.y * joint->rB.y + joint->gamma;
	K.cx.y = -joint->invIB * joint->rB.x * joint->rB.y;
	K.cy.x = K.cx.y;
	K.cy.y = joint->invMassB + joint->invIB * joint->rB.x * joint->rB.x + joint->gamma;

	joint->pivotMass = s2GetInverse22(K);

	joint->C = s2Add(cB, s2Sub(joint->rB, joint->targetA));
	joint->C = s2MulSV(joint->beta, joint->C);

	// Cheat with some damping
	wB *= S2_MAX(0.0f, 1.0f - 0.02f * (60.0f * h));

	vB = s2MulAdd(vB, joint->invMassB, joint->impulse);
	wB += joint->invIB * s2Cross(joint->rB, joint->impulse);

	bodyB->linearVelocity = vB;
	bodyB->angularVelocity = wB;
}

void s2SolveMouse_Soft(s2Joint* base, s2StepContext* context, bool useBias)
{
	S2_MAYBE_UNUSED(useBias);

	s2MouseJoint* joint = &base->mouseJoint;
	s2Body* bodyB = context->bodies + base->edges[1].bodyIndex;

	s2Vec2 vB = bodyB->linearVelocity;
	float wB = bodyB->angularVelocity;

	// Cdot = v + cross(w, r)
	s2Vec2 Cdot = s2Add(vB, s2CrossSV(wB, joint->rB));
	s2Vec2 SoftCdot = s2Add(Cdot, s2MulAdd(joint->C, joint->gamma, joint->impulse));
	s2Vec2 impulse = s2Neg(s2MulMV(joint->pivotMass, SoftCdot));

	s2Vec2 oldImpulse = joint->impulse;
	joint->impulse = s2Add(joint->impulse, impulse);
	float maxImpulse = context->dt * joint->maxForce;
	if (s2LengthSquared(joint->impulse) > maxImpulse * maxImpulse)
	{
		joint->impulse = s2MulSV(maxImpulse / s2Length(joint->impulse), joint->impulse);
	}
	impulse = s2Sub(joint->impulse, oldImpulse);

	vB = s2MulAdd(vB, joint->invMassB, impulse);
	wB += joint->invIB * s2Cross(joint->rB, impulse);

	bodyB->linearVelocity = vB;
	bodyB->angularVelocity = wB;
}

void s2PrepareMouse_XPBD(s2Joint* base, s2StepContext* context)
{
	S2_ASSERT(base->type == s2_mouseJoint);

	int32_t indexB = base->edges[1].bodyIndex;
	S2_ASSERT(0 <= indexB && indexB < context->bodyCapacity);

	s2Body* bodyB = context->bodies + indexB;
	S2_ASSERT(bodyB->object.index == bodyB->object.next);

	s2MouseJoint* joint = &base->mouseJoint;
	joint->localCenterB = bodyB->localCenter;
	joint->invMassB = bodyB->invMass;
	joint->invIB = bodyB->invI;

	s2Vec2 cB = bodyB->position;
	s2Rot qB = bodyB->rot;
	s2Vec2 vB = bodyB->linearVelocity;
	float wB = bodyB->angularVelocity;

	float d = joint->damping;
	float k = joint->stiffness;

	// magic formulas
	// gamma has units of inverse mass.
	// beta has units of inverse time.
	float h = context->dt;
	joint->gamma = h * (d + h * k);
	if (joint->gamma != 0.0f)
	{
		joint->gamma = 1.0f / joint->gamma;
	}
	joint->beta = h * k * joint->gamma;

	// Compute the effective mass matrix.
	joint->rB = s2RotateVector(qB, s2Sub(base->localOriginAnchorB, joint->localCenterB));

	// K    = [(1/m1 + 1/m2) * eye(2) - skew(r1) * invI1 * skew(r1) - skew(r2) * invI2 * skew(r2)]
	//      = [1/m1+1/m2     0    ] + invI1 * [r1.y*r1.y -r1.x*r1.y] + invI2 * [r1.y*r1.y -r1.x*r1.y]
	//        [    0     1/m1+1/m2]           [-r1.x*r1.y r1.x*r1.x]           [-r1.x*r1.y r1.x*r1.x]
	s2Mat22 K;
	K.cx.x = joint->invMassB + joint->invIB * joint->rB.y * joint->rB.y + joint->gamma;
	K.cx.y = -joint->invIB * joint->rB.x * joint->rB.y;
	K.cy.x = K.cx.y;
	K.cy.y = joint->invMassB + joint->invIB * joint->rB.x * joint->rB.x + joint->gamma;

	joint->pivotMass = s2GetInverse22(K);

	joint->C = s2Add(cB, s2Sub(joint->rB, joint->targetA));
	joint->C = s2MulSV(joint->beta, joint->C);

	// Cheat with some damping
	wB *= S2_MAX(0.0f, 1.0f - 0.02f * (60.0f * h));

	vB = s2MulAdd(vB, joint->invMassB, joint->impulse);
	wB += joint->invIB * s2Cross(joint->rB, joint->impulse);

	bodyB->linearVelocity = vB;
	bodyB->angularVelocity = wB;
}

void s2SolveMouse_XPBD(s2Joint* base, s2StepContext* context)
{
	s2MouseJoint* joint = &base->mouseJoint;
	s2Body* bodyB = context->bodies + base->edges[1].bodyIndex;

	s2Vec2 vB = bodyB->linearVelocity;
	float wB = bodyB->angularVelocity;

	// Cdot = v + cross(w, r)
	s2Vec2 Cdot = s2Add(vB, s2CrossSV(wB, joint->rB));
	s2Vec2 SoftCdot = s2Add(Cdot, s2MulAdd(joint->C, joint->gamma, joint->impulse));
	s2Vec2 impulse = s2Neg(s2MulMV(joint->pivotMass, SoftCdot));

	s2Vec2 oldImpulse = joint->impulse;
	joint->impulse = s2Add(joint->impulse, impulse);
	float maxImpulse = context->dt * joint->maxForce;
	if (s2LengthSquared(joint->impulse) > maxImpulse * maxImpulse)
	{
		joint->impulse = s2MulSV(maxImpulse / s2Length(joint->impulse), joint->impulse);
	}
	impulse = s2Sub(joint->impulse, oldImpulse);

	vB = s2MulAdd(vB, joint->invMassB, impulse);
	wB += joint->invIB * s2Cross(joint->rB, impulse);

	bodyB->linearVelocity = vB;
	bodyB->angularVelocity = wB;
}
#endif