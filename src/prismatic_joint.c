// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "body.h"
#include "core.h"
#include "joint.h"
#include "solvers.h"
#include "world.h"

#include "solver2d/debug_draw.h"

// Prismatic joint from Box2D without the motor.

// Linear constraint (point-to-line)
// d = pB - pA = xB + rB - xA - rA
// C = dot(perp, d)
// Cdot = dot(d, cross(wA, perp)) + dot(perp, vB + cross(wB, rB) - vA - cross(wA, rA))
//      = -dot(perp, vA) - dot(cross(rA + d, perp), wA) + dot(perp, vB) + dot(cross(rB, perp), vB)
// J = [-perp, -cross(rA + d, perp), perp, cross(rB, perp)]
//
// Angular constraint
// C = aB - aA - referenceAngle
// Cdot = wB - wA
// J = [0 0 -1 0 0 1]
//
// Block solver for the linear and angular constraints
// J = [-perp -s1 perp s2]
//     [0     -1  0    1]
// s1 = cross(rA + d, perp), s2 = cross(rB, perp)
//
// Spring and limit constraints along the axis
// C = dot(axis, d) - target
// J = [-axis -cross(rA + d, axis) axis cross(rB, axis)]

// Geometry at the current body state. Body A's forces act at the point coincident with the anchor on B.
typedef struct s2PrismaticFrame
{
	s2Vec2 rA;
	s2Vec2 rB;
	s2Vec2 d;
	s2Vec2 axis;
	s2Vec2 perp;
	float a1, a2;
	float s1, s2;
} s2PrismaticFrame;

static s2PrismaticFrame s2GetPrismaticFrame(const s2PrismaticJoint* joint, s2Vec2 dcA, s2Rot qA, s2Vec2 dcB, s2Rot qB)
{
	s2PrismaticFrame f;
	f.rA = s2RotateVector(qA, joint->localAnchorA);
	f.rB = s2RotateVector(qB, joint->localAnchorB);
	f.d = s2Add(s2Add(s2Sub(dcB, dcA), s2Sub(f.rB, f.rA)), joint->centerDiff0);
	f.axis = s2RotateVector(qA, joint->localAxisA);
	f.perp = s2LeftPerp(f.axis);

	s2Vec2 rA = s2Add(f.rA, f.d);
	f.a1 = s2Cross(rA, f.axis);
	f.a2 = s2Cross(f.rB, f.axis);
	f.s1 = s2Cross(rA, f.perp);
	f.s2 = s2Cross(f.rB, f.perp);
	return f;
}

static s2PrismaticFrame s2GetBodyFrame(const s2PrismaticJoint* joint, const s2Body* bodyA, const s2Body* bodyB)
{
	return s2GetPrismaticFrame(joint, bodyA->deltaPosition, bodyA->rot, bodyB->deltaPosition, bodyB->rot);
}

static const s2Softness s2_rigidSoftness = {0.0f, 1.0f, 0.0f};

static float s2GetAxialMass(const s2PrismaticJoint* joint, const s2PrismaticFrame* f)
{
	float k = joint->invMassA + joint->invMassB + joint->invIA * f->a1 * f->a1 + joint->invIB * f->a2 * f->a2;
	return k > 0.0f ? 1.0f / k : 0.0f;
}

void s2PreparePrismatic(s2Joint* base, s2StepContext* context, bool warmStart)
{
	S2_ASSERT(base->type == s2_prismaticJoint);

	s2Body* bodyA = context->bodies + base->edges[0].bodyIndex;
	s2Body* bodyB = context->bodies + base->edges[1].bodyIndex;

	s2PrismaticJoint* joint = &base->prismaticJoint;
	joint->localAnchorA = s2Sub(base->localOriginAnchorA, bodyA->localCenter);
	joint->localAnchorB = s2Sub(base->localOriginAnchorB, bodyB->localCenter);
	joint->centerDiff0 = s2Sub(bodyB->position, bodyA->position);
	joint->invMassA = bodyA->invMass;
	joint->invMassB = bodyB->invMass;
	joint->invIA = bodyA->invI;
	joint->invIB = bodyB->invI;

	joint->softness = s2_rigidSoftness;
	joint->limitSoftness = joint->softness;
	joint->springSoftness = s2MakeSoft(joint->hertz, joint->dampingRatio, context->h);

	if (joint->enableSpring == false || warmStart == false)
	{
		joint->springImpulse = 0.0f;
	}

	if (joint->enableLimit == false || warmStart == false)
	{
		joint->lowerImpulse = 0.0f;
		joint->upperImpulse = 0.0f;
	}

	if (warmStart == false)
	{
		joint->impulse = s2Vec2_zero;
	}
}

// Same softness as the revolute joint
void s2PreparePrismatic_Soft(s2Joint* base, s2StepContext* context, float h, float hertz, bool warmStart)
{
	s2PreparePrismatic(base, context, warmStart);

	s2PrismaticJoint* joint = &base->prismaticJoint;
	joint->softness = s2MakeSoft(hertz, 10.0f, h);
	joint->limitSoftness = joint->softness;
}

// A solver that solves the other rows directly can keep the Gauss-Seidel limits softer
void s2PreparePrismaticLimit_Soft(s2Joint* base, float h, float hertz)
{
	base->prismaticJoint.limitSoftness = s2MakeSoft(hertz, 10.0f, h);
}

void s2PreparePrismatic_XPBD(s2Joint* base, s2StepContext* context)
{
	s2PreparePrismatic(base, context, false);
}

void s2WarmStartPrismatic(s2Joint* base, s2StepContext* context)
{
	s2PrismaticJoint* joint = &base->prismaticJoint;

	s2Body* bodyA = context->bodies + base->edges[0].bodyIndex;
	s2Body* bodyB = context->bodies + base->edges[1].bodyIndex;

	s2PrismaticFrame f = s2GetBodyFrame(joint, bodyA, bodyB);

	float axialImpulse = joint->springImpulse + joint->lowerImpulse - joint->upperImpulse;
	float perpImpulse = joint->impulse.x;
	float angleImpulse = joint->impulse.y;

	s2Vec2 P = s2Add(s2MulSV(axialImpulse, f.axis), s2MulSV(perpImpulse, f.perp));
	float LA = axialImpulse * f.a1 + perpImpulse * f.s1 + angleImpulse;
	float LB = axialImpulse * f.a2 + perpImpulse * f.s2 + angleImpulse;

	bodyA->linearVelocity = s2MulSub(bodyA->linearVelocity, joint->invMassA, P);
	bodyA->angularVelocity -= joint->invIA * LA;
	bodyB->linearVelocity = s2MulAdd(bodyB->linearVelocity, joint->invMassB, P);
	bodyB->angularVelocity += joint->invIB * LB;
}

// A real spring, so it applies in relax passes too
static void s2SolvePrismaticSpring(s2PrismaticJoint* joint, s2Body* bodyA, s2Body* bodyB, const s2PrismaticFrame* f)
{
	if (joint->enableSpring == false)
	{
		return;
	}

	s2Vec2 vA = bodyA->linearVelocity;
	float wA = bodyA->angularVelocity;
	s2Vec2 vB = bodyB->linearVelocity;
	float wB = bodyB->angularVelocity;

	float axialMass = s2GetAxialMass(joint, f);
	s2Softness softness = joint->springSoftness;

	float C = s2Dot(f->axis, f->d) - joint->targetTranslation;
	float Cdot = s2Dot(f->axis, s2Sub(vB, vA)) + f->a2 * wB - f->a1 * wA;
	float impulse = -softness.massScale * axialMass * (Cdot + softness.biasRate * C) - softness.impulseScale * joint->springImpulse;
	joint->springImpulse += impulse;

	s2Vec2 P = s2MulSV(impulse, f->axis);
	bodyA->linearVelocity = s2MulSub(vA, joint->invMassA, P);
	bodyA->angularVelocity = wA - joint->invIA * impulse * f->a1;
	bodyB->linearVelocity = s2MulAdd(vB, joint->invMassB, P);
	bodyB->angularVelocity = wB + joint->invIB * impulse * f->a2;
}

// Speculative inside the limits. Past a limit the softness comes from the solver's formulation.
static void s2SolvePrismaticLimits(s2PrismaticJoint* joint, s2Body* bodyA, s2Body* bodyB, const s2PrismaticFrame* f,
								   float inv_h, s2Softness softness)
{
	if (joint->enableLimit == false)
	{
		return;
	}

	s2Vec2 vA = bodyA->linearVelocity;
	float wA = bodyA->angularVelocity;
	s2Vec2 vB = bodyB->linearVelocity;
	float wB = bodyB->angularVelocity;

	float mA = joint->invMassA, mB = joint->invMassB;
	float iA = joint->invIA, iB = joint->invIB;

	float axialMass = s2GetAxialMass(joint, f);
	float translation = s2Dot(f->axis, f->d);

	// Clamp the speculative distance to a reasonable value
	float speculativeDistance = 0.25f * (joint->upperTranslation - joint->lowerTranslation);

	// Lower limit
	{
		float C = translation - joint->lowerTranslation;

		if (C < speculativeDistance)
		{
			float bias = 0.0f;
			float massScale = 1.0f;
			float impulseScale = 0.0f;
			if (C > 0.0f)
			{
				// speculation
				bias = C * inv_h;
			}
			else
			{
				bias = softness.biasRate * C;
				massScale = softness.massScale;
				impulseScale = softness.impulseScale;
			}

			float oldImpulse = joint->lowerImpulse;
			float Cdot = s2Dot(f->axis, s2Sub(vB, vA)) + f->a2 * wB - f->a1 * wA;
			float impulse = -axialMass * massScale * (Cdot + bias) - impulseScale * oldImpulse;
			joint->lowerImpulse = S2_MAX(oldImpulse + impulse, 0.0f);
			impulse = joint->lowerImpulse - oldImpulse;

			s2Vec2 P = s2MulSV(impulse, f->axis);
			vA = s2MulSub(vA, mA, P);
			wA -= iA * impulse * f->a1;
			vB = s2MulAdd(vB, mB, P);
			wB += iB * impulse * f->a2;
		}
		else
		{
			joint->lowerImpulse = 0.0f;
		}
	}

	// Upper limit
	// Note: signs are flipped to keep C positive when the constraint is satisfied.
	// This also keeps the impulse positive when the limit is active.
	{
		float C = joint->upperTranslation - translation;

		if (C < speculativeDistance)
		{
			float bias = 0.0f;
			float massScale = 1.0f;
			float impulseScale = 0.0f;
			if (C > 0.0f)
			{
				// speculation
				bias = C * inv_h;
			}
			else
			{
				bias = softness.biasRate * C;
				massScale = softness.massScale;
				impulseScale = softness.impulseScale;
			}

			float oldImpulse = joint->upperImpulse;
			float Cdot = s2Dot(f->axis, s2Sub(vA, vB)) + f->a1 * wA - f->a2 * wB;
			float impulse = -axialMass * massScale * (Cdot + bias) - impulseScale * oldImpulse;
			joint->upperImpulse = S2_MAX(oldImpulse + impulse, 0.0f);
			impulse = joint->upperImpulse - oldImpulse;

			s2Vec2 P = s2MulSV(impulse, f->axis);
			vA = s2MulAdd(vA, mA, P);
			wA += iA * impulse * f->a1;
			vB = s2MulSub(vB, mB, P);
			wB -= iB * impulse * f->a2;
		}
		else
		{
			joint->upperImpulse = 0.0f;
		}
	}

	bodyA->linearVelocity = vA;
	bodyA->angularVelocity = wA;
	bodyB->linearVelocity = vB;
	bodyB->angularVelocity = wB;
}

// Perpendicular and angle constraints in block form. Rigid with a zero bias rate.
static void s2SolvePrismaticBlock(s2PrismaticJoint* joint, s2Body* bodyA, s2Body* bodyB, const s2PrismaticFrame* f,
								  s2Softness softness)
{
	s2Vec2 vA = bodyA->linearVelocity;
	float wA = bodyA->angularVelocity;
	s2Vec2 vB = bodyB->linearVelocity;
	float wB = bodyB->angularVelocity;

	float mA = joint->invMassA, mB = joint->invMassB;
	float iA = joint->invIA, iB = joint->invIB;

	s2Vec2 Cdot;
	Cdot.x = s2Dot(f->perp, s2Sub(vB, vA)) + f->s2 * wB - f->s1 * wA;
	Cdot.y = wB - wA;

	s2Vec2 C;
	C.x = s2Dot(f->perp, f->d);
	C.y = s2RelativeAngle(bodyB->rot, bodyA->rot) - joint->referenceAngle;

	float k11 = mA + mB + iA * f->s1 * f->s1 + iB * f->s2 * f->s2;
	float k12 = iA * f->s1 + iB * f->s2;
	float k22 = iA + iB;
	if (k22 == 0.0f)
	{
		// For bodies with fixed rotation
		k22 = 1.0f;
	}

	s2Mat22 K = {{k11, k12}, {k12, k22}};
	s2Vec2 b = s2Solve22(K, s2MulAdd(Cdot, softness.biasRate, C));

	s2Vec2 impulse;
	impulse.x = -softness.massScale * b.x - softness.impulseScale * joint->impulse.x;
	impulse.y = -softness.massScale * b.y - softness.impulseScale * joint->impulse.y;
	joint->impulse = s2Add(joint->impulse, impulse);

	s2Vec2 P = s2MulSV(impulse.x, f->perp);
	float LA = impulse.x * f->s1 + impulse.y;
	float LB = impulse.x * f->s2 + impulse.y;

	bodyA->linearVelocity = s2MulSub(vA, mA, P);
	bodyA->angularVelocity = wA - iA * LA;
	bodyB->linearVelocity = s2MulAdd(vB, mB, P);
	bodyB->angularVelocity = wB + iB * LB;
}

// Position error is left to the position solver, as the revolute joint does
void s2SolvePrismatic(s2Joint* base, s2StepContext* context, float h)
{
	s2PrismaticJoint* joint = &base->prismaticJoint;

	s2Body* bodyA = context->bodies + base->edges[0].bodyIndex;
	s2Body* bodyB = context->bodies + base->edges[1].bodyIndex;

	s2PrismaticFrame f = s2GetBodyFrame(joint, bodyA, bodyB);
	float inv_h = h > 0.0f ? 1.0f / h : 0.0f;
	s2SolvePrismaticSpring(joint, bodyA, bodyB, &f);
	s2SolvePrismaticLimits(joint, bodyA, bodyB, &f, inv_h, s2_rigidSoftness);
	s2SolvePrismaticBlock(joint, bodyA, bodyB, &f, s2_rigidSoftness);
}

// Spring and limits only
void s2SolvePrismaticAxial_Soft(s2Joint* base, s2StepContext* context, float inv_h, bool useBias)
{
	s2PrismaticJoint* joint = &base->prismaticJoint;

	s2Body* bodyA = context->bodies + base->edges[0].bodyIndex;
	s2Body* bodyB = context->bodies + base->edges[1].bodyIndex;

	s2PrismaticFrame f = s2GetBodyFrame(joint, bodyA, bodyB);
	s2SolvePrismaticSpring(joint, bodyA, bodyB, &f);
	s2SolvePrismaticLimits(joint, bodyA, bodyB, &f, inv_h, useBias ? joint->limitSoftness : s2_rigidSoftness);
}

void s2SolvePrismatic_Soft(s2Joint* base, s2StepContext* context, float inv_h, bool useBias)
{
	s2SolvePrismaticAxial_Soft(base, context, inv_h, useBias);

	s2PrismaticJoint* joint = &base->prismaticJoint;

	s2Body* bodyA = context->bodies + base->edges[0].bodyIndex;
	s2Body* bodyB = context->bodies + base->edges[1].bodyIndex;

	s2PrismaticFrame f = s2GetBodyFrame(joint, bodyA, bodyB);
	s2SolvePrismaticBlock(joint, bodyA, bodyB, &f, useBias ? joint->softness : s2_rigidSoftness);
}

void s2SolvePrismatic_Baumgarte(s2Joint* base, s2StepContext* context, float inv_h, bool useBias)
{
	s2PrismaticJoint* joint = &base->prismaticJoint;

	s2Body* bodyA = context->bodies + base->edges[0].bodyIndex;
	s2Body* bodyB = context->bodies + base->edges[1].bodyIndex;

	s2PrismaticFrame f = s2GetBodyFrame(joint, bodyA, bodyB);
	s2Softness softness = {useBias ? s2_baumgarte * inv_h : 0.0f, 1.0f, 0.0f};
	s2SolvePrismaticSpring(joint, bodyA, bodyB, &f);
	s2SolvePrismaticLimits(joint, bodyA, bodyB, &f, inv_h, softness);
	s2SolvePrismaticBlock(joint, bodyA, bodyB, &f, softness);
}

// The limit goes first, then the block, as the revolute joint does. The spring is a force and has no position
// correction.
void s2SolvePrismaticPosition(s2Joint* base, s2StepContext* context)
{
	s2PrismaticJoint* joint = &base->prismaticJoint;

	s2Body* bodyA = context->bodies + base->edges[0].bodyIndex;
	s2Body* bodyB = context->bodies + base->edges[1].bodyIndex;

	s2Vec2 dcA = bodyA->deltaPosition;
	s2Rot qA = bodyA->rot;
	s2Vec2 dcB = bodyB->deltaPosition;
	s2Rot qB = bodyB->rot;

	float mA = joint->invMassA, mB = joint->invMassB;
	float iA = joint->invIA, iB = joint->invIB;

	if (joint->enableLimit)
	{
		s2PrismaticFrame f = s2GetPrismaticFrame(joint, dcA, qA, dcB, qB);
		float translation = s2Dot(f.axis, f.d);
		float lower = joint->lowerTranslation;
		float upper = joint->upperTranslation;
		float C = 0.0f;

		if (upper - lower < 2.0f * s2_linearSlop)
		{
			// Prevent large corrections
			C = S2_CLAMP(translation - lower, -s2_maxLinearCorrection, s2_maxLinearCorrection);
		}
		else if (translation <= lower)
		{
			// Prevent large corrections and allow some slop
			C = S2_CLAMP(translation - lower + s2_linearSlop, -s2_maxLinearCorrection, 0.0f);
		}
		else if (translation >= upper)
		{
			// Prevent large corrections and allow some slop
			C = S2_CLAMP(translation - upper - s2_linearSlop, 0.0f, s2_maxLinearCorrection);
		}

		float impulse = -s2GetAxialMass(joint, &f) * C;
		s2Vec2 P = s2MulSV(impulse, f.axis);

		dcA = s2MulSub(dcA, mA, P);
		qA = s2IntegrateRot(qA, -iA * impulse * f.a1);
		dcB = s2MulAdd(dcB, mB, P);
		qB = s2IntegrateRot(qB, iB * impulse * f.a2);
	}

	{
		s2PrismaticFrame f = s2GetPrismaticFrame(joint, dcA, qA, dcB, qB);

		s2Vec2 C;
		C.x = s2Dot(f.perp, f.d);
		C.y = s2RelativeAngle(qB, qA) - joint->referenceAngle;

		float k11 = mA + mB + iA * f.s1 * f.s1 + iB * f.s2 * f.s2;
		float k12 = iA * f.s1 + iB * f.s2;
		float k22 = iA + iB;
		if (k22 == 0.0f)
		{
			// For bodies with fixed rotation
			k22 = 1.0f;
		}

		s2Mat22 K = {{k11, k12}, {k12, k22}};
		s2Vec2 impulse = s2Solve22(K, s2Neg(C));

		s2Vec2 P = s2MulSV(impulse.x, f.perp);
		float LA = impulse.x * f.s1 + impulse.y;
		float LB = impulse.x * f.s2 + impulse.y;

		dcA = s2MulSub(dcA, mA, P);
		qA = s2IntegrateRot(qA, -iA * LA);
		dcB = s2MulAdd(dcB, mB, P);
		qB = s2IntegrateRot(qB, iB * LB);
	}

	bodyA->deltaPosition = dcA;
	bodyA->rot = qA;
	bodyB->deltaPosition = dcB;
	bodyB->rot = qB;
}

// Positional correction at the attachment points, XPBD paper section 3.3.1
static void s2ApplyPositionCorrection_XPBD(s2PrismaticJoint* joint, s2Vec2 rA, s2Vec2 rB, s2Vec2 n, float lambda,
										   s2Vec2* dcA, s2Rot* qA, s2Vec2* dcB, s2Rot* qB)
{
	s2Vec2 p = s2MulSV(lambda, n);
	*dcA = s2MulSub(*dcA, joint->invMassA, p);
	*qA = s2IntegrateRot(*qA, -joint->invIA * s2Cross(rA, p));
	*dcB = s2MulAdd(*dcB, joint->invMassB, p);
	*qB = s2IntegrateRot(*qB, joint->invIB * s2Cross(rB, p));
}

// Generalized inverse mass of a positional correction along the normal, w1 + w2 in the XPBD paper
static float s2GetInverseMass_XPBD(const s2PrismaticJoint* joint, s2Vec2 rA, s2Vec2 rB, s2Vec2 n)
{
	float rnA = s2Cross(rA, n);
	float rnB = s2Cross(rB, n);
	return joint->invMassA + joint->invIA * rnA * rnA + joint->invMassB + joint->invIB * rnB * rnB;
}

void s2SolvePrismatic_XPBD(s2Joint* base, s2StepContext* context, float inv_h)
{
	s2PrismaticJoint* joint = &base->prismaticJoint;

	s2Body* bodyA = context->bodies + base->edges[0].bodyIndex;
	s2Body* bodyB = context->bodies + base->edges[1].bodyIndex;

	float iA = joint->invIA, iB = joint->invIB;

	if (joint->invMassA == 0.0f && joint->invMassB == 0.0f)
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

	// One correction for the offset off the axis and any limit violation along it, applied at the attachment
	// points. Follows the XPBD paper, section 3.4.2.
	{
		s2PrismaticFrame f = s2GetPrismaticFrame(joint, dcA, qA, dcB, qB);

		s2Vec2 separation = s2MulSV(s2Dot(f.perp, f.d), f.perp);
		if (joint->enableLimit)
		{
			float translation = s2Dot(f.axis, f.d);
			if (translation < joint->lowerTranslation)
			{
				separation = s2MulAdd(separation, translation - joint->lowerTranslation, f.axis);
			}
			else if (translation > joint->upperTranslation)
			{
				separation = s2MulAdd(separation, translation - joint->upperTranslation, f.axis);
			}
		}

		float c = s2Length(separation);
		if (c > 0.0f)
		{
			s2Vec2 n = s2MulSV(1.0f / c, separation);
			float lambda = -c / s2GetInverseMass_XPBD(joint, f.rA, f.rB, n);
			s2ApplyPositionCorrection_XPBD(joint, f.rA, f.rB, n, lambda, &dcA, &qA, &dcB, &qB);
		}
	}

	// The spring is a compliant constraint toward the target, XPBD paper section 3.4.2. Compliance from the
	// effective mass keeps the spring hertz meaning what it does in the other solvers. The damping is in the
	// velocity pass.
	if (joint->enableSpring && joint->hertz > 0.0f)
	{
		s2PrismaticFrame f = s2GetPrismaticFrame(joint, dcA, qA, dcB, qB);
		float C = s2Dot(f.axis, f.d) - joint->targetTranslation;
		float w = s2GetInverseMass_XPBD(joint, f.rA, f.rB, f.axis);
		float omega = 2.0f * s2_pi * joint->hertz;
		float compliance = w * inv_h * inv_h / (omega * omega);
		float lambda = -C / (w + compliance);
		s2ApplyPositionCorrection_XPBD(joint, f.rA, f.rB, f.axis, lambda, &dcA, &qA, &dcB, &qB);
	}

	bodyA->deltaPosition = dcA;
	bodyA->rot = qA;
	bodyB->deltaPosition = dcB;
	bodyB->rot = qB;
}

// Spring damping along the axis, XPBD paper eq. 31 and 33. The damping coefficient is the one the damping ratio
// gives the effective mass at the spring hertz.
void s2SolvePrismaticVelocity_XPBD(s2Joint* base, s2StepContext* context, float h)
{
	s2PrismaticJoint* joint = &base->prismaticJoint;
	if (joint->enableSpring == false || joint->hertz == 0.0f || joint->dampingRatio == 0.0f)
	{
		return;
	}

	s2Body* bodyA = context->bodies + base->edges[0].bodyIndex;
	s2Body* bodyB = context->bodies + base->edges[1].bodyIndex;

	s2PrismaticFrame f = s2GetBodyFrame(joint, bodyA, bodyB);
	float w = s2GetInverseMass_XPBD(joint, f.rA, f.rB, f.axis);
	if (w == 0.0f)
	{
		return;
	}

	s2Vec2 vrA = s2Add(bodyA->linearVelocity, s2CrossSV(bodyA->angularVelocity, f.rA));
	s2Vec2 vrB = s2Add(bodyB->linearVelocity, s2CrossSV(bodyB->angularVelocity, f.rB));

	float omega = 2.0f * s2_pi * joint->hertz;
	float damping = 2.0f * joint->dampingRatio * omega;
	float dv = s2Dot(f.axis, s2Sub(vrB, vrA)) * S2_MIN(damping * h, 1.0f);
	s2Vec2 p = s2MulSV(dv / w, f.axis);

	bodyA->linearVelocity = s2MulAdd(bodyA->linearVelocity, joint->invMassA, p);
	bodyA->angularVelocity += joint->invIA * s2Cross(f.rA, p);
	bodyB->linearVelocity = s2MulSub(bodyB->linearVelocity, joint->invMassB, p);
	bodyB->angularVelocity -= joint->invIB * s2Cross(f.rB, p);
}

void s2DrawPrismatic(s2DebugDraw* draw, s2Joint* base, s2Body* bodyA, s2Body* bodyB)
{
	S2_ASSERT(base->type == s2_prismaticJoint);

	s2PrismaticJoint* joint = &base->prismaticJoint;

	s2Vec2 pA = s2TransformPoint(S2_TRANSFORM(bodyA), base->localOriginAnchorA);
	s2Vec2 pB = s2TransformPoint(S2_TRANSFORM(bodyB), base->localOriginAnchorB);
	s2Vec2 axis = s2RotateVector(bodyA->rot, joint->localAxisA);

	s2Color dimGray = {0.4f, 0.4f, 0.4f, 1.0f};
	s2Color gray = {0.7f, 0.7f, 0.7f, 1.0f};
	s2Color green = {0.3f, 0.9f, 0.3f, 1.0f};
	s2Color red = {0.9f, 0.3f, 0.3f, 1.0f};
	s2Color blue = {0.3f, 0.3f, 0.9f, 1.0f};
	s2Color violet = {0.93f, 0.51f, 0.93f, 1.0f};

	draw->DrawSegment(pA, pB, dimGray, draw->context);

	if (joint->enableLimit)
	{
		float b = 0.25f * base->drawSize;
		s2Vec2 lower = s2MulAdd(pA, joint->lowerTranslation, axis);
		s2Vec2 upper = s2MulAdd(pA, joint->upperTranslation, axis);
		s2Vec2 perp = s2LeftPerp(axis);
		draw->DrawSegment(lower, upper, gray, draw->context);
		draw->DrawSegment(s2MulSub(lower, b, perp), s2MulAdd(lower, b, perp), green, draw->context);
		draw->DrawSegment(s2MulSub(upper, b, perp), s2MulAdd(upper, b, perp), red, draw->context);
	}
	else
	{
		draw->DrawSegment(s2MulSub(pA, base->drawSize, axis), s2MulAdd(pA, base->drawSize, axis), gray, draw->context);
	}

	if (joint->enableSpring)
	{
		draw->DrawPoint(s2MulAdd(pA, joint->targetTranslation, axis), 8.0f, violet, draw->context);
	}

	draw->DrawPoint(pA, 5.0f, gray, draw->context);
	draw->DrawPoint(pB, 5.0f, blue, draw->context);
}
