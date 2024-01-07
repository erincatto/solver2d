// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "joint.h"

#include "body.h"
#include "contact.h"
#include "core.h"
#include "shape.h"
#include "world.h"

#include "solver2d/debug_draw.h"
#include "solver2d/joint_types.h"

void s2LinearStiffness(float* stiffness, float* damping, float frequencyHertz, float dampingRatio, s2BodyId bodyIdA, s2BodyId bodyIdB)
{
	S2_ASSERT(bodyIdA.world == bodyIdB.world);

	s2World* world = s2GetWorldFromIndex(bodyIdA.world);
	S2_ASSERT(0 <= bodyIdA.index && bodyIdA.index < world->bodyPool.capacity);
	S2_ASSERT(0 <= bodyIdB.index && bodyIdB.index < world->bodyPool.capacity);

	s2Body* bodyA = world->bodies + bodyIdA.index;
	s2Body* bodyB = world->bodies + bodyIdB.index;

	float massA = bodyA->mass;
	float massB = bodyB->mass;
	float mass;
	if (massA > 0.0f && massB > 0.0f)
	{
		mass = massA * massB / (massA + massB);
	}
	else if (massA > 0.0f)
	{
		mass = massA;
	}
	else
	{
		mass = massB;
	}

	float omega = 2.0f * s2_pi * frequencyHertz;
	*stiffness = mass * omega * omega;
	*damping = 2.0f * mass * dampingRatio * omega;
}

void s2AngularStiffness(float* stiffness, float* damping, float frequencyHertz, float dampingRatio, s2BodyId bodyIdA, s2BodyId bodyIdB)
{
	S2_ASSERT(bodyIdA.world == bodyIdB.world);

	s2World* world = s2GetWorldFromIndex(bodyIdA.world);
	S2_ASSERT(0 <= bodyIdA.index && bodyIdA.index < world->bodyPool.capacity);
	S2_ASSERT(0 <= bodyIdB.index && bodyIdB.index < world->bodyPool.capacity);

	s2Body* bodyA = world->bodies + bodyIdA.index;
	s2Body* bodyB = world->bodies + bodyIdB.index;

	float IA = bodyA->I;
	float IB = bodyB->I;
	float I;
	if (IA > 0.0f && IB > 0.0f)
	{
		I = IA * IB / (IA + IB);
	}
	else if (IA > 0.0f)
	{
		I = IA;
	}
	else
	{
		I = IB;
	}

	float omega = 2.0f * s2_pi * frequencyHertz;
	*stiffness = I * omega * omega;
	*damping = 2.0f * I * dampingRatio * omega;
}

static s2Joint* s2CreateJoint(s2World* world, s2Body* bodyA, s2Body* bodyB)
{
	s2Joint* joint = (s2Joint*)s2AllocObject(&world->jointPool);
	world->joints = (s2Joint*)world->jointPool.memory;

	int32_t jointIndex = joint->object.index;

	// Doubly linked list on bodyA
	joint->edges[0].bodyIndex = bodyA->object.index;
	joint->edges[0].prevKey = S2_NULL_INDEX;
	joint->edges[0].nextKey = bodyA->jointList;

	int32_t keyA = (jointIndex << 1) | 0;
	if (bodyA->jointList != S2_NULL_INDEX)
	{
		s2Joint* jointA = world->joints + (bodyA->jointList >> 1);
		s2JointEdge* edgeA = jointA->edges + (bodyA->jointList & 1);
		edgeA->prevKey = keyA;
	}
	bodyA->jointList = keyA;
	bodyA->jointCount += 1;

	// Doubly linked list on bodyB
	joint->edges[1].bodyIndex = bodyB->object.index;
	joint->edges[1].prevKey = S2_NULL_INDEX;
	joint->edges[1].nextKey = bodyB->jointList;

	int32_t keyB = (jointIndex << 1) | 1;
	if (bodyB->jointList != S2_NULL_INDEX)
	{
		s2Joint* jointB = world->joints + (bodyB->jointList >> 1);
		s2JointEdge* edgeB = jointB->edges + (bodyB->jointList & 1);
		edgeB->prevKey = keyB;
	}
	bodyB->jointList = keyB;
	bodyB->jointCount += 1;

	return joint;
}

static void s2DestroyContactsBetweenBodies(s2World* world, s2Body* bodyA, s2Body* bodyB)
{
	int32_t contactKey;
	int32_t otherBodyIndex;

	if (bodyA->contactCount < bodyB->contactCount)
	{
		contactKey = bodyA->contactList;
		otherBodyIndex = bodyB->object.index;
	}
	else
	{
		contactKey = bodyB->contactList;
		otherBodyIndex = bodyA->object.index;
	}

	while (contactKey != S2_NULL_INDEX)
	{
		int32_t contactIndex = contactKey >> 1;
		int32_t edgeIndex = contactKey & 1;

		s2Contact* contact = world->contacts + contactIndex;
		contactKey = contact->edges[edgeIndex].nextKey;

		int32_t otherEdgeIndex = edgeIndex ^ 1;
		if (contact->edges[otherEdgeIndex].bodyIndex == otherBodyIndex)
		{
			// Careful, this removes the contact from the current doubly linked list
			s2DestroyContact(world, contact);
		}
	}
}

s2JointId s2World_CreateMouseJoint(s2WorldId worldId, const s2MouseJointDef* def)
{
	s2World* world = s2GetWorldFromId(worldId);

	S2_ASSERT(s2IsBodyIdValid(world, def->bodyIdA));
	S2_ASSERT(s2IsBodyIdValid(world, def->bodyIdB));

	s2Body* bodyA = world->bodies + def->bodyIdA.index;
	s2Body* bodyB = world->bodies + def->bodyIdB.index;

	s2Joint* joint = s2CreateJoint(world, bodyA, bodyB);

	joint->type = s2_mouseJoint;

	joint->localAnchorA = s2InvTransformPoint(bodyA->transform, def->target);
	joint->localAnchorB = s2InvTransformPoint(bodyB->transform, def->target);

	s2MouseJoint empty = {0};
	joint->mouseJoint = empty;
	joint->mouseJoint.targetA = def->target;
	joint->mouseJoint.maxForce = def->maxForce;
	joint->mouseJoint.stiffness = def->stiffness;
	joint->mouseJoint.damping = def->damping;

	s2JointId jointId = {joint->object.index, world->index, joint->object.revision};

	return jointId;
}

s2JointId s2World_CreateRevoluteJoint(s2WorldId worldId, const s2RevoluteJointDef* def)
{
	s2World* world = s2GetWorldFromId(worldId);

	S2_ASSERT(s2IsBodyIdValid(world, def->bodyIdA));
	S2_ASSERT(s2IsBodyIdValid(world, def->bodyIdB));

	s2Body* bodyA = world->bodies + def->bodyIdA.index;
	s2Body* bodyB = world->bodies + def->bodyIdB.index;

	s2Joint* joint = s2CreateJoint(world, bodyA, bodyB);

	joint->type = s2_revoluteJoint;

	joint->localAnchorA = def->localAnchorA;
	joint->localAnchorB = def->localAnchorB;

	s2RevoluteJoint empty = {0};
	joint->revoluteJoint = empty;

	joint->revoluteJoint.referenceAngle = def->referenceAngle;
	joint->revoluteJoint.impulse = s2Vec2_zero;
	joint->revoluteJoint.axialMass = 0.0f;
	joint->revoluteJoint.motorImpulse = 0.0f;
	joint->revoluteJoint.lowerImpulse = 0.0f;
	joint->revoluteJoint.upperImpulse = 0.0f;
	joint->revoluteJoint.lowerAngle = def->lowerAngle;
	joint->revoluteJoint.upperAngle = def->upperAngle;
	joint->revoluteJoint.maxMotorTorque = def->maxMotorTorque;
	joint->revoluteJoint.motorSpeed = def->motorSpeed;
	joint->revoluteJoint.enableLimit = def->enableLimit;
	joint->revoluteJoint.enableMotor = def->enableMotor;
	joint->revoluteJoint.angle = 0.0f;

	// If the joint prevents collisions, then destroy all contacts between attached bodies
	if (def->collideConnected == false)
	{
		s2DestroyContactsBetweenBodies(world, bodyA, bodyB);
	}

	s2JointId jointId = {joint->object.index, world->index, joint->object.revision};

	return jointId;
}

void s2World_DestroyJoint(s2JointId jointId)
{
	s2World* world = s2GetWorldFromIndex(jointId.world);

	S2_ASSERT(0 <= jointId.index && jointId.index < world->jointPool.capacity);

	s2Joint* joint = world->joints + jointId.index;

	S2_ASSERT(0 <= joint->edges[0].bodyIndex && joint->edges[0].bodyIndex < world->bodyPool.capacity);
	S2_ASSERT(0 <= joint->edges[1].bodyIndex && joint->edges[1].bodyIndex < world->bodyPool.capacity);

	s2JointEdge* edgeA = joint->edges + 0;
	s2JointEdge* edgeB = joint->edges + 1;

	s2Body* bodyA = world->bodies + edgeA->bodyIndex;
	s2Body* bodyB = world->bodies + edgeB->bodyIndex;

	// Remove from body A
	if (edgeA->prevKey != S2_NULL_INDEX)
	{
		s2Joint* prevJoint = world->joints + (edgeA->prevKey >> 1);
		s2JointEdge* prevEdge = prevJoint->edges + (edgeA->prevKey & 1);
		prevEdge->nextKey = edgeA->nextKey;
	}

	if (edgeA->nextKey != S2_NULL_INDEX)
	{
		s2Joint* nextJoint = world->joints + (edgeA->nextKey >> 1);
		s2JointEdge* nextEdge = nextJoint->edges + (edgeA->nextKey & 1);
		nextEdge->prevKey = edgeA->prevKey;
	}

	int32_t edgeKeyA = (joint->object.index << 1) | 0;
	if (bodyA->jointList == edgeKeyA)
	{
		bodyA->jointList = edgeA->nextKey;
	}

	bodyA->jointCount -= 1;

	// Remove from body B
	if (edgeB->prevKey != S2_NULL_INDEX)
	{
		s2Joint* prevJoint = world->joints + (edgeB->prevKey >> 1);
		s2JointEdge* prevEdge = prevJoint->edges + (edgeB->prevKey & 1);
		prevEdge->nextKey = edgeB->nextKey;
	}

	if (edgeB->nextKey != S2_NULL_INDEX)
	{
		s2Joint* nextJoint = world->joints + (edgeB->nextKey >> 1);
		s2JointEdge* nextEdge = nextJoint->edges + (edgeB->nextKey & 1);
		nextEdge->prevKey = edgeB->prevKey;
	}

	int32_t edgeKeyB = (joint->object.index << 1) | 1;
	if (bodyB->jointList == edgeKeyB)
	{
		bodyB->jointList = edgeB->nextKey;
	}

	bodyB->jointCount -= 1;

	s2FreeObject(&world->jointPool, &joint->object);
}

extern void s2InitializeMouse(s2Joint* base, s2StepContext* data);
extern void s2InitializeRevolute(s2Joint* base, s2StepContext* data);

void s2PrepareJoints(s2Joint* joint, s2StepContext* data)
{
	switch (joint->type)
	{
		case s2_mouseJoint:
			s2InitializeMouse(joint, data);
			break;

		case s2_revoluteJoint:
			s2InitializeRevolute(joint, data);
			break;

		default:
			S2_ASSERT(false);
	}
}

extern void s2SolveMouseVelocity(s2Joint* base, s2StepContext* data);
extern void s2SolveRevoluteVelocity(s2Joint* base, s2StepContext* data);

void s2SolveJointVelocity(s2Joint* joint, s2StepContext* data)
{
	switch (joint->type)
	{
		case s2_mouseJoint:
			s2SolveMouseVelocity(joint, data);
			break;

		case s2_revoluteJoint:
			s2SolveRevoluteVelocity(joint, data);
			break;

		default:
			S2_ASSERT(false);
	}
}

extern bool s2SolveRevolutePosition(s2Joint* base, s2StepContext* data);

// This returns true if the position errors are within tolerance.
bool s2SolveJointPosition(s2Joint* joint, s2StepContext* data)
{
	switch (joint->type)
	{
		case s2_revoluteJoint:
			return s2SolveRevolutePosition(joint, data);

		default:
			return true;
	}
}

extern void s2DrawRevolute(s2DebugDraw* draw, s2Joint* base, s2Body* bodyA, s2Body* bodyB);

void s2DrawJoint(s2DebugDraw* draw, s2World* world, s2Joint* joint)
{
	s2Body* bodyA = world->bodies + joint->edges[0].bodyIndex;
	s2Body* bodyB = world->bodies + joint->edges[1].bodyIndex;

	s2Transform xfA = bodyA->transform;
	s2Transform xfB = bodyB->transform;
	s2Vec2 pA = s2TransformPoint(bodyA->transform, joint->localAnchorA);
	s2Vec2 pB = s2TransformPoint(bodyB->transform, joint->localAnchorB);

	s2Color color = {0.5f, 0.8f, 0.8f, 1.0f};

	switch (joint->type)
	{
		case s2_mouseJoint:
		{
			s2Vec2 target = joint->mouseJoint.targetA;

			s2Color c1 = {0.0f, 1.0f, 0.0f, 1.0f};
			draw->DrawPoint(target, 4.0f, c1, draw->context);
			draw->DrawPoint(pB, 4.0f, c1, draw->context);

			s2Color c2 = {0.8f, 0.8f, 0.8f, 1.0f};
			draw->DrawSegment(target, pB, c2, draw->context);
		}
		break;

		case s2_revoluteJoint:
			s2DrawRevolute(draw, joint, bodyA, bodyB);
			break;

		default:
			draw->DrawSegment(xfA.p, pA, color, draw->context);
			draw->DrawSegment(pA, pB, color, draw->context);
			draw->DrawSegment(xfB.p, pB, color, draw->context);
	}
}
