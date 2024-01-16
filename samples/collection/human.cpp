// SPDX-FileCopyrightText: 2023 Erin Catto
// SPDX-License-Identifier: MIT

#include "human.h"

#include "solver2d/solver2d.h"
#include "solver2d/geometry.h"
#include "solver2d/math.h"

#include <assert.h>

Human::Human()
{
	for (int i = 0; i < Bone::e_count; ++i)
	{
		m_bones[i].bodyId = s2_nullBodyId;
		m_bones[i].jointId = s2_nullJointId;
		m_bones[i].parentIndex = -1;
	}

	m_isSpawned = false;
}

void Human::Spawn(s2WorldId worldId, s2Vec2 position, float scale, int groupIndex, void* userData)
{
	assert(m_isSpawned == false);

	s2BodyDef bodyDef = s2_defaultBodyDef;
	bodyDef.type = s2_dynamicBody;
	bodyDef.userData = userData;

	s2ShapeDef shapeDef = s2_defaultShapeDef;
	shapeDef.friction = 0.2f;
	//shapeDef.friction = 0.0f;
	shapeDef.filter.groupIndex = -groupIndex;

	s2ShapeDef footShapeDef = shapeDef;
	footShapeDef.friction = 0.05f;
	//footShapeDef.friction = 0.0f;

	float s = scale;
	float maxTorque = 0.05f * s;
	bool enableMotor = false;
	bool enableLimit = true;
	float drawSize = 0.05f;

	// hip
	{
		Bone* bone = m_bones + Bone::e_hip;
		bone->parentIndex = -1;
		
		bodyDef.position = s2Add({0.0f, 0.95f * s}, position);
		bone->bodyId = s2CreateBody(worldId, &bodyDef);

		s2Capsule capsule = {{0.0f, -0.02f * s}, {0.0f, 0.02f * s}, 0.095f * s};
		s2CreateCapsuleShape(bone->bodyId, &shapeDef, &capsule);
	}

	// torso
	{
		Bone* bone = m_bones + Bone::e_torso;
		bone->parentIndex = Bone::e_hip;
		
		bodyDef.position = s2Add({0.0f, 1.2f * s}, position);
		//bodyDef.type = s2_staticBody;
		bone->bodyId = s2CreateBody(worldId, &bodyDef);
		bodyDef.type = s2_dynamicBody;

		s2Capsule capsule = {{0.0f, -0.135f * s}, {0.0f, 0.135f * s}, 0.09f * s};
		s2CreateCapsuleShape(bone->bodyId, &shapeDef, &capsule);

		s2Vec2 pivot = s2Add({0.0f, 1.0f * s}, position);
		s2RevoluteJointDef jointDef = s2DefaultRevoluteJointDef();
		jointDef.bodyIdA = m_bones[bone->parentIndex].bodyId;
		jointDef.bodyIdB = bone->bodyId;
		jointDef.localAnchorA = s2Body_GetLocalPoint(jointDef.bodyIdA, pivot);
		jointDef.localAnchorB = s2Body_GetLocalPoint(jointDef.bodyIdB, pivot);
		jointDef.enableLimit = enableLimit;
		jointDef.lowerAngle = -0.25f * s2_pi;
		jointDef.upperAngle = 0.0f;
		jointDef.enableMotor = enableMotor;
		jointDef.maxMotorTorque = 0.5f * maxTorque;
		jointDef.drawSize = drawSize;

		bone->jointId = s2CreateRevoluteJoint(worldId, &jointDef);
	}

	// head
	{
		Bone* bone = m_bones + Bone::e_head;
		bone->parentIndex = Bone::e_torso;
		
		bodyDef.position = s2Add({0.0f * s, 1.5f * s}, position);
		bone->bodyId = s2CreateBody(worldId, &bodyDef);

		s2Capsule capsule = {{0.0f, -0.0325f * s}, {0.0f, 0.0325f * s}, 0.08f * s};
		s2CreateCapsuleShape(bone->bodyId, &shapeDef, &capsule);

		// neck
		capsule = {{0.0f, -0.12f * s}, {0.0f, -0.08f * s}, 0.05f * s};
		s2CreateCapsuleShape(bone->bodyId, &shapeDef, &capsule);

		s2Vec2 pivot = s2Add({0.0f, 1.4f * s}, position);
		s2RevoluteJointDef jointDef = s2DefaultRevoluteJointDef();
		jointDef.bodyIdA = m_bones[bone->parentIndex].bodyId;
		jointDef.bodyIdB = bone->bodyId;
		jointDef.localAnchorA = s2Body_GetLocalPoint(jointDef.bodyIdA, pivot);
		jointDef.localAnchorB = s2Body_GetLocalPoint(jointDef.bodyIdB, pivot);
		jointDef.enableLimit = enableLimit;
		jointDef.lowerAngle = -0.3f * s2_pi;
		jointDef.upperAngle = 0.1f * s2_pi;
		jointDef.enableMotor = enableMotor;
		jointDef.maxMotorTorque = 0.25f * maxTorque;
		jointDef.drawSize = drawSize;

		bone->jointId = s2CreateRevoluteJoint(worldId, &jointDef);
	}

	#if 1
	// upper left leg
	{
		Bone* bone = m_bones + Bone::e_upperLeftLeg;
		bone->parentIndex = Bone::e_hip;
		
		bodyDef.position = s2Add({0.0f, 0.775f * s}, position);
		bone->bodyId = s2CreateBody(worldId, &bodyDef);

		s2Capsule capsule = {{0.0f, -0.125f * s}, {0.0f, 0.125f * s}, 0.06f * s};
		s2CreateCapsuleShape(bone->bodyId, &shapeDef, &capsule);

		s2Vec2 pivot = s2Add({0.0f, 0.9f * s}, position);
		s2RevoluteJointDef jointDef = s2DefaultRevoluteJointDef();
		jointDef.bodyIdA = m_bones[bone->parentIndex].bodyId;
		jointDef.bodyIdB = bone->bodyId;
		jointDef.localAnchorA = s2Body_GetLocalPoint(jointDef.bodyIdA, pivot);
		jointDef.localAnchorB = s2Body_GetLocalPoint(jointDef.bodyIdB, pivot);
		jointDef.enableLimit = enableLimit;
		jointDef.lowerAngle = -0.05f * s2_pi;
		jointDef.upperAngle = 0.4f * s2_pi;
		jointDef.enableMotor = enableMotor;
		jointDef.maxMotorTorque = maxTorque;
		jointDef.drawSize = drawSize;

		bone->jointId = s2CreateRevoluteJoint(worldId, &jointDef);
	}

	// lower left leg
	{
		Bone* bone = m_bones + Bone::e_lowerLeftLeg;
		bone->parentIndex = Bone::e_upperLeftLeg;
		
		bodyDef.position = s2Add({0.0f, 0.475f * s}, position);
		bone->bodyId = s2CreateBody(worldId, &bodyDef);

		s2Capsule capsule = {{0.0f, -0.14f * s}, {0.0f, 0.125f * s}, 0.05f * s};
		s2CreateCapsuleShape(bone->bodyId, &shapeDef, &capsule);

		//s2Polygon box = s2MakeOffsetBox(0.1f * s, 0.03f * s, {0.05f * s, -0.175f * s}, 0.0f);
		//s2CreatePolygonShape(bone->bodyId, &shapeDef, &box);

		capsule = {{-0.02f * s, -0.175f * s}, {0.13f * s, -0.175f * s}, 0.03f * s};
		s2CreateCapsuleShape(bone->bodyId, &footShapeDef, &capsule);

		s2Vec2 pivot = s2Add({0.0f, 0.625f * s}, position);
		s2RevoluteJointDef jointDef = s2DefaultRevoluteJointDef();
		jointDef.bodyIdA = m_bones[bone->parentIndex].bodyId;
		jointDef.bodyIdB = bone->bodyId;
		jointDef.localAnchorA = s2Body_GetLocalPoint(jointDef.bodyIdA, pivot);
		jointDef.localAnchorB = s2Body_GetLocalPoint(jointDef.bodyIdB, pivot);
		jointDef.enableLimit = enableLimit;
		jointDef.lowerAngle = -0.5f * s2_pi;
		jointDef.upperAngle = -0.02f * s2_pi;
		jointDef.enableMotor = enableMotor;
		jointDef.maxMotorTorque = 0.5f * maxTorque;
		jointDef.drawSize = drawSize;

		bone->jointId = s2CreateRevoluteJoint(worldId, &jointDef);
	}

	// upper right leg
	{
		Bone* bone = m_bones + Bone::e_upperRightLeg;
		bone->parentIndex = Bone::e_hip;

		bodyDef.position = s2Add({0.0f, 0.775f * s}, position);
		bone->bodyId = s2CreateBody(worldId, &bodyDef);

		s2Capsule capsule = {{0.0f, -0.125f * s}, {0.0f, 0.125f * s}, 0.06f * s};
		s2CreateCapsuleShape(bone->bodyId, &shapeDef, &capsule);

		s2Vec2 pivot = s2Add({0.0f, 0.9f * s}, position);
		s2RevoluteJointDef jointDef = s2DefaultRevoluteJointDef();
		jointDef.bodyIdA = m_bones[bone->parentIndex].bodyId;
		jointDef.bodyIdB = bone->bodyId;
		jointDef.localAnchorA = s2Body_GetLocalPoint(jointDef.bodyIdA, pivot);
		jointDef.localAnchorB = s2Body_GetLocalPoint(jointDef.bodyIdB, pivot);
		jointDef.enableLimit = enableLimit;
		jointDef.lowerAngle = -0.05f * s2_pi;
		jointDef.upperAngle = 0.4f * s2_pi;
		jointDef.enableMotor = enableMotor;
		jointDef.maxMotorTorque = maxTorque;
		jointDef.drawSize = drawSize;

		bone->jointId = s2CreateRevoluteJoint(worldId, &jointDef);
	}

	// lower right leg
	{
		Bone* bone = m_bones + Bone::e_lowerRightLeg;
		bone->parentIndex = Bone::e_upperRightLeg;

		bodyDef.position = s2Add({0.0f, 0.475f * s}, position);
		bone->bodyId = s2CreateBody(worldId, &bodyDef);

		s2Capsule capsule = {{0.0f, -0.14f * s}, {0.0f, 0.125f * s}, 0.05f * s};
		s2CreateCapsuleShape(bone->bodyId, &shapeDef, &capsule);

		//s2Polygon box = s2MakeOffsetBox(0.1f * s, 0.03f * s, {0.05f * s, -0.175f * s}, 0.0f);
		//s2CreatePolygonShape(bone->bodyId, &shapeDef, &box);

		capsule = {{-0.02f * s, -0.175f * s}, {0.13f * s, -0.175f * s}, 0.03f * s};
		s2CreateCapsuleShape(bone->bodyId, &footShapeDef, &capsule);

		s2Vec2 pivot = s2Add({0.0f, 0.625f * s}, position);
		s2RevoluteJointDef jointDef = s2DefaultRevoluteJointDef();
		jointDef.bodyIdA = m_bones[bone->parentIndex].bodyId;
		jointDef.bodyIdB = bone->bodyId;
		jointDef.localAnchorA = s2Body_GetLocalPoint(jointDef.bodyIdA, pivot);
		jointDef.localAnchorB = s2Body_GetLocalPoint(jointDef.bodyIdB, pivot);
		jointDef.enableLimit = enableLimit;
		jointDef.lowerAngle = -0.5f * s2_pi;
		jointDef.upperAngle = -0.02f * s2_pi;
		jointDef.enableMotor = enableMotor;
		jointDef.maxMotorTorque = 0.5f * maxTorque;
		jointDef.drawSize = drawSize;

		bone->jointId = s2CreateRevoluteJoint(worldId, &jointDef);
	}

	// upper left arm
	{
		Bone* bone = m_bones + Bone::e_upperLeftArm;
		bone->parentIndex = Bone::e_torso;

		bodyDef.position = s2Add({0.0f, 1.225f * s}, position);
		bone->bodyId = s2CreateBody(worldId, &bodyDef);

		s2Capsule capsule = {{0.0f, -0.125f * s}, {0.0f, 0.125f * s}, 0.035f * s};
		s2CreateCapsuleShape(bone->bodyId, &shapeDef, &capsule);

		s2Vec2 pivot = s2Add({0.0f, 1.35f * s}, position);
		s2RevoluteJointDef jointDef = s2DefaultRevoluteJointDef();
		jointDef.bodyIdA = m_bones[bone->parentIndex].bodyId;
		jointDef.bodyIdB = bone->bodyId;
		jointDef.localAnchorA = s2Body_GetLocalPoint(jointDef.bodyIdA, pivot);
		jointDef.localAnchorB = s2Body_GetLocalPoint(jointDef.bodyIdB, pivot);
		jointDef.enableLimit = enableLimit;
		jointDef.lowerAngle = -0.05f * s2_pi;
		jointDef.upperAngle = 0.8f * s2_pi;
		jointDef.enableMotor = enableMotor;
		jointDef.maxMotorTorque = 0.25f * maxTorque;
		jointDef.drawSize = drawSize;

		bone->jointId = s2CreateRevoluteJoint(worldId, &jointDef);
	}

	// lower left arm
	{
		Bone* bone = m_bones + Bone::e_lowerLeftArm;
		bone->parentIndex = Bone::e_upperLeftArm;

		bodyDef.position = s2Add({0.0f, 0.975f * s}, position);
		bone->bodyId = s2CreateBody(worldId, &bodyDef);

		s2Capsule capsule = {{0.0f, -0.125f * s}, {0.0f, 0.125f * s}, 0.03f * s};
		s2CreateCapsuleShape(bone->bodyId, &shapeDef, &capsule);

		s2Vec2 pivot = s2Add({0.0f, 1.1f * s}, position);
		s2RevoluteJointDef jointDef = s2DefaultRevoluteJointDef();
		jointDef.bodyIdA = m_bones[bone->parentIndex].bodyId;
		jointDef.bodyIdB = bone->bodyId;
		jointDef.localAnchorA = s2Body_GetLocalPoint(jointDef.bodyIdA, pivot);
		jointDef.localAnchorB = s2Body_GetLocalPoint(jointDef.bodyIdB, pivot);
		jointDef.enableLimit = enableLimit;
		jointDef.lowerAngle = 0.01f * s2_pi;
		jointDef.upperAngle = 0.5f * s2_pi;
		jointDef.enableMotor = enableMotor;
		jointDef.maxMotorTorque = 0.1f * maxTorque;
		jointDef.drawSize = drawSize;

		bone->jointId = s2CreateRevoluteJoint(worldId, &jointDef);
	}

	// upper right arm
	{
		Bone* bone = m_bones + Bone::e_upperRightArm;
		bone->parentIndex = Bone::e_torso;

		bodyDef.position = s2Add({0.0f, 1.225f * s}, position);
		bone->bodyId = s2CreateBody(worldId, &bodyDef);

		s2Capsule capsule = {{0.0f, -0.125f * s}, {0.0f, 0.125f * s}, 0.035f * s};
		s2CreateCapsuleShape(bone->bodyId, &shapeDef, &capsule);

		s2Vec2 pivot = s2Add({0.0f, 1.35f * s}, position);
		s2RevoluteJointDef jointDef = s2DefaultRevoluteJointDef();
		jointDef.bodyIdA = m_bones[bone->parentIndex].bodyId;
		jointDef.bodyIdB = bone->bodyId;
		jointDef.localAnchorA = s2Body_GetLocalPoint(jointDef.bodyIdA, pivot);
		jointDef.localAnchorB = s2Body_GetLocalPoint(jointDef.bodyIdB, pivot);
		jointDef.enableLimit = enableLimit;
		jointDef.lowerAngle = -0.05f * s2_pi;
		jointDef.upperAngle = 0.8f * s2_pi;
		jointDef.enableMotor = enableMotor;
		jointDef.maxMotorTorque = 0.25f * maxTorque;
		jointDef.drawSize = drawSize;

		bone->jointId = s2CreateRevoluteJoint(worldId, &jointDef);
	}

	// lower right arm
	{
		Bone* bone = m_bones + Bone::e_lowerRightArm;
		bone->parentIndex = Bone::e_upperRightArm;

		bodyDef.position = s2Add({0.0f, 0.975f * s}, position);
		bone->bodyId = s2CreateBody(worldId, &bodyDef);

		s2Capsule capsule = {{0.0f, -0.125f * s}, {0.0f, 0.125f * s}, 0.03f * s};
		s2CreateCapsuleShape(bone->bodyId, &shapeDef, &capsule);

		s2Vec2 pivot = s2Add({0.0f, 1.1f * s}, position);
		s2RevoluteJointDef jointDef = s2DefaultRevoluteJointDef();
		jointDef.bodyIdA = m_bones[bone->parentIndex].bodyId;
		jointDef.bodyIdB = bone->bodyId;
		jointDef.localAnchorA = s2Body_GetLocalPoint(jointDef.bodyIdA, pivot);
		jointDef.localAnchorB = s2Body_GetLocalPoint(jointDef.bodyIdB, pivot);
		jointDef.enableLimit = enableLimit;
		jointDef.lowerAngle = 0.01f * s2_pi;
		jointDef.upperAngle = 0.5f * s2_pi;
		jointDef.enableMotor = enableMotor;
		jointDef.maxMotorTorque = 0.1f * maxTorque;
		jointDef.drawSize = drawSize;

		bone->jointId = s2CreateRevoluteJoint(worldId, &jointDef);
	}
	#endif

	m_isSpawned = true;
}

void Human::Despawn()
{
	assert(m_isSpawned == true);

	for (int i = 0; i < Bone::e_count; ++i)
	{
		if (S2_IS_NULL(m_bones[i].jointId))
		{
			continue;
		}

		s2DestroyJoint(m_bones[i].jointId);
		m_bones[i].jointId = s2_nullJointId;
	}

	for (int i = 0; i < Bone::e_count; ++i)
	{
		if (S2_IS_NULL(m_bones[i].bodyId))
		{
			continue;
		}

		s2DestroyBody(m_bones[i].bodyId);
		m_bones[i].bodyId = s2_nullBodyId;
	}

	m_isSpawned = false;
}

s2Vec2 Human::GetBonePosition(int index)
{
	assert(index < Bone::e_count);

	if (S2_IS_NULL(m_bones[index].bodyId))
	{
		return s2Vec2_zero;
	}

	return s2Body_GetPosition(m_bones[index].bodyId);
}
