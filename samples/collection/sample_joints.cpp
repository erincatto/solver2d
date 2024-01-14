// SPDX-FileCopyrightText: 2022 Erin Catto
// SPDX-License-Identifier: MIT

//#include "human.h"
#include "sample.h"
#include "settings.h"

#include "solver2d/solver2d.h"
#include "solver2d/geometry.h"
#include "solver2d/math.h"

#include <assert.h>

// A suspension bridge
class Bridge : public Sample
{
public:
	enum
	{
		//e_count = 10
		e_count = 160
	};

	Bridge(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		if (settings.restart == false)
		{
			g_camera.m_zoom = 2.5f;
		}

		s2BodyId groundId = s2_nullBodyId;
		{
			s2BodyDef bodyDef = s2_defaultBodyDef;
			groundId = s2CreateBody(m_worldId, &bodyDef);
		}

		{
			s2Polygon box = s2MakeBox(0.5f, 0.125f);

			s2ShapeDef shapeDef = s2_defaultShapeDef;
			shapeDef.density = 20.0f;

			s2RevoluteJointDef jointDef = s2DefaultRevoluteJointDef();
			int32_t jointIndex = 0;
			m_frictionTorque = 200.0f;

			float xbase = -80.0f;

			s2BodyId prevBodyId = groundId;
			for (int32_t i = 0; i < e_count; ++i)
			{
				s2BodyDef bodyDef = s2_defaultBodyDef;
				bodyDef.type = s2_dynamicBody;
				bodyDef.position = {xbase + 0.5f + 1.0f * i, 20.0f};
				bodyDef.linearDamping = 0.1f;
				bodyDef.angularDamping = 0.1f;

				s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
				s2CreatePolygonShape(bodyId, &shapeDef, &box);

				s2Vec2 pivot = {xbase + 1.0f * i, 20.0f};
				jointDef.bodyIdA = prevBodyId;
				jointDef.bodyIdB = bodyId;
				jointDef.localAnchorA = s2Body_GetLocalPoint(jointDef.bodyIdA, pivot);
				jointDef.localAnchorB = s2Body_GetLocalPoint(jointDef.bodyIdB, pivot);
				//jointDef.enableMotor = true;
				//jointDef.maxMotorTorque = m_frictionTorque;
				s2CreateRevoluteJoint(m_worldId, &jointDef);

				prevBodyId = bodyId;
			}

			s2Vec2 pivot = {xbase + 1.0f * e_count, 20.0f};
			jointDef.bodyIdA = prevBodyId;
			jointDef.bodyIdB = groundId;
			jointDef.localAnchorA = s2Body_GetLocalPoint(jointDef.bodyIdA, pivot);
			jointDef.localAnchorB = s2Body_GetLocalPoint(jointDef.bodyIdB, pivot);
			//jointDef.enableMotor = true;
			//jointDef.maxMotorTorque = m_frictionTorque;
			s2CreateRevoluteJoint(m_worldId, &jointDef);
		}
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new Bridge(settings, solverType);
	}

	float m_frictionTorque;
};

static int sampleBridgeIndex = RegisterSample("Joints", "Bridge", Bridge::Create);

class BallAndChain : public Sample
{
public:
	enum
	{
		e_count = 30
	};

	BallAndChain(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		if (settings.restart == false)
		{
			g_camera.m_center = {0.0f, -5.0f};
		}

		s2BodyId groundId = s2_nullBodyId;
		{
			s2BodyDef bodyDef = s2_defaultBodyDef;
			groundId = s2CreateBody(m_worldId, &bodyDef);
		}

		{
			float hx = 0.5f;
			s2Capsule capsule = {{-hx, 0.0f}, {hx, 0.0f}, 0.125f};

			s2ShapeDef shapeDef = s2_defaultShapeDef;
			shapeDef.density = 20.0f;

			s2RevoluteJointDef jointDef = s2DefaultRevoluteJointDef();

			s2BodyId prevBodyId = groundId;
			for (int32_t i = 0; i < e_count; ++i)
			{
				s2BodyDef bodyDef = s2_defaultBodyDef;
				bodyDef.type = s2_dynamicBody;
				bodyDef.position = {(1.0f + 2.0f * i) * hx, e_count * hx};
				bodyDef.linearDamping = 0.1f;
				bodyDef.angularDamping = 0.1f;

				s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
				s2CreateCapsuleShape(bodyId, &shapeDef, &capsule);

				s2Vec2 pivot = {(2.0f * i) * hx, e_count * hx};
				jointDef.bodyIdA = prevBodyId;
				jointDef.bodyIdB = bodyId;
				jointDef.localAnchorA = s2Body_GetLocalPoint(jointDef.bodyIdA, pivot);
				jointDef.localAnchorB = s2Body_GetLocalPoint(jointDef.bodyIdB, pivot);
				// jointDef.enableMotor = true;
				//jointDef.maxMotorTorque = 100.0f;
				s2CreateRevoluteJoint(m_worldId, &jointDef);

				prevBodyId = bodyId;
			}

			#if 1
			s2Circle circle = {{0.0f, 0.0f}, 4.0f};

			s2BodyDef bodyDef = s2_defaultBodyDef;
			bodyDef.type = s2_dynamicBody;
			bodyDef.position = {(1.0f + 2.0f * e_count) * hx + circle.radius - hx, e_count * hx};
			bodyDef.linearDamping = 0.1f;
			bodyDef.angularDamping = 0.1f;

			s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
			s2CreateCircleShape(bodyId, &shapeDef, &circle);

			s2Vec2 pivot = {(2.0f * e_count) * hx, e_count * hx};
			jointDef.bodyIdA = prevBodyId;
			jointDef.bodyIdB = bodyId;
			jointDef.localAnchorA = s2Body_GetLocalPoint(jointDef.bodyIdA, pivot);
			jointDef.localAnchorB = s2Body_GetLocalPoint(jointDef.bodyIdB, pivot);
			//jointDef.enableMotor = true;
			//jointDef.maxMotorTorque = m_frictionTorque;
			s2CreateRevoluteJoint(m_worldId, &jointDef);
			#endif
		}
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new BallAndChain(settings, solverType);
	}
};

static int sampleBallAndChainIndex = RegisterSample("Joints", "Ball & Chain", BallAndChain::Create);

#if 0
class Ragdoll : public Sample
{
public:
	Ragdoll(const Settings& settings)
		: Sample(settings)
	{
		if (settings.restart == false)
		{
			g_camera.m_zoom = 0.25f;
			g_camera.m_center = {0.0f, 5.0f};
		}

		s2BodyId groundId;
		{
			groundId = s2CreateBody(m_worldId, &s2_defaultBodyDef);
			s2Segment segment = {{-20.0f, 0.0f}, {20.0f, 0.0f}};
			s2CreateSegmentShape(groundId, &s2_defaultShapeDef, &segment);
		}

		m_human.Spawn(m_worldId, {0.0f, 10.0f}, 1.0f, 1, nullptr);
	}

	static Sample* Create(const Settings& settings)
	{
		return new Ragdoll(settings);
	}

	Human m_human;
};

static int sampleRagdoll = RegisterSample("Joints", "Ragdoll", Ragdoll::Create);
#endif
