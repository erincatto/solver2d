// SPDX-FileCopyrightText: 2022 Erin Catto
// SPDX-License-Identifier: MIT

#include "human.h"
#include "sample.h"
#include "settings.h"

#include "solver2d/geometry.h"
#include "solver2d/math.h"
#include "solver2d/solver2d.h"

#include <assert.h>

// A suspension bridge
class Bridge : public Sample
{
public:
	enum
	{
		// e_count = 10
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
			jointDef.drawSize = 0.1f;
			// jointDef.enableMotor = true;
			// jointDef.maxMotorTorque = 200.0f;

			float xbase = -80.0f;

			s2BodyId prevBodyId = groundId;
			for (int i = 0; i < e_count; ++i)
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
				s2CreateRevoluteJoint(m_worldId, &jointDef);

				prevBodyId = bodyId;
			}

			s2Vec2 pivot = {xbase + 1.0f * e_count, 20.0f};
			jointDef.bodyIdA = prevBodyId;
			jointDef.bodyIdB = groundId;
			jointDef.localAnchorA = s2Body_GetLocalPoint(jointDef.bodyIdA, pivot);
			jointDef.localAnchorB = s2Body_GetLocalPoint(jointDef.bodyIdB, pivot);
			// jointDef.enableMotor = true;
			// jointDef.maxMotorTorque = frictionTorque;
			s2CreateRevoluteJoint(m_worldId, &jointDef);
		}
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new Bridge(settings, solverType);
	}
};

static int sampleBridgeIndex = RegisterSample("Joints", "Bridge", Bridge::Create);

#if 1
class BallAndChain : public Sample
{
public:
	enum
	{
		e_count = 40
	};

	BallAndChain(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		if (settings.restart == false)
		{
			g_camera.m_center = {0.0f, -20.0f};
			g_camera.m_zoom = 2.0f;
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
			jointDef.drawSize = 0.1f;

			s2BodyId prevBodyId = groundId;
			for (int i = 0; i < e_count; ++i)
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
				s2CreateRevoluteJoint(m_worldId, &jointDef);

				prevBodyId = bodyId;
			}

			s2Circle circle = {{0.0f, 0.0f}, 8.0f};

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
			s2CreateRevoluteJoint(m_worldId, &jointDef);
		}
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new BallAndChain(settings, solverType);
	}
};

#else

// temp vertical configuration for stability analysis
class BallAndChain : public Sample
{
public:
	enum
	{
		e_count = 2
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
			float hy = 0.5f;
			s2Capsule capsule = {{0.0f, -hy}, {0.0f, hy}, 0.125f};

			s2ShapeDef shapeDef = s2_defaultShapeDef;
			shapeDef.density = 20.0f;

			s2RevoluteJointDef jointDef = s2DefaultRevoluteJointDef();
			jointDef.drawSize = 0.1f;

			s2BodyDef bodyDef = s2_defaultBodyDef;
			bodyDef.type = s2_dynamicBody;

			s2BodyId prevBodyId = groundId;
			for (int i = 0; i < e_count; ++i)
			{
				bodyDef.position = {0.0f, -hy - 2.0f * i * hy};
				s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
				s2CreateCapsuleShape(bodyId, &shapeDef, &capsule);

				s2Vec2 pivot = {0.0f, -2.0f * i * hy};
				jointDef.bodyIdA = prevBodyId;
				jointDef.bodyIdB = bodyId;
				jointDef.localAnchorA = s2Body_GetLocalPoint(jointDef.bodyIdA, pivot);
				jointDef.localAnchorB = s2Body_GetLocalPoint(jointDef.bodyIdB, pivot);
				s2CreateRevoluteJoint(m_worldId, &jointDef);

				prevBodyId = bodyId;
			}

			#if 0
			s2Circle circle = {{0.0f, 0.0f}, 4.0f};
			bodyDef.position = {0.0f, (1.0f - 2.0f * e_count) * hy - circle.radius - hy};
			s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);

			s2CreateCircleShape(bodyId, &shapeDef, &circle);

			s2Vec2 pivot = {0.0f, -2.0f * e_count * hy};
			jointDef.bodyIdA = prevBodyId;
			jointDef.bodyIdB = bodyId;
			jointDef.localAnchorA = s2Body_GetLocalPoint(jointDef.bodyIdA, pivot);
			jointDef.localAnchorB = s2Body_GetLocalPoint(jointDef.bodyIdB, pivot);
			s2CreateRevoluteJoint(m_worldId, &jointDef);
			#endif
		}
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new BallAndChain(settings, solverType);
	}
};
#endif

static int sampleBallAndChainIndex = RegisterSample("Joints", "Ball & Chain", BallAndChain::Create);

class Ragdoll : public Sample
{
public:
	Ragdoll(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		if (settings.restart == false)
		{
			g_camera.m_center = {0.0f, 2.5f};
			g_camera.m_zoom = 0.125f;
		}

		{
			s2BodyDef bodyDef = s2_defaultBodyDef;
			bodyDef.position = {0.0f, -1.0f};
			s2BodyId groundId = s2CreateBody(m_worldId, &bodyDef);
			s2Polygon box = s2MakeBox(20.0f, 1.0f);
			s2CreatePolygonShape(groundId, &s2_defaultShapeDef, &box);
		}

		m_human.Spawn(m_worldId, {0.0f, 4.0f}, 1.0f, 1, nullptr);
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new Ragdoll(settings, solverType);
	}

	Human m_human;
};

static int sampleRagdoll = RegisterSample("Joints", "Ragdoll", Ragdoll::Create);

class RagdollStress : public Sample
{
public:
	enum
	{
		e_count = 32
	};

	RagdollStress(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		if (settings.restart == false)
		{
			g_camera.m_center = {0.0f, 0.0f};
			g_camera.m_zoom = 1.333f;
		}

		{
			s2BodyId groundId = s2CreateBody(m_worldId, &s2_defaultBodyDef);

			s2Vec2 points[] = {
				{-16.8672504, 31.088623},	 {16.8672485, 31.088623},	 {16.8672485, 17.1978741}, {8.26824951, 11.906374},
				{16.8672485, 11.906374},	 {16.8672485, -0.661376953}, {8.26824951, -5.953125},  {16.8672485, -5.953125},
				{16.8672485, -13.229126},	 {3.63799858, -23.151123},	 {3.63799858, -31.088623}, {-3.63800049, -31.088623},
				{-3.63800049, -23.151123},	 {-16.8672504, -13.229126},	 {-16.8672504, -5.953125}, {-8.26825142, -5.953125},
				{-16.8672504, -0.661376953}, {-16.8672504, 11.906374},	 {-8.26825142, 11.906374}, {-16.8672504, 17.1978741},
			};

			int count = sizeof(points) / sizeof(points[0]);

			s2ShapeDef shapeDef = s2_defaultShapeDef;
			shapeDef.friction = 0.2f;
			for (int i = 0; i < count; ++i)
			{
				int i1 = i;
				int i2 = (i + 1) % count;
				s2Capsule capsule = {points[i1], points[i2], 0.5f};
				s2CreateCapsuleShape(groundId, &shapeDef, &capsule);
			}

			float sign = 1.0f;
			float y = 14.0f;
			for (int i = 0; i < 3; ++i)
			{
				s2BodyDef bodyDef = s2_defaultBodyDef;
				bodyDef.position = {0.0f, y};
				bodyDef.type = s2_dynamicBody;

				s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);

				s2Polygon box = s2MakeBox(6.0f, 0.5f);
				s2ShapeDef shapeDef = s2_defaultShapeDef;
				shapeDef.friction = 0.1f;
				shapeDef.restitution = 1.0f;
				shapeDef.density = 1.0f;

				s2CreatePolygonShape(bodyId, &shapeDef, &box);

				s2RevoluteJointDef revoluteDef = s2DefaultRevoluteJointDef();
				revoluteDef.bodyIdA = groundId;
				revoluteDef.bodyIdB = bodyId;
				revoluteDef.localAnchorA = bodyDef.position;
				revoluteDef.localAnchorB = s2Vec2_zero;
				revoluteDef.maxMotorTorque = 200.0f;
				revoluteDef.motorSpeed = 5.0f * sign;
				revoluteDef.enableMotor = true;
				revoluteDef.drawSize = 0.2f;

				s2CreateRevoluteJoint(m_worldId, &revoluteDef);

				y -= 14.0f;
				sign = -sign;
			}
		}

		m_wait = 0.5f;
		m_side = -15.0f;

		for (int i = 0; i < e_count; ++i)
		{
			m_isSpawned[i] = false;
		}

		CreateElement();
	}

	void CreateElement()
	{
		int index = -1;
		for (int i = 0; i < e_count; ++i)
		{
			if (m_isSpawned[i] == false)
			{
				index = i;
				break;
			}
		}

		if (index == -1)
		{
			return;
		}

		s2Vec2 center = {m_side, 28.0f};

		Human* human = m_humans + index;
		human->Spawn(m_worldId, center, 2.0f, index + 1, human);

		m_isSpawned[index] = true;
		m_side = -m_side;
	}

	void Step(Settings& settings, s2Color bodyColor) override
	{
		Sample::Step(settings, bodyColor);

		for (int i = 0; i < e_count; ++i)
		{
			if (m_isSpawned[i] == false)
			{
				continue;
			}

			s2Vec2 p = m_humans[i].GetBonePosition(Bone::e_torso);

			if (p.y < -25.0f)
			{
				m_humans[i].Despawn();
				m_isSpawned[i] = false;
			}
		}

		if (settings.hertz > 0.0f && settings.pause == false)
		{
			m_wait -= 1.0f / settings.hertz;
			if (m_wait < 0.0f)
			{
				CreateElement();
				m_wait += 0.5f;
			}
		}
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new RagdollStress(settings, solverType);
	}

	Human m_humans[e_count];
	bool m_isSpawned[e_count];
	float m_wait;
	float m_side;
};

static int sampleRagdollStress = RegisterSample("Joints", "Ragdoll Stress", RagdollStress::Create);

// This stresses the joint solver. It may blow up and/or crash for some solvers.
class JointGrid : public Sample
{
public:
	JointGrid(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		if (settings.restart == false)
		{
			g_camera.m_center = {55.0f, -55.0f};
			g_camera.m_zoom = 2.5f;
		}

		float rad = 0.4f;
#ifdef NDEBUG
		int numi = 100;
		int numk = 100;
#else
		int numi = 1;
		int numk = 20;
#endif

		float shift = 1.0f;

		// Allocate to avoid huge stack usage
		s2BodyId* bodies = static_cast<s2BodyId*>(malloc(numi * numk * sizeof(s2BodyId)));
		int index = 0;

		float density = 1.0f;
		s2ShapeDef shapeDef = s2_defaultShapeDef;
		shapeDef.filter.categoryBits = 2;
		shapeDef.filter.maskBits = ~2u;

		s2Circle circle = {0};
		circle.radius = rad;

		s2RevoluteJointDef jointDef = s2DefaultRevoluteJointDef();
		jointDef.drawSize = 0.2f;

		for (int k = 0; k < numk; ++k)
		{
			for (int i = 0; i < numi; ++i)
			{
				s2BodyDef bodyDef = s2_defaultBodyDef;
				//if (k == 0 && i == 0)
				//if ((k == 0 || k == numk - 1) && i == 0)
				if (k >= numk / 2 - 3 && k <= numk / 2 + 3 && i == 0)
				{
					bodyDef.type = s2_staticBody;
				}
				else
				{
					bodyDef.type = s2_dynamicBody;
				}

				bodyDef.position = {k * shift, -i * shift};
				bodyDef.gravityScale = 2.0f;

				s2BodyId body = s2CreateBody(m_worldId, &bodyDef);

				s2CreateCircleShape(body, &shapeDef, &circle);

				if (i > 0)
				{
					jointDef.bodyIdA = bodies[index - 1];
					jointDef.bodyIdB = body;
					jointDef.localAnchorA = {0.0f, -0.5f * shift};
					jointDef.localAnchorB = {0.0f, 0.5f * shift};
					s2CreateRevoluteJoint(m_worldId, &jointDef);
				}

				if (k > 0)
				{
					jointDef.bodyIdA = bodies[index - numi];
					jointDef.bodyIdB = body;
					jointDef.localAnchorA = {0.5f * shift, 0.0f};
					jointDef.localAnchorB = {-0.5f * shift, 0.0f};
					s2CreateRevoluteJoint(m_worldId, &jointDef);
				}

				bodies[index++] = body;
			}
		}

		free(bodies);
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new JointGrid(settings, solverType);
	}
};

static int sampleJointGrid = RegisterSample("Joints", "Joint Grid", JointGrid::Create);
