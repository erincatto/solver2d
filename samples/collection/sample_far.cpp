// SPDX-FileCopyrightText: 2022 Erin Catto
// SPDX-License-Identifier: MIT

#include "human.h"
#include "sample.h"
#include "settings.h"

#include "solver2d/geometry.h"
#include "solver2d/math.h"
#include "solver2d/solver2d.h"

#include <math.h>
#include <stdio.h>

// A pyramid far from the origin
class FarPyramid : public Sample
{
public:
	FarPyramid(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		// s2Vec2 origin = {0.0f, 0.0f};
		//s2Vec2 origin = {50000.0f, -45000.0f};
		s2Vec2 origin = {100000.0f, -80000.0f};
		// s2Vec2 origin = {900000.0f, -800000.0f};

		//float originx = 65000.0f;
		//float nextfloat = nextafterf(originx, 2.0f * originx);
		//float ulp = nextfloat - originx;
		//printf("ulp = %g\n", ulp);

		if (settings.restart == false)
		{
			g_camera.m_center = s2Add({0.0f, 6.0f}, origin);
			g_camera.m_zoom = 0.3f;
		}

		{
			s2BodyDef bodyDef = s2_defaultBodyDef;
			bodyDef.position = s2Add({0.0f, -1.0f}, origin);
			s2BodyId groundId = s2CreateBody(m_worldId, &bodyDef);

			s2Polygon box = s2MakeBox(100.0f, 1.0f);
			s2ShapeDef shapeDef = s2_defaultShapeDef;
			s2CreatePolygonShape(groundId, &shapeDef, &box);
		}

		s2BodyDef bodyDef = s2_defaultBodyDef;
		bodyDef.type = s2_dynamicBody;

		s2ShapeDef shapeDef = s2_defaultShapeDef;
		shapeDef.density = 1.0f;

		int baseCount = 10;

		float h = 0.5f;
		s2Polygon box = s2MakeSquare(h);

		float shift = 1.25f * h;

		for (int i = 0; i < baseCount; ++i)
		{
			float y = (2.0f * i + 1.0f) * shift + 0.5f;

			for (int j = i; j < baseCount; ++j)
			{
				float x = (i + 1.0f) * shift + 2.0f * (j - i) * shift - h * baseCount;

				bodyDef.position = s2Add({x, y}, origin);

				s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
				s2CreatePolygonShape(bodyId, &shapeDef, &box);
			}
		}
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new FarPyramid(settings, solverType);
	}
};

static int sampleFarPyramid = RegisterSample("Far", "Pyramid", FarPyramid::Create);

class FarStack : public Sample
{
public:
	FarStack(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		// s2Vec2 origin = {0.0f, 0.0f};
		s2Vec2 origin = {40000.0f, -25000.0f};
		// s2Vec2 origin = {900000.0f, -800000.0f};

		if (settings.restart == false)
		{
			g_camera.m_center = s2Add({0.0f, 1.0f}, origin);
			g_camera.m_zoom = 0.08f;
		}

		{
			s2BodyDef bodyDef = s2_defaultBodyDef;
			bodyDef.position = s2Add({0.0f, -1.0f}, origin);
			s2BodyId groundId = s2CreateBody(m_worldId, &bodyDef);

			s2Polygon box = s2MakeBox(10.0f, 1.0f);
			s2ShapeDef shapeDef = s2_defaultShapeDef;
			s2CreatePolygonShape(groundId, &shapeDef, &box);
		}

		s2BodyDef bodyDef = s2_defaultBodyDef;
		bodyDef.type = s2_dynamicBody;

		s2ShapeDef shapeDef = s2_defaultShapeDef;
		shapeDef.density = 1.0f;

		{
			bodyDef.position = s2Add({1.875f, 0.125f}, origin);
			s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
			s2Circle circle = {{0.0f, 0.0f}, 0.1f};
			s2CreateCircleShape(bodyId, &shapeDef, &circle);
		}

		{
			bodyDef.position = s2Add({-1.875f, 0.15f}, origin);
			s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
			s2Polygon box = s2MakeBox(0.1f, 0.125f);
			s2CreatePolygonShape(bodyId, &shapeDef, &box);
		}

		{
			bodyDef.position = s2Add({0.0f, 0.325f}, origin);
			s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
			s2Polygon box = s2MakeBox(2.0f, 0.05f);
			s2CreatePolygonShape(bodyId, &shapeDef, &box);
		}

		{
			bodyDef.position = s2Add({-0.5f, 0.9f}, origin);
			s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
			s2Polygon box = s2MakeSquare(0.25f);
			s2CreatePolygonShape(bodyId, &shapeDef, &box);
		}

		{
			bodyDef.position = s2Add({-0.55f, 1.7f}, origin);
			s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
			s2Polygon box = s2MakeSquare(0.5f);
			s2CreatePolygonShape(bodyId, &shapeDef, &box);
		}
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new FarStack(settings, solverType);
	}
};

static int sampleFarStack = RegisterSample("Far", "Stack", FarStack::Create);

class FarRecovery : public Sample
{
public:
	FarRecovery(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		// s2Vec2 origin = {0.0f, 0.0f};
		//s2Vec2 origin = {40000.0f, -35000.0f};
		s2Vec2 origin = {80000.0f, -70000.0f};
		//s2Vec2 origin = {900000.0f, -800000.0f};

		if (settings.restart == false)
		{
			g_camera.m_center = s2Add({0.0f, 5.0f}, origin);
			g_camera.m_zoom = 0.25f;
		}

		int baseCount = 4;
		float overlap = 0.25f;
		float extent = 0.5f;

		s2BodyDef bodyDef = s2_defaultBodyDef;
		bodyDef.position = origin;
		s2BodyId groundId = s2CreateBody(m_worldId, &bodyDef);
		s2Segment segment = {{-40.0f, 0.0f}, {40.0f, 0.0f}};
		s2CreateSegmentShape(groundId, &s2_defaultShapeDef, &segment);

		bodyDef.type = s2_dynamicBody;

		s2Polygon box = s2MakeSquare(extent);

		float fraction = 1.0f - overlap;
		float y = extent;
		for (int i = 0; i < baseCount; ++i)
		{
			float x = fraction * extent * (i - baseCount);
			for (int j = i; j < baseCount; ++j)
			{
				bodyDef.position = s2Add({x, y}, origin);
				s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);

				s2CreatePolygonShape(bodyId, &s2_defaultShapeDef, &box);

				x += 2.0f * fraction * extent;
			}

			y += 2.0f * fraction * extent;
		}
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new FarRecovery(settings, solverType);
	}
};

static int sampleFarRecovery = RegisterSample("Far", "Recovery", FarRecovery::Create);

class FarRagdollPile : public Sample
{
public:
	FarRagdollPile(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		//s2Vec2 origin = {0.0f, 0.0f};
		s2Vec2 origin = {10000.0f, -7500.0f};
		//s2Vec2 origin = {20000.0f, -15000.0f};
		//s2Vec2 origin = {80000.0f, -60000.0f};

		if (settings.restart == false)
		{
			g_camera.m_center = s2Add({0.0f, 0.0f}, origin);
			g_camera.m_zoom = 0.055f;
		}

		{
			s2BodyDef bodyDef = s2_defaultBodyDef;
			bodyDef.position = s2Add({0.0f, -1.0f}, origin);
			s2BodyId groundId = s2CreateBody(m_worldId, &bodyDef);
			s2Polygon box;

			box = s2MakeOffsetBox(10.0f, 0.5f, {-5.0f, 2.0f}, -0.15f * s2_pi);
			s2CreatePolygonShape(groundId, &s2_defaultShapeDef, &box);

			box = s2MakeOffsetBox(10.0f, 0.5f, {5.0f, 2.0f}, 0.15f * s2_pi);
			s2CreatePolygonShape(groundId, &s2_defaultShapeDef, &box);
		}

		m_human1.Spawn(m_worldId, s2Add({0.0f, 0.5f}, origin), 1.0f, 1, nullptr);
		m_human2.Spawn(m_worldId, s2Add({-0.2f, 1.0f}, origin), 1.0f, 2, nullptr);
		m_human3.Spawn(m_worldId, s2Add({0.2f, 1.0f}, origin), 1.0f, 3, nullptr);
		m_human4.Spawn(m_worldId, s2Add({-0.4f, 1.5f}, origin), 1.0f, 4, nullptr);
		m_human5.Spawn(m_worldId, s2Add({0.4f, 1.5f}, origin), 1.0f, 5, nullptr);
		m_human6.Spawn(m_worldId, s2Add({0.0f, 2.0f}, origin), 1.0f, 6, nullptr);
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new FarRagdollPile(settings, solverType);
	}

	Human m_human1;
	Human m_human2;
	Human m_human3;
	Human m_human4;
	Human m_human5;
	Human m_human6;
};

static int sampleFarRagdollPile = RegisterSample("Far", "Ragdoll Pile", FarRagdollPile::Create);

class FarChain : public Sample
{
public:
	enum
	{
		e_count = 40
	};

	FarChain(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		// s2Vec2 origin = {0.0f, 0.0f};
		s2Vec2 origin = {40000.0f, -35000.0f};
		// s2Vec2 origin = {900000.0f, -800000.0f};

		if (settings.restart == false)
		{
			g_camera.m_center = s2Add({0.0f, -0.5f}, origin);
			g_camera.m_zoom = 0.2f;
		}

		s2BodyId groundId = s2_nullBodyId;
		{
			s2BodyDef bodyDef = s2_defaultBodyDef;
			bodyDef.position = origin;
			groundId = s2CreateBody(m_worldId, &bodyDef);
		}

		{
			float hx = 0.1f;
			s2Capsule capsule = {{-hx, 0.0f}, {hx, 0.0f}, 0.025f};

			s2ShapeDef shapeDef = s2_defaultShapeDef;
			shapeDef.density = 20.0f;

			s2RevoluteJointDef jointDef = s2DefaultRevoluteJointDef();
			jointDef.drawSize = 0.02f;

			s2Vec2 prevLocalPivot = {0.0f, e_count * hx};
			s2BodyId prevBodyId = groundId;
			for (int i = 0; i < e_count; ++i)
			{
				s2BodyDef bodyDef = s2_defaultBodyDef;
				bodyDef.type = s2_dynamicBody;
				bodyDef.position = s2Add({(1.0f + 2.0f * i) * hx, e_count * hx}, origin);
				bodyDef.linearDamping = 0.1f;
				bodyDef.angularDamping = 0.1f;

				s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
				s2CreateCapsuleShape(bodyId, &shapeDef, &capsule);

				s2Vec2 pivot = {(2.0f * i) * hx, e_count * hx};
				jointDef.bodyIdA = prevBodyId;
				jointDef.bodyIdB = bodyId;
				jointDef.localAnchorA = prevLocalPivot;
				jointDef.localAnchorB = {-hx, 0.0f};
				s2CreateRevoluteJoint(m_worldId, &jointDef);

				prevLocalPivot = {hx, 0.0f};
				prevBodyId = bodyId;
			}
		}
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new FarChain(settings, solverType);
	}
};

static int sampleFarChain = RegisterSample("Joints", "Far Chain", FarChain::Create);
