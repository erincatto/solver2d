// SPDX-FileCopyrightText: 2022 Erin Catto
// SPDX-License-Identifier: MIT

#include "sample.h"
#include "settings.h"

#include "solver2d/geometry.h"
#include "solver2d/hull.h"
#include "solver2d/math.h"
#include "solver2d/solver2d.h"

#include <GLFW/glfw3.h>
#include <imgui.h>
#include <math.h>

class SingleBox : public Sample
{
public:
	SingleBox(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		if (settings.restart == false)
		{
			g_camera.m_center = {0.0f, 3.0f};
			g_camera.m_zoom = 0.2f;
		}

		float extent = 1.0f;

		s2BodyDef bodyDef = s2_defaultBodyDef;
		s2BodyId groundId = s2CreateBody(m_worldId, &bodyDef);

		float groundWidth = 66.0f * extent;
		s2ShapeDef shapeDef = s2_defaultShapeDef;
		shapeDef.friction = 0.5f;

		s2Segment segment = {{-0.5f * 2.0f * groundWidth, 0.0f}, {0.5f * 2.0f * groundWidth, 0.0f}};
		s2CreateSegmentShape(groundId, &shapeDef, &segment);
		bodyDef.type = s2_dynamicBody;

		s2Polygon box = s2MakeBox(extent, extent);
		bodyDef.position = {0.0f, 4.0f};
		s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
		s2CreatePolygonShape(bodyId, &shapeDef, &box);
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new SingleBox(settings, solverType);
	}
};

static int sampleSingleBox = RegisterSample("Contact", "Single Box", SingleBox::Create);

class WarmStartEnergy : public Sample
{
public:
	WarmStartEnergy(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		if (settings.restart == false)
		{
			g_camera.m_center = {0.0f, 3.0f};
			g_camera.m_zoom = 0.2f;
		}

		m_topId = s2_nullBodyId;

		s2BodyId groundId = s2CreateBody(m_worldId, &s2_defaultBodyDef);
		s2Segment segment = {{-10.0f, 0.0f}, {10.0f, 0.0f}};
		s2CreateSegmentShape(groundId, &s2_defaultShapeDef, &segment);

		s2BodyDef bodyDef = s2_defaultBodyDef;
		bodyDef.type = s2_dynamicBody;

		s2ShapeDef shapeDef = s2_defaultShapeDef;
		s2Circle circle = {{0.0f, 0.0f}, 0.5f};

		float separation = 0.0f;

		{
			bodyDef.position = {0.0f, 0.5f + separation};
			s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
			shapeDef.density = 1.0f;
			s2CreateCircleShape(bodyId, &shapeDef, &circle);
		}

		{
			bodyDef.position = {0.0f, 1.5f + separation};
			s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
			shapeDef.density = 1.0f;
			s2CreateCircleShape(bodyId, &shapeDef, &circle);
		}

		{
			bodyDef.position = {0.0f, 2.5f + separation};
			m_topId = s2CreateBody(m_worldId, &bodyDef);
			shapeDef.density = 100.0f;
			s2CreateCircleShape(m_topId, &shapeDef, &circle);
		}
	}

	virtual void Step(Settings& settings, s2Color bodyColor) override
	{
		if (m_stepCount == 120 && S2_NON_NULL(m_topId))
		{
			s2DestroyBody(m_topId);
			m_topId = s2_nullBodyId;
		}

		Sample::Step(settings, bodyColor);
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new WarmStartEnergy(settings, solverType);
	}

	s2BodyId m_topId;
};

static int sampleWarmStartEnergy = RegisterSample("Contact", "Warm Start Energy", WarmStartEnergy::Create);

class HighMassRatio1 : public Sample
{
public:
	HighMassRatio1(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		float extent = 1.0f;

		s2BodyDef bodyDef = s2_defaultBodyDef;
		s2BodyId groundId = s2CreateBody(m_worldId, &bodyDef);

		float groundWidth = 66.0f * extent;
		s2ShapeDef shapeDef = s2_defaultShapeDef;
		shapeDef.friction = 0.5f;

		s2Segment segment = {{-0.5f * 2.0f * groundWidth, 0.0f}, {0.5f * 2.0f * groundWidth, 0.0f}};
		s2CreateSegmentShape(groundId, &shapeDef, &segment);

		bodyDef.type = s2_dynamicBody;

		s2Polygon box = s2MakeBox(extent, extent);

#if 0
		//s2Circle circle = {{0.0f, 0.0f}, extent};
		int count = 2;
		for (int i = 0; i < count; ++i)
		{
			bodyDef.position = {0.0f, (2.0f * i + 1.0f) * 1.0f * extent};
			s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);

			shapeDef.density = i == count - 1 ? 300.0f : 1.0f;
			//s2CreateCircleShape(bodyId, &shapeDef, &circle);
			s2CreatePolygonShape(bodyId, &shapeDef, &box);
		}
#else
		for (int j = 0; j < 3; ++j)
		{
			int count = 10;
			float offset = -20.0f * extent + 2.0f * (count + 1.0f) * extent * j;
			float y = extent;
			while (count > 0)
			{
				for (int i = 0; i < count; ++i)
				{
					float coeff = i - 0.5f * count;

					float yy = count == 1 ? y + 2.0f : y;
					bodyDef.position = {2.0f * coeff * extent + offset, yy};
					s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);

					shapeDef.density = count == 1 ? (j + 1.0f) * 100.0f : 1.0f;
					s2CreatePolygonShape(bodyId, &shapeDef, &box);
				}

				--count;
				y += 2.0f * extent;
			}
		}
#endif
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new HighMassRatio1(settings, solverType);
	}
};

static int sampleHighMassRatio1 = RegisterSample("Contact", "High Mass Ratio 1", HighMassRatio1::Create);

// Big box on small boxes
class HighMassRatio2 : public Sample
{
public:
	HighMassRatio2(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		if (settings.restart == false)
		{
			g_camera.m_center = {0.0f, 3.0f};
			g_camera.m_zoom = 1.0f;
		}

		float extent = 1.0f;

		s2BodyDef bodyDef = s2_defaultBodyDef;
		s2BodyId groundId = s2CreateBody(m_worldId, &bodyDef);

		float groundWidth = 66.0f * extent;
		s2ShapeDef shapeDef = s2_defaultShapeDef;
		shapeDef.density = 1.0f;

		s2Segment segment = {{-0.5f * 2.0f * groundWidth, 0.0f}, {0.5f * 2.0f * groundWidth, 0.0f}};
		s2CreateSegmentShape(groundId, &shapeDef, &segment);

		bodyDef.type = s2_dynamicBody;

		s2Vec2 points[3] = {{-0.5f * extent, 0.0f}, {0.5f * extent, 0.0f}, {0.0f, 1.0f * extent}};
		s2Hull hull = s2ComputeHull(points, 3);
		s2Polygon smallTriangle = s2MakePolygon(&hull);
		s2Polygon smallBox = s2MakeBox(0.5f * extent, 0.5f * extent);
		s2Polygon bigBox = s2MakeBox(10.0f * extent, 10.0f * extent);

		{
			bodyDef.position = {-9.0f * extent, 0.5f * extent};
			s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
			s2CreatePolygonShape(bodyId, &shapeDef, &smallBox);
		}

		{
			bodyDef.position = {9.0f * extent, 0.5f * extent};
			s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
			s2CreatePolygonShape(bodyId, &shapeDef, &smallBox);
		}

		{
			bodyDef.position = {0.0f, (10.0f + 16.0f) * extent};
			s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
			s2CreatePolygonShape(bodyId, &shapeDef, &bigBox);
		}
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new HighMassRatio2(settings, solverType);
	}
};

static int sampleHighMassRatio2 = RegisterSample("Contact", "HighMassRatio2", HighMassRatio2::Create);

class Friction : public Sample
{
public:
	Friction(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		if (settings.restart == false)
		{
			g_camera.m_center = {0.0f, 14.0f};
			g_camera.m_zoom = 0.6f;
		}

		{
			s2BodyDef bodyDef = s2_defaultBodyDef;
			s2BodyId groundId = s2CreateBody(m_worldId, &bodyDef);

			s2ShapeDef shapeDef = s2_defaultShapeDef;
			shapeDef.friction = 0.2f;

			s2Segment segment = {{-40.0f, 0.0f}, {40.0f, 0.0f}};
			s2CreateSegmentShape(groundId, &shapeDef, &segment);

			s2Polygon box = s2MakeOffsetBox(13.0f, 0.25f, {-4.0f, 22.0f}, -0.25f);
			s2CreatePolygonShape(groundId, &shapeDef, &box);

			box = s2MakeOffsetBox(0.25f, 1.0f, {10.5f, 19.0f}, 0.0f);
			s2CreatePolygonShape(groundId, &shapeDef, &box);

			box = s2MakeOffsetBox(13.0f, 0.25f, {4.0f, 14.0f}, 0.25f);
			s2CreatePolygonShape(groundId, &shapeDef, &box);

			box = s2MakeOffsetBox(0.25f, 1.0f, {-10.5f, 11.0f}, 0.0f);
			s2CreatePolygonShape(groundId, &shapeDef, &box);

			box = s2MakeOffsetBox(13.0f, 0.25f, {-4.0f, 6.0f}, -0.25f);
			s2CreatePolygonShape(groundId, &shapeDef, &box);
		}

		{
			s2Polygon box = s2MakeBox(0.5f, 0.5f);

			s2ShapeDef shapeDef = s2_defaultShapeDef;
			shapeDef.density = 25.0f;

			float friction[5] = {0.75f, 0.5f, 0.35f, 0.1f, 0.0f};

			for (int i = 0; i < 5; ++i)
			{
				s2BodyDef bodyDef = s2_defaultBodyDef;
				bodyDef.type = s2_dynamicBody;
				bodyDef.position = {-15.0f + 4.0f * i, 28.0f};
				s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);

				shapeDef.friction = friction[i];
				s2CreatePolygonShape(bodyId, &shapeDef, &box);
			}
		}
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new Friction(settings, solverType);
	}
};

static int sampleFriction = RegisterSample("Contact", "Friction", Friction::Create);

class OverlapRecovery : public Sample
{
public:
	OverlapRecovery(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		if (settings.restart == false)
		{
			g_camera.m_zoom = 0.25f;
			g_camera.m_center = {0.0f, 5.0f};
		}

		int baseCount = 4;
		float overlap = 0.25f;
		float extent = 0.5f;

		s2BodyId groundId = s2CreateBody(m_worldId, &s2_defaultBodyDef);
		s2Segment segment = {{-40.0f, 0.0f}, {40.0f, 0.0f}};
		s2CreateSegmentShape(groundId, &s2_defaultShapeDef, &segment);

		s2BodyDef bodyDef = s2_defaultBodyDef;
		bodyDef.type = s2_dynamicBody;

		s2Polygon box = s2MakeSquare(extent);

		float fraction = 1.0f - overlap;
		float y = extent;
		for (int i = 0; i < baseCount; ++i)
		{
			float x = fraction * extent * (i - baseCount);
			for (int j = i; j < baseCount; ++j)
			{
				bodyDef.position = {x, y};
				s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);

				s2CreatePolygonShape(bodyId, &s2_defaultShapeDef, &box);

				x += 2.0f * fraction * extent;
			}

			y += 2.0f * fraction * extent;
		}
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new OverlapRecovery(settings, solverType);
	}
};

static int sampleOverlapRecovery = RegisterSample("Contact", "Overlap Recovery", OverlapRecovery::Create);

class VerticalStack : public Sample
{
public:
	enum ShapeType
	{
		e_circleShape = 0,
		e_boxShape
	};

	VerticalStack(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		if (settings.restart == false)
		{
			g_camera.m_center = {0.0f, 10.0f};
			g_camera.m_zoom = 0.5f;
		}

		{
			s2BodyDef bodyDef = s2_defaultBodyDef;
			bodyDef.position = {0.0f, -1.0f};
			s2BodyId groundId = s2CreateBody(m_worldId, &bodyDef);

			s2Polygon box = s2MakeBox(100.0f, 1.0f);
			s2ShapeDef shapeDef = s2_defaultShapeDef;
			s2CreatePolygonShape(groundId, &shapeDef, &box);
		}

		ShapeType shapeType = e_boxShape;
		int rowCount = 15;

		s2Circle circle = {0};
		circle.radius = 0.5f;

		s2Polygon box = s2MakeBox(0.5f, 0.5f);

		s2ShapeDef shapeDef = s2_defaultShapeDef;
		shapeDef.density = 1.0f;
		shapeDef.friction = 0.3f;

		float offset;

		if (shapeType == e_circleShape)
		{
			offset = 0.0f;
		}
		else
		{
			offset = 0.01f;
		}

		for (int i = 0; i < rowCount; ++i)
		{
			s2BodyDef bodyDef = s2_defaultBodyDef;
			bodyDef.type = s2_dynamicBody;

			float shift = (i % 2 == 0 ? -offset : offset);
			bodyDef.position = {shift, 0.5f + 1.0f * i};
			s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);

			if (shapeType == e_circleShape)
			{
				s2CreateCircleShape(bodyId, &shapeDef, &circle);
			}
			else
			{
				s2CreatePolygonShape(bodyId, &shapeDef, &box);
			}
		}
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new VerticalStack(settings, solverType);
	}
};

static int sampleVerticalStack = RegisterSample("Contact", "Vertical Stack", VerticalStack::Create);

class Pyramid : public Sample
{
public:
	Pyramid(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		if (settings.restart == false)
		{
			g_camera.m_center = {0.0f, 50.0f};
			g_camera.m_zoom = 2.25f;
		}

		{
			s2BodyDef bodyDef = s2_defaultBodyDef;
			bodyDef.position = {0.0f, -1.0f};
			s2BodyId groundId = s2CreateBody(m_worldId, &bodyDef);

			s2Polygon box = s2MakeBox(100.0f, 1.0f);
			s2ShapeDef shapeDef = s2_defaultShapeDef;
			s2CreatePolygonShape(groundId, &shapeDef, &box);
		}

		s2BodyDef bodyDef = s2_defaultBodyDef;
		bodyDef.type = s2_dynamicBody;

		s2ShapeDef shapeDef = s2_defaultShapeDef;
		shapeDef.density = 1.0f;

#ifdef NDEBUG
		int baseCount = 100;
#else
		int baseCount = 20;
#endif

		float h = 0.5f;
		s2Polygon box = s2MakeSquare(h);

		float shift = 1.0f * h;

		for (int i = 0; i < baseCount; ++i)
		{
			float y = (2.0f * i + 1.0f) * shift;

			for (int j = i; j < baseCount; ++j)
			{
				float x = (i + 1.0f) * shift + 2.0f * (j - i) * shift - h * baseCount;

				bodyDef.position = {x, y};

				s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
				s2CreatePolygonShape(bodyId, &shapeDef, &box);
			}
		}
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new Pyramid(settings, solverType);
	}
};

static int samplePyramid = RegisterSample("Contact", "Pyramid", Pyramid::Create);

// A pyramid far from the origin
class FarPyramid : public Sample
{
public:
	FarPyramid(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		// -10k, 5k
		s2Vec2 origin = {-10000.0f, 5000.0f};

		if (settings.restart == false)
		{
			g_camera.m_center = s2Add({0.0f, 50.0f}, origin);
			g_camera.m_zoom = 2.25f;
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

		int baseCount = 40;

		float h = 0.5f;
		s2Polygon box = s2MakeSquare(h);

		float shift = 1.0f * h;

		for (int i = 0; i < baseCount; ++i)
		{
			float y = (2.0f * i + 1.0f) * shift;

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

static int sampleFarPyramid = RegisterSample("Contact", "Far Pyramid", FarPyramid::Create);

// This sample shows an artifact of sub-stepping. The velocity impulse is resolved
// in the first sub-step and then normal impulse drops to zero and so does friction.
// A better result would be had using a force that is applied each substep.
class Rush : public Sample
{
public:
	enum
	{
		e_count = 400
	};

	Rush(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		if (settings.restart == false)
		{
			g_camera.m_center = {0.0f, 0.0f};
			g_camera.m_zoom = 0.5f;
		}

		s2BodyDef bodyDef = s2_defaultBodyDef;
		s2BodyId groundId = s2CreateBody(m_worldId, &bodyDef);

		s2Circle circle = {{0.0f, 0.0f}, 0.5f};

		s2ShapeDef shapeDef = s2_defaultShapeDef;
		shapeDef.friction = 0.2f;
		shapeDef.density = 100.0f;
		s2CreateCircleShape(groundId, &shapeDef, &circle);

		float distance = 5.0f;
		float deltaAngle = 1.0f / distance;
		float deltaDistance = 0.05f;
		float angle = 0.0f;
		bodyDef.type = s2_dynamicBody;
		bodyDef.gravityScale = 0.0f;

		for (int i = 0; i < e_count; ++i)
		{
			bodyDef.position = {distance * cosf(angle), distance * sinf(angle)};
			// bodyDef.linearVelocity = {2.0f * distance * sinf(angle), -1.5f * distance * cosf(angle)};
			m_bodyIds[i] = s2CreateBody(m_worldId, &bodyDef);
			s2CreateCircleShape(m_bodyIds[i], &shapeDef, &circle);

			angle += deltaAngle;
			distance += deltaDistance;
		}
	}

	virtual void Step(Settings& settings, s2Color bodyColor) override
	{
#if 0
		// This approach shows artifacts of sub-stepping
		float speed = 10.0f;
		for (int i = 0; i < e_count; ++i)
		{
			s2Vec2 p = s2Body_GetPosition(m_bodyIds[i]);
			float distance = s2Length(p);
			if (distance < 0.1f)
			{
				continue;
			}

			float scale = speed / distance;
			s2Vec2 v = {-scale * p.x, -scale * p.y};
			s2Body_SetLinearVelocity(m_bodyIds[i], v);
		}
#else
		// forces work better with substepping
		float force = 1000.0f;
		for (int i = 0; i < e_count; ++i)
		{
			s2Vec2 p = s2Body_GetPosition(m_bodyIds[i]);
			float distance = s2Length(p);
			if (distance < 0.1f)
			{
				continue;
			}

			float scale = force / distance;
			s2Vec2 f = {-scale * p.x, -scale * p.y};
			s2Body_ApplyForceToCenter(m_bodyIds[i], f);
		}
#endif

		Sample::Step(settings, bodyColor);
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new Rush(settings, solverType);
	}

	s2BodyId m_bodyIds[e_count];
};

static int sampleRush = RegisterSample("Contact", "Rush", Rush::Create);

class Arch : public Sample
{
public:
	Arch(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		if (settings.restart == false)
		{
			g_camera.m_center = {0.0f, 8.0f};
			g_camera.m_zoom = 0.35f;
		}

		s2Vec2 ps1[9] = {{16.0f, 0.0f},
						 {14.93803712795643f, 5.133601056842984f},
						 {13.79871746027416f, 10.24928069555078f},
						 {12.56252963284711f, 15.34107019122473f},
						 {11.20040987372525f, 20.39856541571217f},
						 {9.66521217819836f, 25.40369899225096f},
						 {7.87179930638133f, 30.3179337000085f},
						 {5.635199558196225f, 35.03820717801641f},
						 {2.405937953536585f, 39.09554102558315f}};

		s2Vec2 ps2[9] = {{24.0f, 0.0f},
						 {22.33619528222415f, 6.02299846205841f},
						 {20.54936888969905f, 12.00964361211476f},
						 {18.60854610798073f, 17.9470321677465f},
						 {16.46769273811807f, 23.81367936585418f},
						 {14.05325025774858f, 29.57079353071012f},
						 {11.23551045834022f, 35.13775818285372f},
						 {7.752568160730571f, 40.30450679009583f},
						 {3.016931552701656f, 44.28891593799322f}};

		float scale = 0.25f;
		for (int i = 0; i < 9; ++i)
		{
			ps1[i] = s2MulSV(scale, ps1[i]);
			ps2[i] = s2MulSV(scale, ps2[i]);
		}

		s2ShapeDef shapeDef = s2_defaultShapeDef;
		shapeDef.friction = 0.6f;

		{
			s2BodyDef bodyDef = s2_defaultBodyDef;
			s2BodyId groundId = s2CreateBody(m_worldId, &bodyDef);
			s2Segment segment = {{-100.0f, 0.0f}, {100.0f, 0.0f}};
			s2CreateSegmentShape(groundId, &shapeDef, &segment);
		}

		s2BodyDef bodyDef = s2_defaultBodyDef;
		bodyDef.type = s2_dynamicBody;

		for (int i = 0; i < 8; ++i)
		{
			s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
			s2Vec2 ps[4] = {ps1[i], ps2[i], ps2[i + 1], ps1[i + 1]};
			s2Hull hull = s2ComputeHull(ps, 4);
			s2Polygon polygon = s2MakePolygon(&hull);
			s2CreatePolygonShape(bodyId, &shapeDef, &polygon);
		}

		for (int i = 0; i < 8; ++i)
		{
			s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
			s2Vec2 ps[4] = {
				{-ps2[i].x, ps2[i].y}, {-ps1[i].x, ps1[i].y}, {-ps1[i + 1].x, ps1[i + 1].y}, {-ps2[i + 1].x, ps2[i + 1].y}};
			s2Hull hull = s2ComputeHull(ps, 4);
			s2Polygon polygon = s2MakePolygon(&hull);
			s2CreatePolygonShape(bodyId, &shapeDef, &polygon);
		}

		{
			s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
			s2Vec2 ps[4] = {ps1[8], ps2[8], {-ps2[8].x, ps2[8].y}, {-ps1[8].x, ps1[8].y}};
			s2Hull hull = s2ComputeHull(ps, 4);
			s2Polygon polygon = s2MakePolygon(&hull);
			s2CreatePolygonShape(bodyId, &shapeDef, &polygon);
		}

		for (int i = 0; i < 4; ++i)
		{
			s2Polygon box = s2MakeBox(2.0f, 0.5f);
			bodyDef.position = {0.0f, 0.5f + ps2[8].y + 1.0f * i};
			s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
			s2CreatePolygonShape(bodyId, &shapeDef, &box);
		}
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new Arch(settings, solverType);
	}
};

static int sampleArch = RegisterSample("Contact", "Arch", Arch::Create);

class DoubleDomino : public Sample
{
public:

	DoubleDomino(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		if (settings.restart == false)
		{
			g_camera.m_center = {0.0f, 4.0f};
			g_camera.m_zoom = 0.25f;
		}

		{
			s2BodyDef bodyDef = s2_defaultBodyDef;
			bodyDef.position = {0.0f, -1.0f};
			s2BodyId groundId = s2CreateBody(m_worldId, &bodyDef);

			s2Polygon box = s2MakeBox(100.0f, 1.0f);
			s2ShapeDef shapeDef = s2_defaultShapeDef;
			s2CreatePolygonShape(groundId, &shapeDef, &box);
		}


		s2Polygon box = s2MakeBox(0.125f, 0.5f);

		s2ShapeDef shapeDef = s2_defaultShapeDef;
		shapeDef.friction = 0.6f;
		s2BodyDef bodyDef = s2_defaultBodyDef;
		bodyDef.type = s2_dynamicBody;

		int count = 15;
		float x = -0.5f * count;
		for (int i = 0; i < count; ++i)
		{
			bodyDef.position = {x, 0.5f};
			s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);
			s2CreatePolygonShape(bodyId, &shapeDef, &box);
			if (i == 0)
			{
				s2Body_ApplyLinearImpulse(bodyId, s2Vec2{0.2f, 0.0f}, s2Vec2{x, 1.0f});
			}

			x += 1.0f;
		}
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new DoubleDomino(settings, solverType);
	}
};

static int sampleDoubleDomino = RegisterSample("Contact", "Double Domino", DoubleDomino::Create);

class Confined : public Sample
{
public:
	enum
	{
		e_gridCount = 25,
		e_maxCount = e_gridCount * e_gridCount
	};

	Confined(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		if (settings.restart == false)
		{
			g_camera.m_center = {0.0f, 10.0f};
			g_camera.m_zoom = 0.5f;
		}

		{
			s2BodyDef bodyDef = s2_defaultBodyDef;
			s2BodyId groundId = s2CreateBody(m_worldId, &bodyDef);

			s2Segment segment;
			segment = {{-10.0f, 0.0f}, {10.0f, 0.0f}};
			s2CreateSegmentShape(groundId, &s2_defaultShapeDef, &segment);
			segment = {{-10.0f, 0.0f}, {-10.0f, 20.0f}};
			s2CreateSegmentShape(groundId, &s2_defaultShapeDef, &segment);
			segment = {{10.0f, 0.0f}, {10.0f, 20.0f}};
			s2CreateSegmentShape(groundId, &s2_defaultShapeDef, &segment);
			segment = {{-10.0f, 20.0f}, {10.0f, 20.0f}};
			s2CreateSegmentShape(groundId, &s2_defaultShapeDef, &segment);
		}

		m_row = 0;
		m_column = 0;
		m_count = 0;
	}

	void CreateCircle()
	{
		float x = -9.0f + m_column * 18.0f / e_gridCount;
		float y = 1.0f + m_row * 18.0f / e_gridCount;

		s2BodyDef bodyDef = s2_defaultBodyDef;
		bodyDef.type = s2_dynamicBody;
		bodyDef.position = {x, y};
		bodyDef.gravityScale = 0.0f;

		s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);

		s2Circle circle = {{0.0f, 0.0f}, 0.5f};
		s2CreateCircleShape(bodyId, &s2_defaultShapeDef, &circle);
	}

	void Step(Settings& settings, s2Color bodyColor) override
	{
		if (m_count < e_maxCount && m_stepCount % 1 == 0)
		{
			CreateCircle();
			m_count += 1;
			m_row += 1;
			if (m_row == e_gridCount)
			{
				m_row = 0;
				m_column += 1;
			}
		}

		Sample::Step(settings, bodyColor);
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new Confined(settings, solverType);
	}

	int m_row;
	int m_column;
	int m_count;
};

static int sampleConfined = RegisterSample("Contact", "Confined", Confined::Create);
