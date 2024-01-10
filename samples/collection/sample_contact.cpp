// SPDX-FileCopyrightText: 2022 Erin Catto
// SPDX-License-Identifier: MIT

#include "sample.h"
#include "settings.h"

#include "solver2d/hull.h"
#include "solver2d/geometry.h"
#include "solver2d/solver2d.h"

#include <GLFW/glfw3.h>
#include <imgui.h>

class SingleBox : public Sample
{
public:
	SingleBox(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		if (settings.m_restart == false)
		{
			g_camera.m_center = {0.0f, 3.0f};
			g_camera.m_zoom = 0.2f;
		}

		float extent = 1.0f;

		s2BodyDef bodyDef = s2DefaultBodyDef();
		s2BodyId groundId = s2CreateBody(m_worldId, &bodyDef);

		float groundWidth = 66.0f * extent;
		s2ShapeDef shapeDef = s2DefaultShapeDef();
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

class HighMassRatio1 : public Sample
{
  public:
	HighMassRatio1(const Settings& settings, s2SolverType solverType)
		  : Sample(settings, solverType)
	{
		float extent = 1.0f;

		s2BodyDef bodyDef = s2DefaultBodyDef();
		s2BodyId groundId = s2CreateBody(m_worldId, &bodyDef);

		float groundWidth = 66.0f * extent;
		s2ShapeDef shapeDef = s2DefaultShapeDef();
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
		if (settings.m_restart == false)
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
		{
			s2BodyDef bodyDef = s2DefaultBodyDef();
			s2BodyId groundId = s2CreateBody(m_worldId, &bodyDef);

			s2ShapeDef shapeDef = s2DefaultShapeDef();
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

			s2ShapeDef shapeDef = s2DefaultShapeDef();
			shapeDef.density = 25.0f;

			float friction[5] = {0.75f, 0.5f, 0.35f, 0.1f, 0.0f};

			for (int i = 0; i < 5; ++i)
			{
				s2BodyDef bodyDef = s2DefaultBodyDef();
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
		if (settings.m_restart == false)
		{
			g_camera.m_zoom = 0.25f;
			g_camera.m_center = {0.0f, 5.0f};
		}

		int baseCount = 4;
		float overlap = 0.5f;
		float extent = 0.3f;
		float pushout = 3.0f;
		float hertz = 30.0f;
		float dampingRatio = 1.0f;

		s2BodyId groundId = s2CreateBody(m_worldId, &s2_defaultBodyDef);

		float groundWidth = 40.0f;

		s2Segment segment = {{-groundWidth, 0.0f}, {groundWidth, 0.0f}};
		s2CreateSegmentShape(groundId, &s2_defaultShapeDef, &segment);

		s2BodyDef bodyDef = s2_defaultBodyDef;
		bodyDef.type = s2_dynamicBody;

		s2Polygon box = s2MakeBox(extent, extent);

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
	enum
	{
		e_maxRows = 30,
	};

	enum ShapeType
	{
		e_circleShape = 0,
		e_boxShape
	};

	VerticalStack(const Settings& settings, s2SolverType solverType)
		: Sample(settings, solverType)
	{
		{
			s2BodyDef bodyDef = s2_defaultBodyDef;
			bodyDef.position = {0.0f, -1.0f};
			s2BodyId groundId = s2CreateBody(m_worldId, &bodyDef);

			s2Polygon box = s2MakeBox(1000.0f, 1.0f);
			s2ShapeDef shapeDef = s2_defaultShapeDef;
			s2CreatePolygonShape(groundId, &shapeDef, &box);
		}

		for (int i = 0; i < e_maxRows; ++i)
		{
			m_bodies[i] = s2_nullBodyId;
		}

		m_shapeType = e_boxShape;
		m_rowCount = 15;

		CreateStack();
	}

	void CreateStack()
	{
		for (int i = 0; i < e_maxRows; ++i)
		{
			if (S2_NON_NULL(m_bodies[i]))
			{
				s2DestroyBody(m_bodies[i]);
				m_bodies[i] = s2_nullBodyId;
			}
		}

		s2Circle circle = {0};
		circle.radius = 0.5f;

		s2Polygon box = s2MakeBox(0.5f, 0.5f);

		s2ShapeDef shapeDef = s2_defaultShapeDef;
		shapeDef.density = 1.0f;
		shapeDef.friction = 0.3f;

		float offset;

		if (m_shapeType == e_circleShape)
		{
			offset = 0.0f;
		}
		else
		{
			offset = 0.01f;
		}

		for (int i = 0; i < m_rowCount; ++i)
		{
			s2BodyDef bodyDef = s2_defaultBodyDef;
			bodyDef.type = s2_dynamicBody;

			float shift = (i % 2 == 0 ? -offset : offset);
			bodyDef.position = {shift, 0.5f + 1.0f * i};
			s2BodyId bodyId = s2CreateBody(m_worldId, &bodyDef);

			m_bodies[i] = bodyId;

			if (m_shapeType == e_circleShape)
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

	s2BodyId m_bodies[e_maxRows];
	int m_rowCount;
	ShapeType m_shapeType;
};

static int sampleVerticalStack = RegisterSample("Contact", "Vertical Stack", VerticalStack::Create);
