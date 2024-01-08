// SPDX-FileCopyrightText: 2022 Erin Catto
// SPDX-License-Identifier: MIT

#include "sample.h"
#include "settings.h"

#include "solver2d/solver2d.h"
#include "solver2d/geometry.h"

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
		s2BodyId groundId = s2World_CreateBody(m_worldId, &bodyDef);

		float groundWidth = 66.0f * extent;
		s2ShapeDef shapeDef = s2DefaultShapeDef();
		shapeDef.friction = 0.5f;

		s2Segment segment = {{-0.5f * 2.0f * groundWidth, 0.0f}, {0.5f * 2.0f * groundWidth, 0.0f}};
		s2Body_CreateSegment(groundId, &shapeDef, &segment);
		bodyDef.type = s2_dynamicBody;

		s2Polygon box = s2MakeBox(extent, extent);
		bodyDef.position = {0.0f, 4.0f};
		s2BodyId bodyId = s2World_CreateBody(m_worldId, &bodyDef);

		s2Body_CreatePolygon(bodyId, &shapeDef, &box);
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new SingleBox(settings, solverType);
	}
};

static int sampleSingleBox = RegisterSample("Behavior", "Single Box", SingleBox::Create);

class HighMassRatio : public Sample
{
  public:
	HighMassRatio(const Settings& settings, s2SolverType solverType)
		  : Sample(settings, solverType)
	{
		float extent = 1.0f;

		s2BodyDef bodyDef = s2DefaultBodyDef();
		s2BodyId groundId = s2World_CreateBody(m_worldId, &bodyDef);

		float groundWidth = 66.0f * extent;
		s2ShapeDef shapeDef = s2DefaultShapeDef();
		shapeDef.friction = 0.5f;

		s2Segment segment = {{-0.5f * 2.0f * groundWidth, 0.0f}, {0.5f * 2.0f * groundWidth, 0.0f}};
		s2Body_CreateSegment(groundId, &shapeDef, &segment);

		bodyDef.type = s2_dynamicBody;

		s2Polygon box = s2MakeBox(extent, extent);

#if 0
		//s2Circle circle = {{0.0f, 0.0f}, extent};
		int count = 2;
		for (int i = 0; i < count; ++i)
		{
			bodyDef.position = {0.0f, (2.0f * i + 1.0f) * 1.0f * extent};
			s2BodyId bodyId = s2World_CreateBody(m_worldId, &bodyDef);

			shapeDef.density = i == count - 1 ? 300.0f : 1.0f;
			//s2Body_CreateCircle(bodyId, &shapeDef, &circle);
			s2Body_CreatePolygon(bodyId, &shapeDef, &box);
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
					s2BodyId bodyId = s2World_CreateBody(m_worldId, &bodyDef);

					shapeDef.density = count == 1 ? (j + 1.0f) * 100.0f : 1.0f;
					s2Body_CreatePolygon(bodyId, &shapeDef, &box);
				}

				--count;
				y += 2.0f * extent;
			}
		}
#endif
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new HighMassRatio(settings, solverType);
	}
};

static int sampleHighMassRatio = RegisterSample("Behavior", "High Mass Ratio", HighMassRatio::Create);


class Friction : public Sample
{
  public:
	Friction(const Settings& settings, s2SolverType solverType)
		  : Sample(settings, solverType)
	{
		{
			s2BodyDef bodyDef = s2DefaultBodyDef();
			s2BodyId groundId = s2World_CreateBody(m_worldId, &bodyDef);

			s2ShapeDef shapeDef = s2DefaultShapeDef();
			shapeDef.friction = 0.2f;

			s2Segment segment = {{-40.0f, 0.0f}, {40.0f, 0.0f}};
			s2Body_CreateSegment(groundId, &shapeDef, &segment);

			s2Polygon box = s2MakeOffsetBox(13.0f, 0.25f, {-4.0f, 22.0f}, -0.25f);
			s2Body_CreatePolygon(groundId, &shapeDef, &box);

			box = s2MakeOffsetBox(0.25f, 1.0f, {10.5f, 19.0f}, 0.0f);
			s2Body_CreatePolygon(groundId, &shapeDef, &box);

			box = s2MakeOffsetBox(13.0f, 0.25f, {4.0f, 14.0f}, 0.25f);
			s2Body_CreatePolygon(groundId, &shapeDef, &box);

			box = s2MakeOffsetBox(0.25f, 1.0f, {-10.5f, 11.0f}, 0.0f);
			s2Body_CreatePolygon(groundId, &shapeDef, &box);

			box = s2MakeOffsetBox(13.0f, 0.25f, {-4.0f, 6.0f}, -0.25f);
			s2Body_CreatePolygon(groundId, &shapeDef, &box);
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
				s2BodyId bodyId = s2World_CreateBody(m_worldId, &bodyDef);

				shapeDef.friction = friction[i];
				s2Body_CreatePolygon(bodyId, &shapeDef, &box);
			}
		}
	}

	static Sample* Create(const Settings& settings, s2SolverType solverType)
	{
		return new Friction(settings, solverType);
	}
};

static int sampleIndex2 = RegisterSample("Behavior", "Friction", Friction::Create);
