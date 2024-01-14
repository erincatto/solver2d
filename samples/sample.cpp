// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "sample.h"

#include "settings.h"
#include "solver2d/callbacks.h"
#include "solver2d/manifold.h"
#include "solver2d/math.h"
#include "solver2d/solver2d.h"
#include "solver2d/timer.h"

#include <GLFW/glfw3.h>
#include <stdio.h>
#include <string.h>

Sample::Sample(const Settings& settings, s2SolverType solverType)
{
	s2Vec2 gravity = {0.0f, -10.0f};

	s2WorldDef worldDef = s2DefaultWorldDef();
	worldDef.solverType = solverType;

	m_solverType = solverType;
	m_worldId = s2CreateWorld(&worldDef);
	m_mouseJointId = s2_nullJointId;

	m_stepCount = 0;
}

Sample::~Sample()
{
	s2DestroyWorld(m_worldId);
}

void Sample::DrawTitle(Settings& settings, const char* string)
{
	g_draw.DrawString(5, 5, string);
	settings.textLine = int32_t(26.0f);
}

struct QueryContext
{
	s2Vec2 point;
	s2BodyId bodyId = s2_nullBodyId;
};

bool QueryCallback(s2ShapeId shapeId, void* context)
{
	QueryContext* queryContext = static_cast<QueryContext*>(context);

	s2BodyId bodyId = s2Shape_GetBody(shapeId);
	s2BodyType bodyType = s2Body_GetType(bodyId);
	if (bodyType != s2_dynamicBody)
	{
		// continue query
		return true;
	}

	bool overlap = s2Shape_TestPoint(shapeId, queryContext->point);
	if (overlap)
	{
		// found shape
		queryContext->bodyId = bodyId;
		return false;
	}

	return true;
}

void Sample::MouseDown(s2Vec2 p, int button, int mod)
{
	if (S2_NON_NULL(m_mouseJointId))
	{
		return;
	}

	if (button == GLFW_MOUSE_BUTTON_1)
	{
		// Make a small box.
		s2Box box;
		s2Vec2 d = {0.001f, 0.001f};
		box.lowerBound = s2Sub(p, d);
		box.upperBound = s2Add(p, d);

		// Query the world for overlapping shapes.
		QueryContext queryContext = {p, s2_nullBodyId};
		s2World_QueryAABB(m_worldId, box, QueryCallback, &queryContext);

		if (S2_NON_NULL(queryContext.bodyId))
		{
			float frequencyHertz = 5.0f;
			float dampingRatio = 0.7f;
			float mass = s2Body_GetMass(queryContext.bodyId);

			s2BodyDef bodyDef = s2DefaultBodyDef();
			m_groundBodyId = s2CreateBody(m_worldId, &bodyDef);

			s2MouseJointDef jd;
			jd.bodyIdA = m_groundBodyId;
			jd.bodyIdB = queryContext.bodyId;
			jd.target = p;
			jd.maxForce = 1000.0f * mass;

			float omega = 2.0f * s2_pi * frequencyHertz;
			jd.stiffness = mass * omega * omega;
			jd.damping = 2.0f * mass * dampingRatio * omega;

			m_mouseJointId = s2CreateMouseJoint(m_worldId, &jd);
		}
	}
}

void Sample::MouseUp(s2Vec2 p, int button)
{
	if (S2_NON_NULL(m_mouseJointId) && button == GLFW_MOUSE_BUTTON_1)
	{
		s2DestroyJoint(m_mouseJointId);
		m_mouseJointId = s2_nullJointId;

		s2DestroyBody(m_groundBodyId);
		m_groundBodyId = s2_nullBodyId;
	}
}

void Sample::MouseMove(s2Vec2 p)
{
	if (S2_NON_NULL(m_mouseJointId))
	{
		s2MouseJoint_SetTarget(m_mouseJointId, p);
	}
}

void Sample::Step(Settings& settings, s2Color bodyColor)
{
	bodyColor.a = 0.6f;
	g_draw.m_debugDraw.dynamicBodyColor = bodyColor;

	for (int i = 0; i < settings.multiSteps; ++i)
	{
		s2World_Step(m_worldId, settings.timeStep, settings.velocityIterations, settings.positionIterations, settings.enableWarmStarting);
	}

	if (settings.enablesSolvers[m_solverType])
	{
		s2World_Draw(m_worldId, &g_draw.m_debugDraw);
	}

	if (settings.timeStep > 0.0f)
	{
		++m_stepCount;
	}

	if (settings.drawStats)
	{
		s2Statistics s = s2World_GetStatistics(m_worldId);

		g_draw.DrawString(5, settings.textLine, "bodies/contacts/joints = %d/%d/%d/%d", s.bodyCount, s.contactCount, s.jointCount);
		settings.textLine += settings.textIncrement;

		g_draw.DrawString(5, settings.textLine, "proxies/height = %d/%d", s.proxyCount, s.treeHeight);
		settings.textLine += settings.textIncrement;

		g_draw.DrawString(5, settings.textLine, "stack allocator capacity/used = %d/%d", s.stackCapacity, s.stackUsed);
		settings.textLine += settings.textIncrement;
	}
}

SampleEntry g_sampleEntries[MAX_SAMPLES] = {{nullptr}};
int g_sampleCount = 0;

int RegisterSample(const char* category, const char* name, SampleCreateFcn* fcn)
{
	int index = g_sampleCount;
	if (index < MAX_SAMPLES)
	{
		g_sampleEntries[index] = {category, name, fcn};
		++g_sampleCount;
		return index;
	}

	return -1;
}
