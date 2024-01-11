// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver2d/types.h"

struct Settings
{
	void Save();
	void Load();

	int m_sampleIndex = 0;
	int m_windowWidth = 1920;
	int m_windowHeight = 1080;
	float m_hertz = 60.0f;
	float m_timeStep = 1.0f / 60.0f;
	int m_velocityIterations = 8;
	int m_positionIterations = 3;
	int m_textLine = 0;
	int m_textIncrement = 18;
	bool m_enablesSolvers[s2_solverTypeCount] = {};
	bool m_enableWarmStarting = true;
	bool m_drawShapes = true;
	bool m_drawJoints = true;
	bool m_drawAABBs = false;
	bool m_drawContactPoints = false;
	bool m_drawContactNormals = false;
	bool m_drawContactImpulse = false;
	bool m_drawFrictionImpulse = false;
	bool m_drawCOMs = false;
	bool m_drawStats = false;
	bool m_pause = false;
	bool m_singleStep = false;
	bool m_restart = false;
};
