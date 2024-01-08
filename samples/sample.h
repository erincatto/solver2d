// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver2d/id.h"
#include "solver2d/manifold.h"
#include "solver2d/timer.h"
#include "solver2d/types.h"

#include "draw.h"

#include <stdlib.h>

struct Settings;
class Test;

#ifdef _DEBUG
constexpr bool g_sampleDebug = true;
#else
constexpr bool g_sampleDebug = false;
#endif

#define RAND_LIMIT 32767

/// Random number in range [-1,1]
inline float RandomFloat()
{
	float r = (float)(rand() & (RAND_LIMIT));
	r /= RAND_LIMIT;
	r = 2.0f * r - 1.0f;
	return r;
}

/// Random floating point number in range [lo, hi]
inline float RandomFloat(float lo, float hi)
{
	float r = (float)(rand() & (RAND_LIMIT));
	r /= RAND_LIMIT;
	r = (hi - lo) * r + lo;
	return r;
}

class Sample
{
  public:
	Sample(const Settings& settings, s2SolverType solverType);
	virtual ~Sample();

	void DrawTitle(const char* string);
	virtual void Step(Settings& settings, s2Color bodyColor);
	virtual void UpdateUI()
	{
	}
	virtual void Keyboard(int)
	{
	}
	virtual void KeyboardUp(int)
	{
	}
	virtual void MouseDown(s2Vec2 p, int button, int mod);
	virtual void MouseUp(s2Vec2 p, int button);
	virtual void MouseMove(s2Vec2 p);

	s2SolverType m_solverType;
	s2BodyId m_groundBodyId;
	int m_textLine;
	s2WorldId m_worldId;
	s2JointId m_mouseJointId;
	int m_stepCount;
	int m_textIncrement;
};

typedef Sample* SampleCreateFcn(const Settings& settings, s2SolverType solverType);

int RegisterSample(const char* category, const char* name, SampleCreateFcn* fcn);

struct SampleEntry
{
	const char* category;
	const char* name;
	SampleCreateFcn* createFcn;
};

#define MAX_SAMPLES 256
extern SampleEntry g_sampleEntries[MAX_SAMPLES];
extern int g_sampleCount;
