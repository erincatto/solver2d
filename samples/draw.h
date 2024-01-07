// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver2d/debug_draw.h"
#include "solver2d/types.h"

//
struct Camera
{
	Camera();

	void ResetView();
	s2Vec2 ConvertScreenToWorld(s2Vec2 screenPoint);
	s2Vec2 ConvertWorldToScreen(s2Vec2 worldPoint);
	void BuildProjectionMatrix(float* m, float zBias);

	s2Vec2 m_center;
	float m_zoom;
	int m_width;
	int m_height;
};

// This class implements Box2D debug drawing callbacks
class Draw
{
public:
	Draw();
	~Draw();

	void Create();
	void Destroy();

	void DrawPolygon(const s2Vec2* vertices, int vertexCount, s2Color color);

	void DrawSolidPolygon(const s2Vec2* vertices, int vertexCount, s2Color color);

	void DrawRoundedPolygon(const s2Vec2* vertices, int vertexCount, float radius, s2Color fillColor, s2Color outlineColor);

	void DrawCircle(s2Vec2 center, float radius, s2Color color);

	void DrawSolidCircle(s2Vec2 center, float radius, s2Vec2 axis, s2Color color);

	void DrawCapsule(s2Vec2 p1, s2Vec2 p2, float radius, s2Color color);

	void DrawSolidCapsule(s2Vec2 p1, s2Vec2 p2, float radius, s2Color color);

	void DrawSegment(s2Vec2 p1, s2Vec2 p2, s2Color color);

	void DrawTransform(s2Transform xf);

	void DrawPoint(s2Vec2 p, float size, s2Color color);

	void DrawString(int x, int y, const char* string, ...);

	void DrawString(s2Vec2 p, const char* string, ...);

	void DrawAABB(s2Box aabb, s2Color color);

	void Flush();

	bool m_showUI;
	struct GLRenderPoints* m_points;
	struct GLRenderLines* m_lines;
	struct GLRenderTriangles* m_triangles;
	struct GLRenderRoundedTriangles* m_roundedTriangles;
	s2DebugDraw m_debugDraw;
};

extern Draw g_draw;
extern Camera g_camera;
extern struct GLFWwindow* g_mainWindow;
