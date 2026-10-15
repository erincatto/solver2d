// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

// Scenes draw through these globals. The runner has no GL context so every call is dropped.

#include "draw.h"

Draw g_draw;
Camera g_camera;
GLFWwindow* g_mainWindow = nullptr;

Camera::Camera()
{
	m_center = {0.0f, 0.0f};
	m_zoom = 1.0f;
	m_width = 1280;
	m_height = 800;
}

void Camera::ResetView()
{
}

s2Vec2 Camera::ConvertScreenToWorld(s2Vec2 screenPoint)
{
	return screenPoint;
}

s2Vec2 Camera::ConvertWorldToScreen(s2Vec2 worldPoint)
{
	return worldPoint;
}

void Camera::BuildProjectionMatrix(float*, float)
{
}

Draw::Draw()
{
	m_showUI = false;
	m_points = nullptr;
	m_lines = nullptr;
	m_triangles = nullptr;
	m_roundedTriangles = nullptr;
	m_debugDraw = {};
}

Draw::~Draw()
{
}

void Draw::Create()
{
}

void Draw::Destroy()
{
}

void Draw::DrawPolygon(const s2Vec2*, int, s2Color)
{
}

void Draw::DrawSolidPolygon(const s2Vec2*, int, s2Color)
{
}

void Draw::DrawRoundedPolygon(const s2Vec2*, int, float, s2Color, s2Color)
{
}

void Draw::DrawCircle(s2Vec2, float, s2Color)
{
}

void Draw::DrawSolidCircle(s2Vec2, float, s2Vec2, s2Color)
{
}

void Draw::DrawCapsule(s2Vec2, s2Vec2, float, s2Color)
{
}

void Draw::DrawSolidCapsule(s2Vec2, s2Vec2, float, s2Color)
{
}

void Draw::DrawSegment(s2Vec2, s2Vec2, s2Color)
{
}

void Draw::DrawTransform(s2Transform)
{
}

void Draw::DrawPoint(s2Vec2, float, s2Color)
{
}

void Draw::DrawString(int, int, const char*, ...)
{
}

void Draw::DrawString(s2Vec2, const char*, ...)
{
}

void Draw::DrawAABB(s2Box, s2Color)
{
}

void Draw::Flush()
{
}
