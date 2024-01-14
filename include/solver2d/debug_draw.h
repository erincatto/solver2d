// SPDX-FileCopyrightText: 2022 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "types.h"

/// This struct holds callbacks you can implement to draw a solver2d world.
typedef struct s2DebugDraw
{
	/// Draw a closed polygon provided in CCW order.
	void (*DrawPolygon)(const s2Vec2* vertices, int vertexCount, s2Color color, void* context);

	/// Draw a solid closed polygon provided in CCW order.
	void (*DrawSolidPolygon)(const s2Vec2* vertices, int vertexCount, s2Color color, void* context);

	/// Draw a rounded polygon provided in CCW order.
	void (*DrawRoundedPolygon)(const s2Vec2* vertices, int vertexCount, float radius, s2Color lineColor, s2Color fillColor, void* context);

	/// Draw a circle.
	void (*DrawCircle)(s2Vec2 center, float radius, s2Color color, void* context);

	/// Draw a solid circle.
	void (*DrawSolidCircle)(s2Vec2 center, float radius, s2Vec2 axis, s2Color color, void* context);

	/// Draw a capsule.
	void (*DrawCapsule)(s2Vec2 p1, s2Vec2 p2, float radius, s2Color color, void* context);

	/// Draw a solid capsule.
	void (*DrawSolidCapsule)(s2Vec2 p1, s2Vec2 p2, float radius, s2Color color, void* context);

	/// Draw a line segment.
	void (*DrawSegment)(s2Vec2 p1, s2Vec2 p2, s2Color color, void* context);

	/// Draw a transform. Choose your own length scale.
	/// @param xf a transform.
	void (*DrawTransform)(s2Transform xf, void* context);

	/// Draw a point.
	void (*DrawPoint)(s2Vec2 p, float size, s2Color color, void* context);

	/// Draw a string.
	void (*DrawString)(s2Vec2 p, const char* s, void* context);

	s2Color dynamicBodyColor;
	bool drawShapes;
	bool drawJoints;
	bool drawAABBs;
	bool drawMass;
	bool drawContactPoints;
	bool drawContactNormals;
	bool drawContactImpulses;
	bool drawFrictionImpulses;
	void* context;
} s2DebugDraw;
