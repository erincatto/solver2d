// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver2d/id.h"
#include "solver2d/types.h"

typedef struct s2Manifold s2Manifold;

typedef bool s2QueryResultFcn(s2ShapeId shapeId);
typedef float s2RayResultFcn(s2ShapeId shape, s2Vec2 point, s2Vec2 normal, float fraction);
