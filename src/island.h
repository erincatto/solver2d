// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include <stdint.h>

typedef struct b2StepContext b2StepContext;
typedef struct b2World b2World;

void b2SolveWorld(b2World* world, b2StepContext* context);
