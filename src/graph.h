// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "array.h"
#include "bitset.h"
#include "table.h"

#include "solver2d/dynamic_tree.h"

typedef struct s2Contact s2Contact;
typedef struct s2StepContext s2StepContext;
typedef struct s2World s2World;

void s2SolveGraphPGS(s2World* world, const s2StepContext* stepContext);
void s2SolveGraphSoftPGS(s2World* world, const s2StepContext* stepContext);
void s2SolveGraphSoftTGS(s2World* world, const s2StepContext* stepContext);
void s2SolveGraphStickyTGS(s2World* world, const s2StepContext* stepContext);
