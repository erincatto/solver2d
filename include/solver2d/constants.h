// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#define s2_pi 3.14159265359f
#define s2_speculativeDistance (16.0f * s2_linearSlop)
#define s2_aabbMargin (0.1f + s2_speculativeDistance)
#define s2_linearSlop 0.005f
#define s2_angularSlop (2.0f / 180.0f * s2_pi)
#define s2_maxPolygonVertices 8
#define s2_maxWorlds 32
#define s2_maxLinearCorrection 0.2f
#define s2_maxAngularCorrection (8.0f / 180.0f * s2_pi)
#define s2_baumgarte 0.2f
#define s2_timeToSleep 0.5f
#define s2_linearSleepTolerance 0.01f
#define s2_angularSleepTolerance (2.0f / 180.0f * s2_pi)
#define s2_huge (100000.0f)
#define s2_maxBaumgarteVelocity 40.0f
