// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once
#include <stdint.h>

void* s2Alloc(uint32_t size);
void s2Free(void* mem, uint32_t size);
