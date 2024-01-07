// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "allocate.h"

#include "core.h"

#if defined(S2_COMPILER_MSVC)
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#include <stdlib.h>
#else
#include <stdlib.h>
#endif

#include <stdint.h>

void* s2Alloc(uint32_t size)
{
	uint32_t size16 = ((size - 1) | 0xF) + 1;
#ifdef S2_PLATFORM_WINDOWS
	void* ptr = _aligned_malloc(size16, 16);
#else
	void* ptr = aligned_alloc(16, size16);
#endif

	return ptr;
}

void s2Free(void* mem, uint32_t size)
{
	if (mem == NULL)
	{
		return;
	}

#ifdef S2_PLATFORM_WINDOWS
		_aligned_free(mem);
#else
		free(mem);
#endif
}
