// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include <stdint.h>

#if defined(_DEBUG)
#include <stdio.h>
#endif

// Define platform
#if defined(_WIN64)
#define S2_PLATFORM_WINDOWS
#elif defined(__ANDROID__)
#define S2_PLATFORM_ANDROID
#elif defined(__linux__)
#define S2_PLATFORM_LINUX
#elif defined(__APPLE__)
#include <TargetConditionals.h>
#if defined(TARGET_OS_IPHONE) && !TARGET_OS_IPHONE
#define S2_PLATFORM_MACOS
#else
#define S2_PLATFORM_IOS
#endif
#elif defined(__EMSCRIPTEN__)
#define S2_PLATFORM_WASM
#else
#error Unsupported platform
#endif

// Define CPU
#if defined(__x86_64__) || defined(_M_X64)
#define S2_CPU_X64
#elif defined(__aarch64__) || defined(_M_ARM64)
#define S2_CPU_ARM
#elif defined(JPH_PLATFORM_WASM)
#define S2_CPU_WASM
#else
#error Unsupported CPU
#endif

// Define compiler
#if defined(__clang__)
#define S2_COMPILER_CLANG
#elif defined(__GNUC__)
#define S2_COMPILER_GCC
#elif defined(_MSC_VER)
#define S2_COMPILER_MSVC
#endif

#if defined(S2_PLATFORM_WINDOWS)
#define S2_BREAKPOINT __debugbreak()
#elif defined(S2_PLATFORM_LINUX) || defined(S2_PLATFORM_ANDROID) || defined(S2_PLATFORM_MACOS) || defined(S2_PLATFORM_IOS)
#if defined(S2_CPU_X64)
#define S2_BREAKPOINT __asm volatile("int $0x3")
#elif defined(S2_CPU_ARM)
#define S2_BREAKPOINT __builtin_trap()
#endif
#elif defined(S2_PLATFORM_WASM)
#define S2_BREAKPOINT                                                                                                                      \
	do                                                                                                                                     \
	{                                                                                                                                      \
	} while (0)
#else
#error Unknown platform
#endif

#if defined(_DEBUG)
#define S2_ASSERT(condition)                                                                                                               \
	do                                                                                                                                     \
	{                                                                                                                                      \
		if (!(condition) && printf("ASSERTION: %s, %s, line %d\n", #condition, __FILE__, (int)__LINE__))                                                      \
			S2_BREAKPOINT;                                                                                                                 \
	} while (0)
#else
#define S2_ASSERT(...) ((void)0)
#endif

#if defined(_DEBUG) || 0
#define S2_VALIDATE 1
#else
#define S2_VALIDATE 0
#endif

