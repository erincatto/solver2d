// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "types.h"

#include <math.h>

#define S2_MIN(A, B) ((A) < (B) ? (A) : (B))
#define S2_MAX(A, B) ((A) > (B) ? (A) : (B))
#define S2_ABS(A) ((A) > 0.0f ? (A) : -(A))
#define S2_CLAMP(A, B, C) S2_MIN(S2_MAX(A, B), C)

static const s2Vec2 s2Vec2_zero = {0.0f, 0.0f};
static const s2Rot s2Rot_identity = {0.0f, 1.0f};
static const s2Transform s2Transform_identity = {{0.0f, 0.0f}, {0.0f, 1.0f}};
static const s2Mat22 s2Mat22_zero = {{0.0f, 0.0f}, {0.0f, 0.0f}};

#ifdef __cplusplus
extern "C"
{
#endif

bool s2IsValid(float a);
bool s2IsValidVec2(s2Vec2 v);

/// Convert this vector into a unit vector
s2Vec2 s2Normalize(s2Vec2 v);

/// This asserts of the vector is too short
s2Vec2 s2NormalizeChecked(s2Vec2 v);

s2Vec2 s2GetLengthAndNormalize(float* length, s2Vec2 v);

#ifdef __cplusplus
}
#endif

/// Make a vector
static inline s2Vec2 s2MakeVec2(float x, float y)
{
	return S2_LITERAL(s2Vec2){x, y};
}

/// Vector dot product
static inline float s2Dot(s2Vec2 a, s2Vec2 b)
{
	return a.x * b.x + a.y * b.y;
}

/// Vector cross product. In 2D this yields a scalar.
static inline float s2Cross(s2Vec2 a, s2Vec2 b)
{
	return a.x * b.y - a.y * b.x;
}

/// Perform the cross product on a vector and a scalar. In 2D this produces
/// a vector.
static inline s2Vec2 s2CrossVS(s2Vec2 v, float s)
{
	return S2_LITERAL(s2Vec2){s * v.y, -s * v.x};
}

/// Perform the cross product on a scalar and a vector. In 2D this produces
/// a vector.
static inline s2Vec2 s2CrossSV(float s, s2Vec2 v)
{
	return S2_LITERAL(s2Vec2){-s * v.y, s * v.x};
}

/// Get a right pointing perpendicular vector. Equivalent to s2CrossVS(v, 1.0f).
static inline s2Vec2 s2RightPerp(s2Vec2 v)
{
	return S2_LITERAL(s2Vec2){v.y, -v.x};
}

/// Vector addition
static inline s2Vec2 s2Add(s2Vec2 a, s2Vec2 b)
{
	return S2_LITERAL(s2Vec2){a.x + b.x, a.y + b.y};
}

/// Vector subtraction
static inline s2Vec2 s2Sub(s2Vec2 a, s2Vec2 b)
{
	return S2_LITERAL(s2Vec2){a.x - b.x, a.y - b.y};
}

/// Vector subtraction
static inline s2Vec2 s2Neg(s2Vec2 a)
{
return S2_LITERAL(s2Vec2){-a.x, -a.y};
}

/// Vector linear interpolation
static inline s2Vec2 s2Lerp(s2Vec2 a, s2Vec2 b, float t)
{
	return S2_LITERAL(s2Vec2){a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)};
}

/// Component-wise multiplication
static inline s2Vec2 s2Mul(s2Vec2 a, s2Vec2 b)
{
	return S2_LITERAL(s2Vec2){a.x * b.x, a.y * b.y};
}

/// Multiply a scalar and vector
static inline s2Vec2 s2MulSV(float s, s2Vec2 v)
{
	return S2_LITERAL(s2Vec2){s * v.x, s * v.y};
}

/// a + s * b
static inline s2Vec2 s2MulAdd(s2Vec2 a, float s, s2Vec2 b)
{
	return S2_LITERAL(s2Vec2){a.x + s * b.x, a.y + s * b.y};
}

/// a - s * b
static inline s2Vec2 s2MulSub(s2Vec2 a, float s, s2Vec2 b)
{
	return S2_LITERAL(s2Vec2){a.x - s * b.x, a.y - s * b.y};
}

/// Component-wise absolute vector
static inline s2Vec2 s2Abs(s2Vec2 a)
{
	s2Vec2 b;
	b.x = S2_ABS(a.x);
	b.y = S2_ABS(a.y);
	return b;
}

/// Component-wise absolute vector
static inline s2Vec2 s2Min(s2Vec2 a, s2Vec2 b)
{
	s2Vec2 c;
	c.x = S2_MIN(a.x, b.x);
	c.y = S2_MIN(a.y, b.y);
	return c;
}

/// Component-wise absolute vector
static inline s2Vec2 s2Max(s2Vec2 a, s2Vec2 b)
{
	s2Vec2 c;
	c.x = S2_MAX(a.x, b.x);
	c.y = S2_MAX(a.y, b.y);
	return c;
}

/// Component-wise clamp vector so v into the range [a, b]
static inline s2Vec2 s2Clamp(s2Vec2 v, s2Vec2 a, s2Vec2 b)
{
	s2Vec2 c;
	c.x = S2_CLAMP(v.x, a.x, b.x);
	c.y = S2_CLAMP(v.y, a.y, b.y);
	return c;
}

/// Get the length of this vector (the norm).
static inline float s2Length(s2Vec2 v)
{
	return sqrtf(v.x * v.x + v.y * v.y);
}

/// Get the length of this vector (the norm).
static inline float s2LengthSquared(s2Vec2 v)
{
	return v.x * v.x + v.y * v.y;
}

static inline float s2Distance(s2Vec2 a, s2Vec2 b)
{
	float dx = b.x - a.x;
	float dy = b.y - a.y;
	return sqrtf(dx * dx + dy * dy);
}

/// Get the length of this vector (the norm).
static inline float s2DistanceSquared(s2Vec2 a, s2Vec2 b)
{
	s2Vec2 c = {b.x - a.x, b.y - a.y};
	return c.x * c.x + c.y * c.y;
}

/// Set using an angle in radians.
static inline s2Rot s2MakeRot(float angle)
{
	s2Rot q = {sinf(0.5f * angle), cosf(0.5f * angle)};
	return q;
}

static inline s2Rot s2NormalizeRot(s2Rot q)
{
	float mag = sqrtf(q.c * q.c + q.s * q.s);
	float invMag = mag > 0.0 ? 1.0f / mag : 0.0f;
	s2Rot qn = {q.c * invMag, q.s * invMag};
	return qn;
}

static inline s2Rot s2IntegrateRot(s2Rot q1, float omegah)
{
	// dc/dt = -omega * sin(t)
	// ds/dt = omega * cos(t)
	// c2 = c1 - omega * h * s1
	// s2 = s2 + omega * h * c1
	s2Rot q2 = {q1.c - omegah * q1.s, q1.s + omegah * q1.c};
	return s2NormalizeRot(q2);
}

/// Get the angle in radians
static inline float s2Rot_GetAngle(s2Rot q)
{
	return 2.0f * atan2f(q.s, q.c);
}

/// Get the x-axis
static inline s2Vec2 s2Rot_GetXAxis(s2Rot q)
{
	s2Vec2 v = {q.c, q.s};
	return v;
}

/// Get the y-axis
static inline s2Vec2 s2Rot_GetYAxis(s2Rot q)
{
	s2Vec2 v = {-q.s, q.c};
	return v;
}

/// Multiply two rotations: q * r
static inline s2Rot s2MulRot(s2Rot q, s2Rot r)
{
	// [qc -qs] * [rc -rs] = [qc*rc-qs*rs -qc*rs-qs*rc]
	// [qs  qc]   [rs  rc]   [qs*rc+qc*rs -qs*rs+qc*rc]
	// s = qs * rc + qc * rs
	// c = qc * rc - qs * rs
	s2Rot qr;
	qr.s = q.s * r.c + q.c * r.s;
	qr.c = q.c * r.c - q.s * r.s;
	return qr;
}

/// Transpose multiply two rotations: qT * r
static inline s2Rot s2InvMulRot(s2Rot q, s2Rot r)
{
	// [ qc qs] * [rc -rs] = [qc*rc+qs*rs -qc*rs+qs*rc]
	// [-qs qc]   [rs  rc]   [-qs*rc+qc*rs qs*rs+qc*rc]
	// s = qc * rs - qs * rc
	// c = qc * rc + qs * rs
	s2Rot qr;
	qr.s = q.c * r.s - q.s * r.c;
	qr.c = q.c * r.c + q.s * r.s;
	return qr;
}

/// Rotate a vector
static inline s2Vec2 s2RotateVector(s2Rot q, s2Vec2 v)
{
	// s2CrossSV(s,v) = {-s * v.y, s * v.x}
	// v + 2.0f * s2CrossSV(q.z, s2CrossSV(q.z, v) + q.w * v)
	// v + 2.0f * s2CrossSV(q.z, {-q.z * v.y, q.z * v.x} + q.w * v)
	// v + 2.0f * s2CrossSV(q.z, {-q.z * v.y + q.w * v.x, q.z * v.x + q.w * v.y})
	// v + 2.0f * {-q.z * (q.z * v.x + q.w * v.y), q.z * (-q.z * v.y + q.w * v.x)}
	// {v.x - 2.0f * q.z * (q.z * v.x + q.w * v.y), v.y - 2.0f * q.z * (q.z * v.y - q.w * v.x)}

	return S2_LITERAL(s2Vec2){q.c * v.x - q.s * v.y, q.s * v.x + q.c * v.y};
}

/// Inverse rotate a vector
static inline s2Vec2 s2InvRotateVector(s2Rot q, s2Vec2 v)
{
	return S2_LITERAL(s2Vec2){q.c * v.x + q.s * v.y, -q.s * v.x + q.c * v.y};
}

/// Transform a point (e.g. local space to world space)
static inline s2Vec2 s2TransformPoint(s2Transform xf, const s2Vec2 p)
{
	float x = (xf.q.c * p.x - xf.q.s * p.y) + xf.p.x;
	float y = (xf.q.s * p.x + xf.q.c * p.y) + xf.p.y;

	return S2_LITERAL(s2Vec2){x, y};
}

// Inverse transform a point (e.g. world space to local space)
static inline s2Vec2 s2InvTransformPoint(s2Transform xf, const s2Vec2 p)
{
	float vx = p.x - xf.p.x;
	float vy = p.y - xf.p.y;
	return S2_LITERAL(s2Vec2){xf.q.c * vx + xf.q.s * vy, -xf.q.s * vx + xf.q.c * vy};
}

// v2 = A.q.Rot(B.q.Rot(v1) + B.p) + A.p
//    = (A.q * B.q).Rot(v1) + A.q.Rot(B.p) + A.p
static inline s2Transform s2MulTransforms(s2Transform A, s2Transform B)
{
	s2Transform C;
	C.q = s2MulRot(A.q, B.q);
	C.p = s2Add(s2RotateVector(A.q, B.p), A.p);
	return C;
}

// v2 = A.q' * (B.q * v1 + B.p - A.p)
//    = A.q' * B.q * v1 + A.q' * (B.p - A.p)
static inline s2Transform s2InvMulTransforms(s2Transform A, s2Transform B)
{
	s2Transform C;
	C.q = s2InvMulRot(A.q, B.q);
	C.p = s2InvRotateVector(A.q, s2Sub(B.p, A.p));
	return C;
}

static inline s2Vec2 s2MulMV(s2Mat22 A, s2Vec2 v)
{
	s2Vec2 u = {A.cx.x * v.x + A.cy.x * v.y, A.cx.y * v.x + A.cy.y * v.y};
	return u;
}

static inline s2Mat22 s2GetInverse22(s2Mat22 A)
{
	float a = A.cx.x, b = A.cy.x, c = A.cx.y, d = A.cy.y;
	s2Mat22 B;
	float det = a * d - b * c;
	if (det != 0.0f)
	{
		det = 1.0f / det;
	}
	B.cx.x = det * d;
	B.cy.x = -det * b;
	B.cx.y = -det * c;
	B.cy.y = det * a;
	return B;
}

/// Solve A * x = b, where b is a column vector. This is more efficient
/// than computing the inverse in one-shot cases.
static inline s2Vec2 s2Solve22(s2Mat22 A, s2Vec2 b)
{
	float a11 = A.cx.x, a12 = A.cy.x, a21 = A.cx.y, a22 = A.cy.y;
	float det = a11 * a22 - a12 * a21;
	if (det != 0.0f)
	{
		det = 1.0f / det;
	}
	s2Vec2 x = {det * (a22 * b.x - a12 * b.y), det * (a11 * b.y - a21 * b.x)};
	return x;
}
