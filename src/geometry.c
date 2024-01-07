// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "solver2d/geometry.h"

#include "core.h"

#include "solver2d/aabb.h"
#include "solver2d/distance.h"
#include "solver2d/hull.h"
#include "solver2d/math.h"

#include <float.h>

bool s2IsValidRay(const s2RayCastInput* input)
{
	bool isValid = s2IsValidVec2(input->p1) && s2IsValidVec2(input->p2) && s2IsValid(input->maxFraction) &&
				   0.0f <= input->maxFraction && input->maxFraction < s2_huge;
	return isValid;
}

s2Polygon s2MakePolygon(const s2Hull* hull)
{
	S2_ASSERT(hull->count >= 3);

	s2Polygon shape;
	shape.count = hull->count;
	shape.radius = 0.0f;

	// Copy vertices
	for (int32_t i = 0; i < shape.count; ++i)
	{
		shape.vertices[i] = hull->points[i];
	}

	// Compute normals. Ensure the edges have non-zero length.
	for (int32_t i = 0; i < shape.count; ++i)
	{
		int32_t i1 = i;
		int32_t i2 = i + 1 < shape.count ? i + 1 : 0;
		s2Vec2 edge = s2Sub(shape.vertices[i2], shape.vertices[i1]);
		S2_ASSERT(s2Dot(edge, edge) > FLT_EPSILON * FLT_EPSILON);
		shape.normals[i] = s2Normalize(s2CrossVS(edge, 1.0f));
	}

	return shape;
}

s2Polygon s2MakeSquare(float h)
{
	return s2MakeBox(h, h);
}

s2Polygon s2MakeBox(float hx, float hy)
{
	S2_ASSERT(s2IsValid(hx) && hx > 0.0f);
	S2_ASSERT(s2IsValid(hy) && hy > 0.0f);

	s2Polygon shape = {0};
	shape.count = 4;
	shape.vertices[0] = (s2Vec2){-hx, -hy};
	shape.vertices[1] = (s2Vec2){hx, -hy};
	shape.vertices[2] = (s2Vec2){hx, hy};
	shape.vertices[3] = (s2Vec2){-hx, hy};
	shape.normals[0] = (s2Vec2){0.0f, -1.0f};
	shape.normals[1] = (s2Vec2){1.0f, 0.0f};
	shape.normals[2] = (s2Vec2){0.0f, 1.0f};
	shape.normals[3] = (s2Vec2){-1.0f, 0.0f};
	shape.radius = 0.0f;
	return shape;
}

s2Polygon s2MakeRoundedBox(float hx, float hy, float radius)
{
	s2Polygon shape = s2MakeBox(hx, hy);
	shape.radius = radius;
	return shape;
}

s2Polygon s2MakeOffsetBox(float hx, float hy, s2Vec2 center, float angle)
{
	s2Transform xf;
	xf.p = center;
	xf.q = s2MakeRot(angle);

	s2Polygon shape = {0};
	shape.count = 4;
	shape.vertices[0] = s2TransformPoint(xf, (s2Vec2){-hx, -hy});
	shape.vertices[1] = s2TransformPoint(xf, (s2Vec2){hx, -hy});
	shape.vertices[2] = s2TransformPoint(xf, (s2Vec2){hx, hy});
	shape.vertices[3] = s2TransformPoint(xf, (s2Vec2){-hx, hy});
	shape.normals[0] = s2RotateVector(xf.q, (s2Vec2){0.0f, -1.0f});
	shape.normals[1] = s2RotateVector(xf.q, (s2Vec2){1.0f, 0.0f});
	shape.normals[2] = s2RotateVector(xf.q, (s2Vec2){0.0f, 1.0f});
	shape.normals[3] = s2RotateVector(xf.q, (s2Vec2){-1.0f, 0.0f});
	shape.radius = 0.0f;
	return shape;
}

s2Polygon s2MakeCapsule(s2Vec2 p1, s2Vec2 p2, float radius)
{
	s2Polygon shape = {0};
	shape.vertices[0] = p1;
	shape.vertices[1] = p2;

	s2Vec2 axis = s2NormalizeChecked(s2Sub(p2, p1));
	s2Vec2 normal = s2RightPerp(axis);

	shape.normals[0] = normal;
	shape.normals[1] = s2Neg(normal);
	shape.count = 2;
	shape.radius = radius;

	return shape;
}

s2MassData s2ComputeCircleMass(const s2Circle* shape, float density)
{
	float rr = shape->radius * shape->radius;

	s2MassData massData;
	massData.mass = density * s2_pi * rr;
	massData.center = shape->point;

	// inertia about the local origin
	massData.I = massData.mass * (0.5f * rr + s2Dot(shape->point, shape->point));
	return massData;
}

s2MassData s2ComputeCapsuleMass(const s2Capsule* shape, float density)
{
	float radius = shape->radius;
	float rr = radius * radius;
	s2Vec2 p1 = shape->point1;
	s2Vec2 p2 = shape->point2;
	float length = s2Length(s2Sub(p2, p1));
	float ll = length * length;

	s2MassData massData;
	massData.mass = density * (s2_pi * radius + 2.0f * length) * radius;
	massData.center.x = 0.5f * (p1.x + p2.x);
	massData.center.y = 0.5f * (p1.y + p2.y);

	// two offset half circles, both halves add up to full circle and each half is offset by half length
	// half circles = 2 * (1/4 * radius*radius + 1/4 * length*length)
	// rectangle = (width*width * length*length)/12
	float circleInertia = 0.5f * (rr + ll);
	float boxInertia = (4.0f * rr + ll) / 12.0f;
	massData.I = massData.mass * (circleInertia + boxInertia);

	return massData;
}

s2MassData s2ComputePolygonMass(const s2Polygon* shape, float density)
{
	// Polygon mass, centroid, and inertia.
	// Let rho be the polygon density in mass per unit area.
	// Then:
	// mass = rho * int(dA)
	// centroid.x = (1/mass) * rho * int(x * dA)
	// centroid.y = (1/mass) * rho * int(y * dA)
	// I = rho * int((x*x + y*y) * dA)
	//
	// We can compute these integrals by summing all the integrals
	// for each triangle of the polygon. To evaluate the integral
	// for a single triangle, we make a change of variables to
	// the (u,v) coordinates of the triangle:
	// x = x0 + e1x * u + e2x * v
	// y = y0 + e1y * u + e2y * v
	// where 0 <= u && 0 <= v && u + v <= 1.
	//
	// We integrate u from [0,1-v] and then v from [0,1].
	// We also need to use the Jacobian of the transformation:
	// D = cross(e1, e2)
	//
	// Simplification: triangle centroid = (1/3) * (p1 + p2 + p3)
	//
	// The rest of the derivation is handled by computer algebra.

	S2_ASSERT(shape->count > 0);

	if (shape->count == 1)
	{
		s2Circle circle;
		circle.point = shape->vertices[0];
		circle.radius = shape->radius;
		return s2ComputeCircleMass(&circle, density);
	}

	if (shape->count == 2)
	{
		s2Capsule capsule;
		capsule.point1 = shape->vertices[0];
		capsule.point2 = shape->vertices[1];
		capsule.radius = shape->radius;
		return s2ComputeCapsuleMass(&capsule, density);
	}

	s2Vec2 vertices[s2_maxPolygonVertices];
	int32_t count = shape->count;
	float radius = shape->radius;

	if (radius > 0.0f)
	{
		// Push out vertices according to radius. This improves
		// the mass accuracy, especially the rotational inertia.
		for (int32_t i = 0; i < count; ++i)
		{
			int32_t j = i == 0 ? count - 1 : i - 1;
			s2Vec2 n1 = shape->normals[j];
			s2Vec2 n2 = shape->normals[i];

			s2Vec2 mid = s2Normalize(s2Add(n1, n2));
			s2Vec2 t1 = {-n1.y, n1.x};
			float sinHalfAngle = s2Cross(mid, t1);

			float offset = radius;
			if (sinHalfAngle > FLT_EPSILON)
			{
				offset = radius / sinHalfAngle;
			}

			vertices[i] = s2MulAdd(shape->vertices[i], offset, mid);
		}
	}
	else
	{
		for (int32_t i = 0; i < count; ++i)
		{
			vertices[i] = shape->vertices[i];
		}
	}

	s2Vec2 center = {0.0f, 0.0f};
	float area = 0.0f;
	float I = 0.0f;

	// Get a reference point for forming triangles.
	// Use the first vertex to reduce round-off errors.
	s2Vec2 r = vertices[0];

	const float inv3 = 1.0f / 3.0f;

	for (int32_t i = 1; i < count - 1; ++i)
	{
		// Triangle edges
		s2Vec2 e1 = s2Sub(vertices[i], r);
		s2Vec2 e2 = s2Sub(vertices[i + 1], r);

		float D = s2Cross(e1, e2);

		float triangleArea = 0.5f * D;
		area += triangleArea;

		// Area weighted centroid, r at origin
		center = s2MulAdd(center, triangleArea * inv3, s2Add(e1, e2));

		float ex1 = e1.x, ey1 = e1.y;
		float ex2 = e2.x, ey2 = e2.y;

		float intx2 = ex1 * ex1 + ex2 * ex1 + ex2 * ex2;
		float inty2 = ey1 * ey1 + ey2 * ey1 + ey2 * ey2;

		I += (0.25f * inv3 * D) * (intx2 + inty2);
	}

	s2MassData massData;

	// Total mass
	massData.mass = density * area;

	// Center of mass, shift back from origin at r
	S2_ASSERT(area > FLT_EPSILON);
	float invArea = 1.0f / area;
	center.x *= invArea;
	center.y *= invArea;
	massData.center = s2Add(r, center);

	// Inertia tensor relative to the local origin (point s).
	massData.I = density * I;

	// Shift to center of mass then to original body origin.
	massData.I += massData.mass * (s2Dot(massData.center, massData.center) - s2Dot(center, center));

	return massData;
}

s2Box s2ComputeCircleAABB(const s2Circle* shape, s2Transform xf)
{
	s2Vec2 p = s2TransformPoint(xf, shape->point);
	float r = shape->radius;

	s2Box aabb = {{p.x - r, p.y - r}, {p.x + r, p.y + r}};
	return aabb;
}

s2Box s2ComputeCapsuleAABB(const s2Capsule* shape, s2Transform xf)
{
	s2Vec2 v1 = s2TransformPoint(xf, shape->point1);
	s2Vec2 v2 = s2TransformPoint(xf, shape->point2);

	s2Vec2 r = {shape->radius, shape->radius};
	s2Vec2 lower = s2Sub(s2Min(v1, v2), r);
	s2Vec2 upper = s2Add(s2Max(v1, v2), r);

	s2Box aabb = {lower, upper};
	return aabb;
}

s2Box s2ComputePolygonAABB(const s2Polygon* shape, s2Transform xf)
{
	S2_ASSERT(shape->count > 0);
	s2Vec2 lower = s2TransformPoint(xf, shape->vertices[0]);
	s2Vec2 upper = lower;

	for (int32_t i = 1; i < shape->count; ++i)
	{
		s2Vec2 v = s2TransformPoint(xf, shape->vertices[i]);
		lower = s2Min(lower, v);
		upper = s2Max(upper, v);
	}

	s2Vec2 r = {shape->radius, shape->radius};
	lower = s2Sub(lower, r);
	upper = s2Add(upper, r);

	s2Box aabb = {lower, upper};
	return aabb;
}

s2Box s2ComputeSegmentAABB(const s2Segment* shape, s2Transform xf)
{
	s2Vec2 v1 = s2TransformPoint(xf, shape->point1);
	s2Vec2 v2 = s2TransformPoint(xf, shape->point2);

	s2Vec2 lower = s2Min(v1, v2);
	s2Vec2 upper = s2Max(v1, v2);

	s2Box aabb = {lower, upper};
	return aabb;
}

bool s2PointInCircle(s2Vec2 point, const s2Circle* shape)
{
	s2Vec2 center = shape->point;
	return s2DistanceSquared(point, center) <= shape->radius * shape->radius;
}

bool s2PointInCapsule(s2Vec2 point, const s2Capsule* shape)
{
	float rr = shape->radius * shape->radius;
	s2Vec2 p1 = shape->point1;
	s2Vec2 p2 = shape->point2;

	s2Vec2 d = s2Sub(p2, p1);
	float dd = s2Dot(d, d);
	if (dd == 0.0f)
	{
		// Capsule is really a circle
		return s2DistanceSquared(point, p1) <= rr;
	}

	// Get closest point on capsule segment
	// c = p1 + t * d
	// dot(point - c, d) = 0
	// dot(point - p1 - t * d, d) = 0
	// t = dot(point - p1, d) / dot(d, d)
	float t = s2Dot(s2Sub(point, p1), d) / dd;
	t = S2_CLAMP(t, 0.0f, 1.0f);
	s2Vec2 c = s2MulAdd(p1, t, d);

	// Is query point within radius around closest point?
	return s2DistanceSquared(point, c) <= rr;
}

bool s2PointInPolygon(s2Vec2 point, const s2Polygon* shape)
{
	s2DistanceInput input = {0};
	input.proxyA = s2MakeProxy(shape->vertices, shape->count, 0.0f);
	input.proxyB = s2MakeProxy(&point, 1, 0.0f);
	input.transformA = s2Transform_identity;
	input.transformB = s2Transform_identity;
	input.useRadii = false;

	s2DistanceCache cache = {0};
	s2DistanceOutput output = s2ShapeDistance(&cache, &input);

	return output.distance <= shape->radius;
}

// Precision Improvements for Ray / Sphere Intersection - Ray Tracing Gems 2019
// http://www.codercorner.com/blog/?p=321
s2RayCastOutput s2RayCastCircle(const s2RayCastInput* input, const s2Circle* shape)
{
	S2_ASSERT(s2IsValidRay(input));

	s2Vec2 p = shape->point;

	s2RayCastOutput output = {0};

	// Shift ray so circle center is the origin
	s2Vec2 s = s2Sub(input->p1, p);
	float length;
	s2Vec2 d = s2GetLengthAndNormalize(&length, s2Sub(input->p2, input->p1));
	if (length == 0.0f)
	{
		// zero length ray
		return output;
	}

	// Find closest point on ray to origin

	// solve: dot(s + t * d, d) = 0
	float t = -s2Dot(s, d);

	// c is the closest point on the line to the origin
	s2Vec2 c = s2MulAdd(s, t, d);

	float cc = s2Dot(c, c);
	float r = shape->radius;
	float rr = r * r;

	if (cc > rr)
	{
		// closest point is outside the circle
		return output;
	}

	// Pythagorus
	float h = sqrtf(rr - cc);

	float fraction = t - h;

	if (fraction < 0.0f || input->maxFraction * length < fraction)
	{
		// outside the range of the ray segment
		return output;
	}

	s2Vec2 hitPoint = s2MulAdd(s, fraction, d);

	output.fraction = fraction / length;
	output.normal = s2Normalize(hitPoint);
	output.point = s2MulAdd(p, shape->radius, output.normal);
	output.hit = true;

	return output;
}

s2RayCastOutput s2RayCastCapsule(const s2RayCastInput* input, const s2Capsule* shape)
{
	S2_ASSERT(s2IsValidRay(input));

	s2RayCastOutput output = {0};

	s2Vec2 v1 = shape->point1;
	s2Vec2 v2 = shape->point2;

	s2Vec2 e = s2Sub(v2, v1);

	float capsuleLength;
	s2Vec2 a = s2GetLengthAndNormalize(&capsuleLength, e);

	if (capsuleLength < FLT_EPSILON)
	{
		// Capsule is really a circle
		s2Circle circle = {v1, shape->radius};
		return s2RayCastCircle(input, &circle);
	}

	s2Vec2 p1 = input->p1;
	s2Vec2 p2 = input->p2;

	// Ray from capsule start to ray start
	s2Vec2 q = s2Sub(p1, v1);
	float qa = s2Dot(q, a);

	// Vector to ray start that is perpendicular to capsule axis
	s2Vec2 qp = s2MulAdd(q, -qa, a);

	float radius = shape->radius;

	// Does the ray start within the infinite length capsule?
	if (s2Dot(qp, qp) < radius * radius)
	{
		if (qa < 0.0f)
		{
			// start point behind capsule segment
			s2Circle circle = {v1, shape->radius};
			return s2RayCastCircle(input, &circle);
		}

		if (qa > 1.0f)
		{
			// start point ahead of capsule segment
			s2Circle circle = {v2, shape->radius};
			return s2RayCastCircle(input, &circle);
		}

		// ray starts inside capsule -> no hit
		return output;
	}

	// Perpendicular to capsule axis, pointing right
	s2Vec2 n = {a.y, -a.x};

	float rayLength;
	s2Vec2 u = s2GetLengthAndNormalize(&rayLength, s2Sub(p2, p1));

	// Intersect ray with infinite length capsule
	// v1 + radius * n + s1 * a = p1 + s2 * u
	// v1 - radius * n + s1 * a = p1 + s2 * u

	// s1 * a - s2 * u = b
	// b = q - radius * ap
	// or
	// b = q + radius * ap

	// Cramer's rule [a -u]
	float den = -a.x * u.y + u.x * a.y;
	if (-FLT_EPSILON < den && den < FLT_EPSILON)
	{
		// Ray is parallel to capsule and outside infinite length capsule
		return output;
	}

	s2Vec2 b1 = s2MulSub(q, radius, n);
	s2Vec2 bz2 = s2MulAdd(q, radius, n);

	float invDen = 1.0f / den;

	// Cramer's rule [a b1]
	float s21 = (a.x * b1.y - b1.x * a.y) * invDen;

	// Cramer's rule [a bz2]
	float s22 = (a.x * bz2.y - bz2.x * a.y) * invDen;

	float s2;
	s2Vec2 b;
	if (s21 < s22)
	{
		s2 = s21;
		b = b1;
	}
	else
	{
		s2 = s22;
		b = bz2;
		n = s2Neg(n);
	}

	if (s2 < 0.0f || input->maxFraction * rayLength < s2)
	{
		return output;
	}

	// Cramer's rule [b -u]
	float s1 = (-b.x * u.y + u.x * b.y) * invDen;

	if (s1 < 0.0f)
	{
		// ray passes behind capsule segment
		s2Circle circle = {v1, shape->radius};
		return s2RayCastCircle(input, &circle);
	}
	else if (capsuleLength < s1)
	{
		// ray passes ahead of capsule segment
		s2Circle circle = {v2, shape->radius};
		return s2RayCastCircle(input, &circle);
	}
	else
	{
		// ray hits capsule side
		output.fraction = s2 / rayLength;
		output.point = s2Add(s2Lerp(v1, v2, s1 / capsuleLength), s2MulSV(shape->radius, n));
		output.normal = n;
		output.hit = true;
		return output;
	}
}

// Ray vs line segment
s2RayCastOutput s2RayCastSegment(const s2RayCastInput* input, const s2Segment* shape)
{
	{
		// Put the ray into the edge's frame of reference.
		s2Vec2 p1 = input->p1;
		s2Vec2 p2 = input->p2;
		s2Vec2 d = s2Sub(p2, p1);

		s2Vec2 v1 = shape->point1;
		s2Vec2 v2 = shape->point2;
		s2Vec2 e = s2Sub(v2, v1);

		s2RayCastOutput output = {0};

		float length;
		s2Vec2 eUnit = s2GetLengthAndNormalize(&length, e);
		if (length == 0.0f)
		{
			return output;
		}

		// Normal points to the right, looking from v1 towards v2
		s2Vec2 normal = {eUnit.y, -eUnit.x};

		// Intersect ray with infinite segment using normal
		// Similar to intersecting a ray with an infinite plane
		// p = p1 + t * d
		// dot(normal, p - v1) = 0
		// dot(normal, p1 - v1) + t * dot(normal, d) = 0
		float numerator = s2Dot(normal, s2Sub(v1, p1));
		float denominator = s2Dot(normal, d);

		if (denominator == 0.0f)
		{
			// parallel
			return output;
		}

		float t = numerator / denominator;
		if (t < 0.0f || input->maxFraction < t)
		{
			// out of ray range
			return output;
		}

		// Intersection point on infinite segment
		s2Vec2 p = s2MulAdd(p1, t, d);

		// Compute position of p along segment
		// p = v1 + s * e
		// s = dot(p - v1, e) / dot(e, e)

		float s = s2Dot(s2Sub(p, v1), eUnit);
		if (s < 0.0f || length < s)
		{
			// out of segment range
			return output;
		}

		if (numerator > 0.0f)
		{
			normal = s2Neg(normal);
		}

		output.fraction = t;
		output.normal = normal;
		output.hit = true;

		return output;
	}
}

s2RayCastOutput s2RayCastPolygon(const s2RayCastInput* input, const s2Polygon* shape)
{
	S2_ASSERT(s2IsValidRay(input));

	{
		// Put the ray into the polygon's frame of reference.
		s2Vec2 p1 = input->p1;
		s2Vec2 p2 = input->p2;
		s2Vec2 d = s2Sub(p2, p1);

		float lower = 0.0f, upper = input->maxFraction;

		int32_t index = -1;

		s2RayCastOutput output = {0};

		for (int32_t i = 0; i < shape->count; ++i)
		{
			// p = p1 + a * d
			// dot(normal, p - v) = 0
			// dot(normal, p1 - v) + a * dot(normal, d) = 0
			float numerator = s2Dot(shape->normals[i], s2Sub(shape->vertices[i], p1));
			float denominator = s2Dot(shape->normals[i], d);

			if (denominator == 0.0f)
			{
				if (numerator < 0.0f)
				{
					return output;
				}
			}
			else
			{
				// Note: we want this predicate without division:
				// lower < numerator / denominator, where denominator < 0
				// Since denominator < 0, we have to flip the inequality:
				// lower < numerator / denominator <==> denominator * lower > numerator.
				if (denominator < 0.0f && numerator < lower * denominator)
				{
					// Increase lower.
					// The segment enters this half-space.
					lower = numerator / denominator;
					index = i;
				}
				else if (denominator > 0.0f && numerator < upper * denominator)
				{
					// Decrease upper.
					// The segment exits this half-space.
					upper = numerator / denominator;
				}
			}

			// The use of epsilon here causes the S2_ASSERT on lower to trip
			// in some cases. Apparently the use of epsilon was to make edge
			// shapes work, but now those are handled separately.
			// if (upper < lower - s2_epsilon)
			if (upper < lower)
			{
				return output;
			}
		}

		S2_ASSERT(0.0f <= lower && lower <= input->maxFraction);

		if (index >= 0)
		{
			output.fraction = lower;
			output.normal = shape->normals[index];
			output.point = s2Lerp(p1, p2, output.fraction);
			output.hit = true;
		}

		return output;
	}
}
