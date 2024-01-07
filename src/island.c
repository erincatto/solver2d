// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#include "island.h"

#include "array.h"
#include "body.h"
#include "contact.h"
#include "contact_solver.h"
#include "core.h"
#include "joint.h"
#include "shape.h"
#include "solver2d/aabb.h"
#include "solver2d/callbacks.h"
#include "solver2d/timer.h"
#include "solver_data.h"
#include "stack_allocator.h"
#include "world.h"

#include <float.h>
#include <stdlib.h>

/*
Position Correction Notes
=========================
I tried the several algorithms for position correction of the 2D revolute joint.
I looked at these systems:
- simple pendulum (1m diameter sphere on massless 5m stick) with initial angular velocity of 100 rad/s.
- suspension bridge with 30 1m long planks of length 1m.
- multi-link chain with 30 1m long links.

Here are the algorithms:

Baumgarte - A fraction of the position error is added to the velocity error. There is no
separate position solver.

Pseudo Velocities - After the velocity solver and position integration,
the position error, Jacobian, and effective mass are recomputed. Then
the velocity constraints are solved with pseudo velocities and a fraction
of the position error is added to the pseudo velocity error. The pseudo
velocities are initialized to zero and there is no warm-starting. After
the position solver, the pseudo velocities are added to the positions.
This is also called the First Order World method or the Position LCP method.

Modified Nonlinear Gauss-Seidel (NGS) - Like Pseudo Velocities except the
position error is re-computed for each constraint and the positions are updated
after the constraint is solved. The radius vectors (aka Jacobians) are
re-computed too (otherwise the algorithm has horrible instability). The pseudo
velocity states are not needed because they are effectively zero at the beginning
of each iteration. Since we have the current position error, we allow the
iterations to terminate early if the error becomes smaller than s2_linearSlop.

Full NGS or just NGS - Like Modified NGS except the effective mass are re-computed
each time a constraint is solved.

Here are the results:
Baumgarte - this is the cheapest algorithm but it has some stability problems,
especially with the bridge. The chain links separate easily close to the root
and they jitter as they struggle to pull together. This is one of the most common
methods in the field. The big drawback is that the position correction artificially
affects the momentum, thus leading to instabilities and false bounce. I used a
bias factor of 0.2. A larger bias factor makes the bridge less stable, a smaller
factor makes joints and contacts more spongy.

Pseudo Velocities - the is more stable than the Baumgarte method. The bridge is
stable. However, joints still separate with large angular velocities. Drag the
simple pendulum in a circle quickly and the joint will separate. The chain separates
easily and does not recover. I used a bias factor of 0.2. A larger value lead to
the bridge collapsing when a heavy cube drops on it.

Modified NGS - this algorithm is better in some ways than Baumgarte and Pseudo
Velocities, but in other ways it is worse. The bridge and chain are much more
stable, but the simple pendulum goes unstable at high angular velocities.

Full NGS - stable in all tests. The joints display good stiffness. The bridge
still sags, but this is better than infinite forces.

Recommendations
Pseudo Velocities are not really worthwhile because the bridge and chain cannot
recover from joint separation. In other cases the benefit over Baumgarte is small.

Modified NGS is not a robust method for the revolute joint due to the violent
instability seen in the simple pendulum. Perhaps it is viable with other constraint
types, especially scalar constraints where the effective mass is a scalar.

This leaves Baumgarte and Full NGS. Baumgarte has small, but manageable instabilities
and is very fast. I don't think we can escape Baumgarte, especially in highly
demanding cases where high constraint fidelity is not needed.

Full NGS is robust and easy on the eyes. I recommend this as an option for
higher fidelity simulation and certainly for suspension bridges and long chains.
Full NGS might be a good choice for ragdolls, especially motorized ragdolls where
joint separation can be problematic. The number of NGS iterations can be reduced
for better performance without harming robustness much.

Each joint in a can be handled differently in the position solver. So I recommend
a system where the user can select the algorithm on a per joint basis. I would
probably default to the slower Full NGS and let the user select the faster
Baumgarte method in performance critical scenarios.
*/

/*
2D Rotation

R = [cos(theta) -sin(theta)]
	[sin(theta) cos(theta) ]

thetaDot = omega

Let q1 = cos(theta), q2 = sin(theta).
R = [q1 -q2]
	[q2  q1]

q1Dot = -thetaDot * q2
q2Dot = thetaDot * q1

q1_new = q1_old - dt * w * q2
q2_new = q2_old + dt * w * q1
then normalize.

This might be faster than computing sin+cos.
However, we can compute sin+cos of the same angle fast.
*/

void s2SolveWorld(s2World* world, s2StepContext* context)
{
	s2Body* bodies = world->bodies;
	s2Joint* joints = world->joints;
	s2Vec2 gravity = world->gravity;

	float h = context->dt;

	// Integrate velocities and apply damping. Initialize the body state.
	int bodyCapacity = world->bodyPool.capacity;
	for (int i = 0; i < bodyCapacity; ++i)
	{
		s2Body* body = bodies + i;
		if (s2ObjectIsFree(&body->object))
		{
			continue;
		}

		float invMass = body->invMass;
		float invI = body->invI;

		if (body->type == s2_dynamicBody)
		{
			s2Vec2 v = body->linearVelocity;
			float w = body->angularVelocity;

			// Integrate velocities
			v = s2Add(v, s2MulSV(h * invMass, s2MulAdd(body->force, body->gravityScale * body->mass, gravity)));
			w = w + h * invI * body->torque;

			// Apply damping.
			// ODE: dv/dt + c * v = 0
			// Solution: v(t) = v0 * exp(-c * t)
			// Time step: v(t + dt) = v0 * exp(-c * (t + dt)) = v0 * exp(-c * t) * exp(-c * dt) = v * exp(-c * dt)
			// v2 = exp(-c * dt) * v1
			// Pade approximation:
			// v2 = v1 * 1 / (1 + c * dt)
			v = s2MulSV(1.0f / (1.0f + h * body->linearDamping), v);
			w *= 1.0f / (1.0f + h * body->angularDamping);

			body->linearVelocity = v;
			body->angularVelocity = w;
		}
	}

	// Solver data
	s2ContactSolver contactSolver = {0};
	s2ContactSolver_Initialize(&contactSolver);

	int jointCapacity = world->jointPool.capacity;
	for (int i = 0; i < jointCapacity; ++i)
	{
		s2Joint* joint = joints + i;
		if (s2ObjectIsFree(&joint->object))
		{
			continue;
		}
		s2InitVelocityConstraints(joint, context);
	}

	// Solve velocity constraints
	for (int i = 0; i < context->velocityIterations; ++i)
	{
		for (int i = 0; i < jointCapacity; ++i)
		{
			s2Joint* joint = joints + i;
			if (s2ObjectIsFree(&joint->object))
			{
				continue;
			}

			s2SolveVelocityConstraints(joint, context);
		}

		s2ContactSolver_SolveVelocityConstraints(&contactSolver);
	}

	// Special handling for restitution
	s2ContactSolver_ApplyRestitution(&contactSolver);

	// Store impulses for warm starting
	s2ContactSolver_StoreImpulses(&contactSolver);

	// Integrate positions
	for (int i = 0; i < bodyCapacity; ++i)
	{
		s2Body* body = bodies + i;
		if (s2ObjectIsFree(&body->object))
		{
			continue;
		}

		s2Vec2 c = body->position;
		float a = body->angle;
		s2Vec2 v = body->linearVelocity;
		float w = body->angularVelocity;

		// Clamp large velocities
		s2Vec2 translation = s2MulSV(h, v);
		if (s2Dot(translation, translation) > s2_maxTranslation * s2_maxTranslation)
		{
			float ratio = s2_maxTranslation / s2Length(translation);
			v = s2MulSV(ratio, v);
		}

		float rotation = h * w;
		if (rotation * rotation > s2_maxRotation * s2_maxRotation)
		{
			float ratio = s2_maxRotation / S2_ABS(rotation);
			w *= ratio;
		}

		// Integrate
		c = s2MulAdd(c, h, v);
		a += h * w;

		body->position = c;
		body->angle = a;
		body->linearVelocity = v;
		body->angularVelocity = w;

		body->position0 = body->position;
		body->angle0 = body->angle;
	}

	// Solve position constraints
	bool positionSolved = false;
	for (int i = 0; i < context->positionIterations; ++i)
	{
		bool contactsOkay = s2ContactSolver_SolvePositionConstraintsBlock(&contactSolver);

		bool jointsOkay = true;
		for (int i = 0; i < jointCapacity; ++i)
		{
			s2Joint* joint = joints + i;
			if (s2ObjectIsFree(&joint->object))
			{
				continue;
			}

			bool jointOkay = s2SolvePositionConstraints(joint, context);
			jointsOkay = jointsOkay && jointOkay;
		}

		if (contactsOkay && jointsOkay)
		{
			// Exit early if the position errors are small.
			positionSolved = true;
			break;
		}
	}

	// Update transform
	for (int i = 0; i < bodyCapacity; ++i)
	{
		s2Body* body = bodies + i;
		if (s2ObjectIsFree(&body->object))
		{
			continue;
		}
		body->transform.q = s2MakeRot(body->angle);
		body->transform.p = s2Sub(body->position, s2RotateVector(body->transform.q, body->localCenter));
	}

	s2Contact* contacts = world->contacts;
	const s2Vec2 aabbMargin = {s2_aabbMargin, s2_aabbMargin};

	for (int i = 0; i < bodyCapacity; ++i)
	{
		s2Body* body = bodies + i;
		if (s2ObjectIsFree(&body->object))
		{
			continue;
		}

		body->force = s2Vec2_zero;
		body->torque = 0.0f;

		// Update shapes AABBs
		int shapeIndex = body->shapeList;
		while (shapeIndex != S2_NULL_INDEX)
		{
			s2Shape* shape = world->shapes + shapeIndex;

			shape->aabb = s2Shape_ComputeAABB(shape, body->transform);

			if (s2AABB_Contains(shape->fatAABB, shape->aabb) == false)
			{
				shape->fatAABB.lowerBound = s2Sub(shape->aabb.lowerBound, aabbMargin);
				shape->fatAABB.upperBound = s2Add(shape->aabb.upperBound, aabbMargin);
				shape->enlargedAABB = true;
			}

			shapeIndex = shape->nextShapeIndex;
		}
	}
}
