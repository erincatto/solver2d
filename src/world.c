// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#define _CRT_SECURE_NO_WARNINGS

#include "world.h"

#include "allocate.h"
#include "array.h"
#include "block_allocator.h"
#include "body.h"
#include "broad_phase.h"
#include "contact.h"
#include "core.h"
#include "joint.h"
#include "pool.h"
#include "shape.h"
#include "solvers.h"
#include "stack_allocator.h"

#include "solver2d/aabb.h"
#include "solver2d/constants.h"
#include "solver2d/debug_draw.h"
#include "solver2d/distance.h"
#include "solver2d/solver2d.h"
#include "solver2d/timer.h"

#include <stdio.h>
#include <string.h>

s2World s2_worlds[s2_maxWorlds];

s2World* s2GetWorldFromId(s2WorldId id)
{
	S2_ASSERT(0 <= id.index && id.index < s2_maxWorlds);
	s2World* world = s2_worlds + id.index;
	S2_ASSERT(id.revision == world->revision);
	return world;
}

s2World* s2GetWorldFromIndex(int16_t index)
{
	S2_ASSERT(0 <= index && index < s2_maxWorlds);
	s2World* world = s2_worlds + index;
	S2_ASSERT(world->blockAllocator != NULL);
	return world;
}

s2WorldId s2CreateWorld(const s2WorldDef* def)
{
	s2WorldId id = s2_nullWorldId;
	for (int16_t i = 0; i < s2_maxWorlds; ++i)
	{
		if (s2_worlds[i].blockAllocator == NULL)
		{
			id.index = i;
			break;
		}
	}

	if (id.index == s2_nullWorldId.index)
	{
		return id;
	}

	s2InitializeContactRegisters();

	s2World empty = {0};
	s2World* world = s2_worlds + id.index;
	*world = empty;

	world->index = id.index;

	world->blockAllocator = s2CreateBlockAllocator();
	world->stackAllocator = s2CreateStackAllocator(1024 * 1024);

	world->solverType = def->solverType;

	s2CreateBroadPhase(&world->broadPhase);

	// pools
	world->bodyPool = s2CreatePool(sizeof(s2Body), S2_MAX(4, 1));
	world->bodies = (s2Body*)world->bodyPool.memory;

	world->shapePool = s2CreatePool(sizeof(s2Shape), S2_MAX(4, 1));
	world->shapes = (s2Shape*)world->shapePool.memory;

	world->contactPool = s2CreatePool(sizeof(s2Contact), S2_MAX(4, 1));
	world->contacts = (s2Contact*)world->contactPool.memory;

	world->jointPool = s2CreatePool(sizeof(s2Joint), S2_MAX(4, 1));
	world->joints = (s2Joint*)world->jointPool.memory;

	world->stepId = 0;

	// Globals start at 0. It should be fine for this to roll over.
	world->revision += 1;

	world->gravity = (s2Vec2){0.0f, -10.0f};
	world->restitutionThreshold = 1.0f;

	id.revision = world->revision;

	return id;
}

void s2DestroyWorld(s2WorldId id)
{
	s2World* world = s2GetWorldFromId(id);

	s2DestroyPool(&world->jointPool);
	s2DestroyPool(&world->contactPool);
	s2DestroyPool(&world->shapePool);
	s2DestroyPool(&world->bodyPool);
	s2DestroyBroadPhase(&world->broadPhase);
	s2DestroyBlockAllocator(world->blockAllocator);
	s2DestroyStackAllocator(world->stackAllocator);

	memset(world, 0, sizeof(s2World));
}

void s2World_Step(s2WorldId worldId, float timeStep, int velocityIterations, int positionIterations)
{
	s2World* world = s2GetWorldFromId(worldId);
	world->stepId += 1;

	// Stage 1: Update collision pairs and create contacts
	s2UpdateBroadPhasePairs(world);

	// Stage 2: Optimize broad-phase
	s2BroadPhase* broadPhase = &world->broadPhase;
	s2BroadPhase_RebuildTrees(broadPhase);

	// Stage 3: Update contacts
	s2Shape* shapes = world->shapes;
	s2Body* bodies = world->bodies;
	s2Contact* contacts = world->contacts;
	int contactCapacity = world->contactPool.capacity;

	for (int i = 0; i < contactCapacity; ++i)
	{
		s2Contact* contact = world->contacts + i;
		if (s2IsFree(&contact->object))
		{
			continue;
		}

		s2Shape* shapeA = shapes + contact->shapeIndexA;
		s2Shape* shapeB = shapes + contact->shapeIndexB;

		// Do proxies still overlap?
		bool overlap = s2AABB_Overlaps(shapeA->fatAABB, shapeB->fatAABB);
		if (overlap == true)
		{
			// Update contact respecting shape/body order (A,B)
			s2Body* bodyA = bodies + shapeA->bodyIndex;
			s2Body* bodyB = bodies + shapeB->bodyIndex;
			s2UpdateContact(world, contact, shapeA, bodyA, shapeB, bodyB);
		}
		else
		{
			// this does not change the contact capacity or invalidate this loop
			s2DestroyContact(world, contact);
		}
	}

	// Stage 3: Integrate velocities, solve velocity constraints, and integrate positions.
	s2StepContext context = {0};
	context.dt = timeStep;
	context.velocityIterations = velocityIterations;
	context.positionIterations = positionIterations;
	if (timeStep > 0.0f)
	{
		context.inv_dt = 1.0f / timeStep;
	}
	else
	{
		context.inv_dt = 0.0f;
	}

	context.restitutionThreshold = world->restitutionThreshold;
	context.bodies = world->bodies;
	context.bodyCapacity = world->bodyPool.capacity;

	if (timeStep > 0.0f)
	{
		switch (world->solverType)
		{
			case s2_solverPGS_NGS_Block:
				s2Solve_PGS_NGS_Block(world, &context);
				break;

			case s2_solverPGS_NGS:
				s2Solve_PGS_NGS(world, &context);
				break;

			case s2_solverPGS_Soft:
				s2Solve_PGS_Soft(world, &context);
				break;

			case s2_solverXPBD:
				s2Solve_XPDB(world, &context);
				break;

			case s2_solverTGS_Soft:
				s2Solve_TGS_Soft(world, &context);
				break;

			case s2_solverTGS_Sticky:
				s2Solve_TGS_Sticky(world, &context);
				break;

			case s2_solverTGS_NGS:
				s2Solve_TGS_NGS(world, &context);
				break;

			default:
				break;
		}
	}

	// Stage 4: Update transforms and broad-phase
	int bodyCapacity = world->bodyPool.capacity;
	const s2Vec2 aabbMargin = {s2_aabbMargin, s2_aabbMargin};
	for (int i = 0; i < bodyCapacity; ++i)
	{
		s2Body* body = bodies + i;
		if (s2IsFree(&body->object))
		{
			continue;
		}

		if (body->type == s2_staticBody)
		{
			continue;
		}

		body->transform.q = s2MakeRot(body->angle);
		body->transform.p = s2Sub(body->position, s2RotateVector(body->transform.q, body->localCenter));
		body->force = s2Vec2_zero;
		body->torque = 0.0f;

		int shapeIndex = body->shapeList;
		while (shapeIndex != S2_NULL_INDEX)
		{
			s2Shape* shape = world->shapes + shapeIndex;

			shape->aabb = s2Shape_ComputeAABB(shape, body->transform);

			if (s2AABB_Contains(shape->fatAABB, shape->aabb) == false)
			{
				shape->fatAABB.lowerBound = s2Sub(shape->aabb.lowerBound, aabbMargin);
				shape->fatAABB.upperBound = s2Add(shape->aabb.upperBound, aabbMargin);
				s2BroadPhase_EnlargeProxy(broadPhase, shape->proxyKey, shape->fatAABB);
			}

			shapeIndex = shape->nextShapeIndex;
		}
	}

	s2ValidateBroadphase(&world->broadPhase);

	s2GrowStack(world->stackAllocator);
}

static void s2DrawShape(s2DebugDraw* draw, s2Shape* shape, s2Transform xf, s2Color color)
{
	switch (shape->type)
	{
		case s2_capsuleShape:
		{
			s2Capsule* capsule = &shape->capsule;
			s2Vec2 p1 = s2TransformPoint(xf, capsule->point1);
			s2Vec2 p2 = s2TransformPoint(xf, capsule->point2);
			draw->DrawSolidCapsule(p1, p2, capsule->radius, color, draw->context);
		}
		break;

		case s2_circleShape:
		{
			s2Circle* circle = &shape->circle;
			s2Vec2 center = s2TransformPoint(xf, circle->point);
			s2Vec2 axis = s2RotateVector(xf.q, (s2Vec2){1.0f, 0.0f});
			draw->DrawSolidCircle(center, circle->radius, axis, color, draw->context);
		}
		break;

		case s2_polygonShape:
		{
			s2Color fillColor = {0.5f * color.r, 0.5f * color.g, 0.5f * color.b, 0.5f};

			s2Polygon* poly = &shape->polygon;
			int count = poly->count;
			S2_ASSERT(count <= s2_maxPolygonVertices);
			s2Vec2 vertices[s2_maxPolygonVertices];

			for (int i = 0; i < count; ++i)
			{
				vertices[i] = s2TransformPoint(xf, poly->vertices[i]);
			}

			if (poly->radius > 0.0f)
			{
				draw->DrawRoundedPolygon(vertices, count, poly->radius, fillColor, color, draw->context);
			}
			else
			{
				draw->DrawSolidPolygon(vertices, count, color, draw->context);
			}
		}
		break;

		case s2_segmentShape:
		{
			s2Segment* segment = &shape->segment;
			s2Vec2 p1 = s2TransformPoint(xf, segment->point1);
			s2Vec2 p2 = s2TransformPoint(xf, segment->point2);
			draw->DrawSegment(p1, p2, color, draw->context);
		}
		break;

		default:
			break;
	}
}

void s2World_Draw(s2WorldId worldId, s2DebugDraw* draw)
{
	s2World* world = s2GetWorldFromId(worldId);

	if (draw->drawShapes)
	{
		int count = world->bodyPool.capacity;
		for (int i = 0; i < count; ++i)
		{
			s2Body* b = world->bodies + i;
			if (b->object.next != i)
			{
				continue;
			}

			s2Transform xf = b->transform;
			int shapeIndex = b->shapeList;
			while (shapeIndex != S2_NULL_INDEX)
			{
				s2Shape* shape = world->shapes + shapeIndex;
				if (b->type == s2_dynamicBody && b->mass == 0.0f)
				{
					// Bad body
					s2DrawShape(draw, shape, xf, (s2Color){1.0f, 0.0f, 0.0f, 1.0f});
				}
				else if (b->type == s2_staticBody)
				{
					s2DrawShape(draw, shape, xf, (s2Color){0.5f, 0.9f, 0.5f, 1.0f});
				}
				else if (b->type == s2_kinematicBody)
				{
					s2DrawShape(draw, shape, xf, (s2Color){0.5f, 0.5f, 0.9f, 1.0f});
				}
				else
				{
					s2DrawShape(draw, shape, xf, draw->dynamicBodyColor);
				}

				shapeIndex = shape->nextShapeIndex;
			}
		}
	}

	if (draw->drawJoints)
	{
		int count = world->jointPool.capacity;
		for (int i = 0; i < count; ++i)
		{
			s2Joint* joint = world->joints + i;
			if (joint->object.next != i)
			{
				continue;
			}

			s2DrawJoint(draw, world, joint);
		}
	}

	// if (debugDraw->drawPi & s2Draw::e_pairBit)
	//{
	//		s2Color color(0.3f, 0.9f, 0.9f);
	//		for (s2Contact* c = m_contactManager.m_contactList; c; c = c->GetNext())
	//		{
	//		s2Shape* fixtureA = c->GetFixtureA();
	//		s2Shape* fixtureB = c->GetFixtureB();
	//		int32 indexA = c->GetChildIndexA();
	//		int32 indexB = c->GetChildIndexB();
	//		s2Vec2 cA = fixtureA->GetAABB(indexA).GetCenter();
	//		s2Vec2 cB = fixtureB->GetAABB(indexB).GetCenter();

	//		m_debugDraw->DrawSegment(cA, cB, color);
	//		}
	//}

	if (draw->drawAABBs)
	{
		s2Color color = {0.9f, 0.3f, 0.9f, 1.0f};

		int count = world->bodyPool.capacity;
		for (int i = 0; i < count; ++i)
		{
			s2Body* b = world->bodies + i;
			if (b->object.next != i)
			{
				continue;
			}

			char buffer[32];
			sprintf(buffer, "%d", b->object.index);
			draw->DrawString(b->position, buffer, draw->context);

			int shapeIndex = b->shapeList;
			while (shapeIndex != S2_NULL_INDEX)
			{
				s2Shape* shape = world->shapes + shapeIndex;
				s2Box aabb = shape->fatAABB;

				s2Vec2 vs[4] = {{aabb.lowerBound.x, aabb.lowerBound.y},
								{aabb.upperBound.x, aabb.lowerBound.y},
								{aabb.upperBound.x, aabb.upperBound.y},
								{aabb.lowerBound.x, aabb.upperBound.y}};

				draw->DrawPolygon(vs, 4, color, draw->context);

				shapeIndex = shape->nextShapeIndex;
			}
		}

		// for (s2Shape* f = b->GetFixtureList(); f; f = f->GetNext())
		//{
		//	for (int32 i = 0; i < f->m_proxyCount; ++i)
		//	{
		//		s2FixtureProxy* proxy = f->m_proxies + i;
		//		s2Box aabb = bp->GetFatAABB(proxy->proxyId);
		//		s2Vec2 vs[4];
		//		vs[0].Set(aabb.lowerBound.x, aabb.lowerBound.y);
		//		vs[1].Set(aabb.upperBound.x, aabb.lowerBound.y);
		//		vs[2].Set(aabb.upperBound.x, aabb.upperBound.y);
		//		vs[3].Set(aabb.lowerBound.x, aabb.upperBound.y);

		//		m_debugDraw->DrawPolygon(vs, 4, color);
		//	}
		//}
	}

	// if (flags & s2Draw::e_centerOfMassBit)
	//{
	//		for (s2Body* b = m_bodyList; b; b = b->GetNext())
	//		{
	//		s2Transform xf = b->GetTransform();
	//		xf.p = b->GetWorldCenter();
	//		m_debugDraw->DrawTransform(xf);
	//		}
	// }
}

s2Statistics s2World_GetStatistics(s2WorldId worldId)
{
	s2World* world = s2GetWorldFromId(worldId);
	s2Statistics s = {0};
	s.bodyCount = world->bodyPool.count;
	s.contactCount = world->contactPool.count;
	s.jointCount = world->jointPool.count;

	s2DynamicTree* tree = world->broadPhase.trees + s2_dynamicBody;
	s.proxyCount = tree->nodeCount;
	s.treeHeight = s2DynamicTree_GetHeight(tree);
	s.stackCapacity = s2GetStackCapacity(world->stackAllocator);
	s.stackUsed = s2GetMaxStackAllocation(world->stackAllocator);
	return s;
}

#if 0
struct s2WorldRayCastWrapper
{
	float RayCastCallback(const s2RayCastInput& input, int32 proxyId)
	{
		void* userData = broadPhase->GetUserData(proxyId);
		s2FixtureProxy* proxy = (s2FixtureProxy*)userData;
		s2Shape* shape = proxy->shape;
		int32 index = proxy->childIndex;
		s2RayCastOutput output;
		bool hit = shape->RayCast(&output, input, index);

		if (hit)
		{
			float fraction = output.fraction;
			s2Vec2 point = (1.0f - fraction) * input.p1 + fraction * input.p2;
			return callback->ReportFixture(shape, point, output.normal, fraction);
		}

		return input.maxFraction;
	}

	const s2BroadPhase* broadPhase;
	s2RayCastCallback* callback;
};

void s2World::RayCast(s2RayCastCallback* callback, const s2Vec2& point1, const s2Vec2& point2) const
{
	s2WorldRayCastWrapper wrapper;
	wrapper.broadPhase = &m_contactManager.m_broadPhase;
	wrapper.callback = callback;
	s2RayCastInput input;
	input.maxFraction = 1.0f;
	input.p1 = point1;
	input.p2 = point2;
	m_contactManager.m_broadPhase.RayCast(&wrapper, input);
}

int32 s2World::GetProxyCount() const
{
	return m_contactManager.m_broadPhase.GetProxyCount();
}

int32 s2World::GetTreeHeight() const
{
	return m_contactManager.m_broadPhase.GetTreeHeight();
}

int32 s2World::GetTreeBalance() const
{
	return m_contactManager.m_broadPhase.GetTreeBalance();
}

float s2World::GetTreeQuality() const
{
	return m_contactManager.m_broadPhase.GetTreeQuality();
}

void s2World::ShiftOrigin(const s2Vec2& newOrigin)
{
	S2_ASSERT(m_locked == false);
	if (m_locked)
	{
		return;
	}

	for (s2Body* b = m_bodyList; b; b = b->m_next)
	{
		b->m_xf.p -= newOrigin;
		b->m_sweep.c0 -= newOrigin;
		b->m_sweep.c -= newOrigin;
	}

	for (s2Joint* j = m_jointList; j; j = j->m_next)
	{
		j->ShiftOrigin(newOrigin);
	}

	m_contactManager.m_broadPhase.ShiftOrigin(newOrigin);
}

void s2World::Dump()
{
	if (m_locked)
	{
		return;
	}

	s2OpenDump("box2d_dump.inl");

	s2Dump("s2Vec2 g(%.9g, %.9g);\n", m_gravity.x, m_gravity.y);
	s2Dump("m_world->SetGravity(g);\n");

	s2Dump("s2Body** bodies = (s2Body**)s2Alloc(%d * sizeof(s2Body*));\n", m_bodyCount);
	s2Dump("s2Joint** joints = (s2Joint**)s2Alloc(%d * sizeof(s2Joint*));\n", m_jointCount);

	int32 i = 0;
	for (s2Body* b = m_bodyList; b; b = b->m_next)
	{
		b->m_islandIndex = i;
		b->Dump();
		++i;
	}

	i = 0;
	for (s2Joint* j = m_jointList; j; j = j->m_next)
	{
		j->m_index = i;
		++i;
	}

	// First pass on joints, skip gear joints.
	for (s2Joint* j = m_jointList; j; j = j->m_next)
	{
		if (j->m_type == e_gearJoint)
		{
			continue;
		}

		s2Dump("{\n");
		j->Dump();
		s2Dump("}\n");
	}

	// Second pass on joints, only gear joints.
	for (s2Joint* j = m_jointList; j; j = j->m_next)
	{
		if (j->m_type != e_gearJoint)
		{
			continue;
		}

		s2Dump("{\n");
		j->Dump();
		s2Dump("}\n");
	}

	s2Dump("s2Free(joints);\n");
	s2Dump("s2Free(bodies);\n");
	s2Dump("joints = nullptr;\n");
	s2Dump("bodies = nullptr;\n");

	s2CloseDump();
}
#endif

typedef struct WorldQueryContext
{
	s2World* world;
	s2QueryCallbackFcn* fcn;
	void* userContext;
} WorldQueryContext;

static bool TreeQueryCallback(int proxyId, int shapeIndex, void* context)
{
	S2_MAYBE_UNUSED(proxyId);

	WorldQueryContext* worldContext = (WorldQueryContext*)context;
	s2World* world = worldContext->world;

	S2_ASSERT(0 <= shapeIndex && shapeIndex < world->shapePool.capacity);

	s2Shape* shape = world->shapes + shapeIndex;
	S2_ASSERT(shape->object.index == shape->object.next);

	s2ShapeId shapeId = {shape->object.index, world->index, shape->object.revision};
	bool result = worldContext->fcn(shapeId, worldContext->userContext);
	return result;
}

void s2World_QueryAABB(s2WorldId worldId, s2Box aabb, s2QueryCallbackFcn* fcn, void* context)
{
	s2World* world = s2GetWorldFromId(worldId);

	WorldQueryContext worldContext = {world, fcn, context};

	for (int i = 0; i < s2_bodyTypeCount; ++i)
	{
		s2DynamicTree_Query(world->broadPhase.trees + i, aabb, TreeQueryCallback, &worldContext);
	}
}

bool s2IsBodyIdValid(s2World* world, s2BodyId id)
{
	if (id.world != world->index)
	{
		return false;
	}

	if (id.index >= world->bodyPool.capacity)
	{
		return false;
	}

	s2Body* body = world->bodies + id.index;
	if (body->object.index != body->object.next)
	{
		return false;
	}

	if (body->object.revision != id.revision)
	{
		return false;
	}

	return true;
}
