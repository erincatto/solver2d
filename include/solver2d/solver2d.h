// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver2d/id.h"
#include "solver2d/joint_types.h"
#include "solver2d/timer.h"
#include "solver2d/types.h"

typedef struct s2Capsule s2Capsule;
typedef struct s2Circle s2Circle;
typedef struct s2Polygon s2Polygon;
typedef struct s2DebugDraw s2DebugDraw;
typedef struct s2Segment s2Segment;

#ifdef __cplusplus
extern "C"
{
#endif

s2WorldId s2CreateWorld(const s2WorldDef* def);
void s2DestroyWorld(s2WorldId worldId);

void s2World_Step(s2WorldId worldId, float timeStep, int32_t velIters, int32_t posIters, bool warmStart);
void s2World_Draw(s2WorldId worldId, s2DebugDraw* debugDraw);

/// Joint stiffness for the Soft Step solvers. Soft Step always caps it at a quarter of the substep rate.
/// The default is 60 Hz, or 480 Hz with a cap fraction of 2 for Soft Step D.
/// @param hertz joint stiffness in cycles per second
/// @param capFraction Soft Step D caps the stiffness at this fraction of the substep rate, zero for no cap
void s2World_SetJointHertz(s2WorldId worldId, float hertz, float capFraction);

/// Soft Step D solves its direct joints with conjugate residual instead of the exact factorization.
/// @param iterations conjugate residual iterations per solve, zero for the exact factorization (default)
void s2World_SetJointCRIterations(s2WorldId worldId, int iterations);

/// When Soft Step Seed runs its exact joint solve. The default is before the first substep.
void s2World_SetJointSeedMode(s2WorldId worldId, s2SeedMode mode);

/// Soft Step D relaxes its direct joints with the direct solve (default), or with Soft Step's Gauss-Seidel relax.
void s2World_EnableDirectRelax(s2WorldId worldId, bool flag);

struct s2Statistics s2World_GetStatistics(s2WorldId worldId);

s2BodyId s2CreateBody(s2WorldId worldId, const s2BodyDef* def);

void s2DestroyBody(s2BodyId bodyId);

s2Vec2 s2Body_GetPosition(s2BodyId bodyId);
float s2Body_GetAngle(s2BodyId bodyId);
s2Vec2 s2Body_GetLocalPoint(s2BodyId bodyId, s2Vec2 globalPoint);

void s2Body_SetTransform(s2BodyId bodyId, s2Vec2 position, float angle);
void s2Body_SetLinearVelocity(s2BodyId bodyId, s2Vec2 linearVelocity);
void s2Body_SetAngularVelocity(s2BodyId bodyId, float angularVelocity);
void s2Body_ApplyForceToCenter(s2BodyId bodyId, s2Vec2 force);
void s2Body_ApplyLinearImpulse(s2BodyId bodyId, s2Vec2 impulse, s2Vec2 point);

s2BodyType s2Body_GetType(s2BodyId bodyId);
float s2Body_GetMass(s2BodyId bodyId);

s2ShapeId s2CreateCircleShape(s2BodyId bodyId, const s2ShapeDef* def, const s2Circle* circle);
s2ShapeId s2CreateSegmentShape(s2BodyId bodyId, const s2ShapeDef* def, const s2Segment* segment);
s2ShapeId s2CreateCapsuleShape(s2BodyId bodyId, const s2ShapeDef* def, const s2Capsule* capsule);
s2ShapeId s2CreatePolygonShape(s2BodyId bodyId, const s2ShapeDef* def, const s2Polygon* polygon);

s2BodyId s2Shape_GetBody(s2ShapeId shapeId);
bool s2Shape_TestPoint(s2ShapeId shapeId, s2Vec2 point);

s2JointId s2CreateMouseJoint(s2WorldId worldId, const s2MouseJointDef* def);
s2JointId s2CreateRevoluteJoint(s2WorldId worldId, const s2RevoluteJointDef* def);
s2JointId s2CreateWeldJoint(s2WorldId worldId, const s2WeldJointDef* def);
s2JointId s2CreatePrismaticJoint(s2WorldId worldId, const s2PrismaticJointDef* def);
void s2DestroyJoint(s2JointId jointId);

void s2MouseJoint_SetTarget(s2JointId jointId, s2Vec2 target);

void s2RevoluteJoint_EnableLimit(s2JointId jointId, bool enableLimit);
void s2RevoluteJoint_EnableMotor(s2JointId jointId, bool enableMotor);
void s2RevoluteJoint_SetMotorSpeed(s2JointId jointId, float motorSpeed);
float s2RevoluteJoint_GetMotorTorque(s2JointId jointId, float inverseTimeStep);

/// This function receives shapes found in the AABB query.
/// @return true if the query should continue
typedef bool s2QueryCallbackFcn(s2ShapeId shapeId, void* context);

void s2World_QueryAABB(s2WorldId worldId, s2Box aabb, s2QueryCallbackFcn* fcn, void* context);

#ifdef __cplusplus
}
#endif
