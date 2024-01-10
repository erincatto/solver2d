// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver2d/id.h"
#include "solver2d/joint_types.h"
#include "solver2d/timer.h"
#include "solver2d/types.h"

typedef struct Capsule Capsule;
typedef struct Circle Circle;
typedef struct Polygon Polygon;
typedef struct DebugDraw DebugDraw;
typedef struct Segment Segment;

#ifdef __cplusplus
extern "C"
{
#endif

s2WorldId s2CreateWorld(const s2WorldDef* def);
void s2DestroyWorld(s2WorldId worldId);

void s2World_Step(s2WorldId worldId, float timeStep, int32_t velocityIterations, int32_t positionIterations);
void s2World_Draw(s2WorldId worldId, s2DebugDraw* debugDraw);

struct s2Statistics s2World_GetStatistics(s2WorldId worldId);

s2BodyId s2CreateBody(s2WorldId worldId, const s2BodyDef* def);

void s2DestroyBody(s2BodyId bodyId);

s2Vec2 s2Body_GetPosition(s2BodyId bodyId);
float s2Body_GetAngle(s2BodyId bodyId);
s2Vec2 s2Body_GetLocalPoint(s2BodyId bodyId, s2Vec2 globalPoint);

void s2Body_SetTransform(s2BodyId bodyId, s2Vec2 position, float angle);
void s2Body_SetLinearVelocity(s2BodyId bodyId, s2Vec2 linearVelocity);
void s2Body_SetAngularVelocity(s2BodyId bodyId, float angularVelocity);

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
