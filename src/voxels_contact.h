// SPDX-FileCopyrightText: 2026 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "arena_allocator.h"
#include "voxels.h"

#include "box3d/types.h"

#include <stdbool.h>

typedef struct b3World b3World;
typedef struct b3Contact b3Contact;
typedef struct b3Shape b3Shape;

// Compute contact manifolds between a voxel shape (shapeA) and a convex shape (shapeB:
// sphere/capsule/hull), using the canonical-cuboid internal-edges fix. Mirrors the role of
// b3ComputeMeshManifolds. Returns true if touching.
bool b3ComputeVoxelsManifolds( b3World* world, int workerIndex, b3Contact* contact, const b3Shape* shapeA,
							   b3WorldTransform xfA, const b3Shape* shapeB, b3WorldTransform xfB, bool isFast, b3Arena arena );
