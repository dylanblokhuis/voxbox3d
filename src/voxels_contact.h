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

// The integer voxel range [low, high] a canonical cuboid spans (grows along filled axes).
typedef struct b3CanonicalVoxel
{
	b3IVec3 low;
	b3IVec3 high;
} b3CanonicalVoxel;

// Build the canonical voxel range for an exposed voxel (grows to the domain limit along each axis
// whose face is FILLED, not free). See voxels_contact.c for the rationale.
b3CanonicalVoxel b3CanonicalFromVoxel( const b3Voxels* v, b3IVec3 key, uint8_t freeFaces );

// Materialize the canonical cuboid: expanded axes clamped to `domain` (the other shape's loosened
// AABB in this voxel grid's local frame). Outputs local-space center and half extents.
void b3CanonicalCuboid( const b3Voxels* v, b3IVec3 key, b3CanonicalVoxel canon, b3AABB domain, b3Vec3* outCenter,
						b3Vec3* outHalf );

// Compute contact manifolds between a voxel shape (shapeA) and a convex shape (shapeB:
// sphere/capsule/hull), using the canonical-cuboid internal-edges fix. Mirrors the role of
// b3ComputeMeshManifolds. Returns true if touching.
bool b3ComputeVoxelsManifolds( b3World* world, int workerIndex, b3Contact* contact, const b3Shape* shapeA,
							   b3WorldTransform xfA, const b3Shape* shapeB, b3WorldTransform xfB, bool isFast, b3Arena arena );

// Compute contact manifolds between two voxel shapes. Returns true if touching.
bool b3ComputeVoxelsVoxelsManifolds( b3World* world, int workerIndex, b3Contact* contact, const b3Shape* shapeA,
									 b3WorldTransform xfA, const b3Shape* shapeB, b3WorldTransform xfB, bool isFast,
									 b3Arena arena );
