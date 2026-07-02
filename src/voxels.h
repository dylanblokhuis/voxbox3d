// SPDX-FileCopyrightText: 2025 Erin Catto
// SPDX-License-Identifier: MIT
//
// Voxel collision shape — dense fixed-grid port of parry's `Voxels`.
// See G:/parry/VOXEL_COLLISIONS_GUIDE.md. The runtime type `b3Voxels` is owned by the shape
// union (src/shape.h) as a heap pointer, mirroring how `b3CompoundData` is stored.

#pragma once

#include "voxels_consts.h"

#include "box3d/types.h"

#include <math.h>
#include <stdbool.h>

// A voxel's neighborhood state, encoded as the 6-neighbor byte described in voxels_consts.h.
typedef uint8_t b3VoxelState;

static inline bool b3VoxelStateIsEmpty( b3VoxelState s )
{
	return s == B3_VOXEL_EMPTY_FACE_MASK;
}

// Which faces have NO adjacent non-empty voxel (the exposed faces).
static inline uint8_t b3VoxelStateFreeFaces( b3VoxelState s )
{
	if ( s == B3_VOXEL_INTERIOR_FACE_MASK || s == B3_VOXEL_EMPTY_FACE_MASK )
	{
		return 0;
	}
	return ( ~s ) & b3_voxelAllFaces;
}

static inline b3VoxelType b3VoxelStateType( b3VoxelState s )
{
	B3_ASSERT( s < 65 );
	return (b3VoxelType)b3_facesToVoxelTypes[s];
}

static inline uint16_t b3VoxelStateFeatureMask( b3VoxelState s )
{
	B3_ASSERT( s < 65 );
	return b3_facesToFeatureMasks[s];
}

static inline uint32_t b3VoxelStateOctantMask( b3VoxelState s )
{
	B3_ASSERT( s < 65 );
	return b3_facesToOctantMasks[s];
}

// The integer domain of the grid: voxel keys span [0, dim) per axis.
typedef struct b3IVec3
{
	int x, y, z;
} b3IVec3;

// Dense fixed-grid voxel storage. State index = x + y * cx + z * cx * cy.
// Local-space: voxel (x,y,z) min corner = origin + (x,y,z) * voxelSize; center adds half voxelSize.
typedef struct b3Voxels
{
	b3VoxelState* states;	// cx*cy*cz entries
	int cx, cy, cz;			// grid dimensions (each >= 1)
	b3Vec3 voxelSize;		// per-axis edge length (each component > 0)
	b3Vec3 origin;			// local-space min corner of voxel (0,0,0)
	b3AABB localAABB;		// cached bounds of the whole grid
	b3Vec3 proxyCorners[8]; // 8 corners of localAABB, for b3MakeShapeProxy coarse queries
} b3Voxels;

// ---- Geometry helpers ----

static inline int b3VoxelsIndex( const b3Voxels* v, int x, int y, int z )
{
	return x + y * v->cx + z * v->cx * v->cy;
}

// Center of voxel with integer key k.
static inline b3Vec3 b3VoxelsCenterKey( const b3Voxels* v, b3IVec3 k )
{
	b3Vec3 c;
	c.x = v->origin.x + ( (float)k.x + 0.5f ) * v->voxelSize.x;
	c.y = v->origin.y + ( (float)k.y + 0.5f ) * v->voxelSize.y;
	c.z = v->origin.z + ( (float)k.z + 0.5f ) * v->voxelSize.z;
	return c;
}

static inline b3Vec3 b3VoxelsVoxelCenter( const b3Voxels* v, int x, int y, int z )
{
	return b3VoxelsCenterKey( v, (b3IVec3){ x, y, z } );
}

// Integer key of the voxel containing local point p (floor division, clamped to domain).
static inline b3IVec3 b3VoxelsAtPoint( const b3Voxels* v, b3Vec3 p )
{
	b3Vec3 q;
	q.x = ( p.x - v->origin.x ) / v->voxelSize.x;
	q.y = ( p.y - v->origin.y ) / v->voxelSize.y;
	q.z = ( p.z - v->origin.z ) / v->voxelSize.z;
	b3IVec3 k;
	k.x = (int)floorf( q.x );
	k.y = (int)floorf( q.y );
	k.z = (int)floorf( q.z );
	return k;
}

static inline bool b3VoxelsInBounds( const b3Voxels* v, int x, int y, int z )
{
	return x >= 0 && y >= 0 && z >= 0 && x < v->cx && y < v->cy && z < v->cz;
}

// Bounds-checked state lookup. Returns B3_VOXEL_EMPTY_FACE_MASK for out-of-domain keys.
static inline b3VoxelState b3VoxelsState( const b3Voxels* v, int x, int y, int z )
{
	if ( !b3VoxelsInBounds( v, x, y, z ) )
	{
		return B3_VOXEL_EMPTY_FACE_MASK;
	}
	return v->states[b3VoxelsIndex( v, x, y, z )];
}

// The integer voxel range [mins, maxs] (semi-open) intersecting a local AABB, NOT clamped.
typedef struct b3VoxelRange
{
	b3IVec3 mins;
	b3IVec3 maxs;
} b3VoxelRange;

static inline b3VoxelRange b3VoxelsRangeIntersectingAABB( const b3Voxels* v, b3AABB aabb )
{
	b3Vec3 lo, hi;
	lo.x = ( aabb.lowerBound.x - v->origin.x ) / v->voxelSize.x;
	lo.y = ( aabb.lowerBound.y - v->origin.y ) / v->voxelSize.y;
	lo.z = ( aabb.lowerBound.z - v->origin.z ) / v->voxelSize.z;
	hi.x = ( aabb.upperBound.x - v->origin.x ) / v->voxelSize.x;
	hi.y = ( aabb.upperBound.y - v->origin.y ) / v->voxelSize.y;
	hi.z = ( aabb.upperBound.z - v->origin.z ) / v->voxelSize.z;
	b3VoxelRange r;
	r.mins.x = (int)floorf( lo.x );
	r.mins.y = (int)floorf( lo.y );
	r.mins.z = (int)floorf( lo.z );
	r.maxs.x = (int)ceilf( hi.x );
	r.maxs.y = (int)ceilf( hi.y );
	r.maxs.z = (int)ceilf( hi.z );
	return r;
}

// ---- Lifecycle ----

// Build a voxel grid from a creation definition. Clones the states array. Returns NULL on bad input.
b3Voxels* b3CreateVoxels( const b3VoxelsDef* def );
void b3DestroyVoxels( b3Voxels* v );

// Recompute every voxel's neighborhood byte from the grid (called at the end of construction).
void b3VoxelsRecomputeAllStates( b3Voxels* v );

// Number of non-empty voxels.
int b3VoxelsFilledCount( const b3Voxels* v );

// Mass properties: each filled voxel is a uniform cuboid; combined via parallel-axis theorem.
b3MassData b3ComputeVoxelsMass( const b3Voxels* v, float density );

// Local-space centroid (unweighted mean of filled voxel centers; falls back to AABB center).
b3Vec3 b3VoxelsCentroid( const b3Voxels* v );

// Total exposed surface area (sum over free faces of each voxel).
float b3VoxelsSurfaceArea( const b3Voxels* v );

// Exposed surface area projected onto a plane with the given normal.
float b3VoxelsProjectedArea( const b3Voxels* v, b3Vec3 planeNormal );

// Flip a voxel filled/empty and update its 6 neighbors' bytes. Returns the previous state.
b3VoxelState b3VoxelsSetVoxel( b3Voxels* v, int x, int y, int z, bool filled );

// ---- World queries (voxels treated as plain sparse cubes) ----

// Ray cast against the grid (local space). Amanatides-Woo DDA march.
b3CastOutput b3RayCastVoxels( const b3Voxels* v, const b3RayCastInput* input );

// Shape cast a convex proxy (in the grid's local frame) against the grid. Returns the earliest hit.
b3CastOutput b3ShapeCastVoxels( const b3Voxels* v, const b3ShapeCastInput* input );

// Overlap test against a convex proxy already transformed into the grid's local frame.
bool b3OverlapVoxels( const b3Voxels* v, b3Transform transform, const b3ShapeProxy* proxy );

// Character mover collision. Emits outward planes for voxels near the mover capsule.
int b3CollideMoverAndVoxels( b3PlaneResult* planes, int capacity, const b3Voxels* v, const b3Capsule* mover );
