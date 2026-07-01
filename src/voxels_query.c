// SPDX-FileCopyrightText: 2026 Erin Catto
// SPDX-License-Identifier: MIT
//
// Voxel world queries: ray cast (DDA), overlap, and character mover. For these non-manifold
// queries voxels are treated purely as sparse axis-aligned cubes (guide section 11) — the
// feature/octant machinery is irrelevant here. Per-voxel work delegates to the existing box-hull
// query routines so we reuse well-tested convex code.

#include "voxels.h"

#include "aabb.h"
#include "core.h"
#include "math_internal.h"
#include "shape.h"

#include "box3d/collision.h"

#include <float.h>
#include <math.h>

// Build a box hull for voxel (x,y,z) at its center in local space.
static b3BoxHull b3VoxelBoxHull( const b3Voxels* v, int x, int y, int z )
{
	b3Vec3 c = b3VoxelsVoxelCenter( v, x, y, z );
	b3Transform xf = { c, b3Quat_identity };
	return b3MakeTransformedBoxHull( 0.5f * v->voxelSize.x, 0.5f * v->voxelSize.y, 0.5f * v->voxelSize.z, xf );
}

b3CastOutput b3RayCastVoxels( const b3Voxels* v, const b3RayCastInput* input )
{
	b3CastOutput output = { 0 };

	b3Vec3 p1 = input->origin;
	b3Vec3 d = input->translation;

	// Clip the ray to the whole-grid AABB to find the entry parameter.
	float tmin, tmax;
	b3Vec3 p2 = b3Add( p1, d );
	if ( !b3RayCastAABB( v->localAABB, p1, p2, &tmin, &tmax ) )
	{
		return output; // misses the grid entirely
	}

	float tEnter = tmin > 0.0f ? tmin : 0.0f;
	float maxFraction = input->maxFraction;
	if ( tEnter > maxFraction )
	{
		return output;
	}

	// Starting cell (clamped into the grid).
	b3Vec3 entryPoint = b3MulAdd( p1, tEnter, d );
	b3IVec3 cell = b3VoxelsAtPoint( v, entryPoint );
	cell.x = cell.x < 0 ? 0 : ( cell.x >= v->cx ? v->cx - 1 : cell.x );
	cell.y = cell.y < 0 ? 0 : ( cell.y >= v->cy ? v->cy - 1 : cell.y );
	cell.z = cell.z < 0 ? 0 : ( cell.z >= v->cz ? v->cz - 1 : cell.z );

	// Amanatides-Woo DDA march.
	for ( int guard = 0; guard < ( v->cx + v->cy + v->cz ) + 3; ++guard )
	{
		if ( !b3VoxelStateIsEmpty( b3VoxelsState( v, cell.x, cell.y, cell.z ) ) )
		{
			// Exact hit against this voxel's box.
			b3BoxHull box = b3VoxelBoxHull( v, cell.x, cell.y, cell.z );
			b3CastOutput hit = b3RayCastHull( &box.base, input );
			if ( hit.hit && hit.fraction <= maxFraction )
			{
				return hit;
			}
		}

		// Per-axis exit parameter from the current cell.
		b3Vec3 cmin, cmax;
		cmin.x = v->origin.x + (float)cell.x * v->voxelSize.x;
		cmin.y = v->origin.y + (float)cell.y * v->voxelSize.y;
		cmin.z = v->origin.z + (float)cell.z * v->voxelSize.z;
		cmax.x = cmin.x + v->voxelSize.x;
		cmax.y = cmin.y + v->voxelSize.y;
		cmax.z = cmin.z + v->voxelSize.z;

		float tExit[3];
		int step[3];
		for ( int i = 0; i < 3; ++i )
		{
			float di = b3GetByIndex( d, i );
			float mn = b3GetByIndex( cmin, i );
			float mx = b3GetByIndex( cmax, i );
			float oi = b3GetByIndex( p1, i );
			if ( di > 0.0f )
			{
				tExit[i] = ( mx - oi ) / di;
				step[i] = 1;
			}
			else if ( di < 0.0f )
			{
				tExit[i] = ( mn - oi ) / di;
				step[i] = -1;
			}
			else
			{
				tExit[i] = FLT_MAX;
				step[i] = 0;
			}
			if ( tExit[i] < 0.0f )
			{
				tExit[i] = FLT_MAX;
			}
		}

		// Advance along the axis with the nearest exit.
		int axis = 0;
		if ( tExit[1] < tExit[axis] ) axis = 1;
		if ( tExit[2] < tExit[axis] ) axis = 2;

		if ( tExit[axis] > maxFraction )
		{
			break; // ray leaves the grid before hitting anything
		}

		if ( axis == 0 )
			cell.x += step[0];
		else if ( axis == 1 )
			cell.y += step[1];
		else
			cell.z += step[2];

		if ( !b3VoxelsInBounds( v, cell.x, cell.y, cell.z ) )
		{
			break;
		}
	}

	return output;
}

bool b3OverlapVoxels( const b3Voxels* v, b3Transform transform, const b3ShapeProxy* proxy )
{
	// The proxy is in world space and `transform` maps voxel-local -> world (same convention as
	// b3OverlapHull). For the broad-phase range we need the proxy AABB in the grid's local frame.
	b3AABB worldAabb = b3ComputeProxyAABB( proxy );
	b3Transform invTransform = b3InvertTransform( transform );
	b3AABB qaabb = b3AABB_Transform( invTransform, worldAabb );
	b3VoxelRange range = b3VoxelsRangeIntersectingAABB( v, qaabb );
	int x0 = range.mins.x < 0 ? 0 : range.mins.x;
	int y0 = range.mins.y < 0 ? 0 : range.mins.y;
	int z0 = range.mins.z < 0 ? 0 : range.mins.z;
	int x1 = range.maxs.x > v->cx ? v->cx : range.maxs.x;
	int y1 = range.maxs.y > v->cy ? v->cy : range.maxs.y;
	int z1 = range.maxs.z > v->cz ? v->cz : range.maxs.z;

	for ( int z = z0; z < z1; ++z )
	{
		for ( int y = y0; y < y1; ++y )
		{
			for ( int x = x0; x < x1; ++x )
			{
				if ( b3VoxelStateIsEmpty( b3VoxelsState( v, x, y, z ) ) )
				{
					continue;
				}
				b3BoxHull box = b3VoxelBoxHull( v, x, y, z );
				if ( b3OverlapHull( &box.base, transform, proxy ) )
				{
					return true;
				}
			}
		}
	}
	return false;
}

int b3CollideMoverAndVoxels( b3PlaneResult* planes, int capacity, const b3Voxels* v, const b3Capsule* mover )
{
	if ( capacity == 0 )
	{
		return 0;
	}

	// AABB of the mover capsule (local space) expanded by its radius.
	b3Vec3 lo = b3Min( mover->center1, mover->center2 );
	b3Vec3 hi = b3Max( mover->center1, mover->center2 );
	b3Vec3 r = { mover->radius, mover->radius, mover->radius };
	b3AABB maabb = { b3Sub( lo, r ), b3Add( hi, r ) };

	b3VoxelRange range = b3VoxelsRangeIntersectingAABB( v, maabb );
	int x0 = range.mins.x < 0 ? 0 : range.mins.x;
	int y0 = range.mins.y < 0 ? 0 : range.mins.y;
	int z0 = range.mins.z < 0 ? 0 : range.mins.z;
	int x1 = range.maxs.x > v->cx ? v->cx : range.maxs.x;
	int y1 = range.maxs.y > v->cy ? v->cy : range.maxs.y;
	int z1 = range.maxs.z > v->cz ? v->cz : range.maxs.z;

	int count = 0;
	for ( int z = z0; z < z1 && count < capacity; ++z )
	{
		for ( int y = y0; y < y1 && count < capacity; ++y )
		{
			for ( int x = x0; x < x1 && count < capacity; ++x )
			{
				b3VoxelState s = b3VoxelsState( v, x, y, z );
				if ( b3VoxelStateIsEmpty( s ) )
				{
					continue;
				}
				// Skip interior voxels: they can't produce an exposed contact plane.
				if ( b3VoxelStateType( s ) == b3_voxelTypeInterior )
				{
					continue;
				}
				b3BoxHull box = b3VoxelBoxHull( v, x, y, z );
				int n = b3CollideMoverAndHull( planes + count, &box.base, mover );
				count += n;
			}
		}
	}
	return count;
}
