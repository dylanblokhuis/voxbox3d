// SPDX-FileCopyrightText: 2025 Erin Catto
// SPDX-License-Identifier: MIT
//
// Voxel storage and neighborhood bookkeeping — port of parry's
// `voxels_edition.rs` / `voxels_neighborhood.rs` for a dense fixed grid.
// See G:/parry/VOXEL_COLLISIONS_GUIDE.md sections 3 and 5.

#include "voxels.h"

#include "core.h"
#include "math_internal.h"

#include "box3d/base.h"

#include <math.h>

// Read-only: sample the 6 axis-neighbors of (x,y,z) and build that voxel's neighborhood byte.
// Bit k*2 = +1 neighbor filled, bit k*2+1 = -1 neighbor filled (k = axis 0,1,2 = x,y,z).
static uint8_t b3VoxelsComputeNeighborhoodBits( const b3Voxels* v, int x, int y, int z )
{
	uint8_t bits = 0;
	const int coord[3] = { x, y, z };
	for ( int k = 0; k < 3; ++k )
	{
		int c = coord[k];
		// +1 neighbor
		int xp[3] = { x, y, z };
		xp[k] = c + 1;
		if ( !b3VoxelStateIsEmpty( b3VoxelsState( v, xp[0], xp[1], xp[2] ) ) )
		{
			bits |= (uint8_t)( 1 << ( k * 2 ) );
		}
		// -1 neighbor
		int xm[3] = { x, y, z };
		xm[k] = c - 1;
		if ( !b3VoxelStateIsEmpty( b3VoxelsState( v, xm[0], xm[1], xm[2] ) ) )
		{
			bits |= (uint8_t)( 1 << ( k * 2 + 1 ) );
		}
	}
	return bits;
}

void b3VoxelsRecomputeAllStates( b3Voxels* v )
{
	for ( int z = 0; z < v->cz; ++z )
	{
		for ( int y = 0; y < v->cy; ++y )
		{
			for ( int x = 0; x < v->cx; ++x )
			{
				int idx = b3VoxelsIndex( v, x, y, z );
				if ( b3VoxelStateIsEmpty( v->states[idx] ) )
				{
					continue; // empties keep the EMPTY sentinel
				}
				v->states[idx] = b3VoxelsComputeNeighborhoodBits( v, x, y, z );
			}
		}
	}
}

// Mutate the 6 axis-neighbors of (x,y,z) to reflect that it just became filled/empty, and return
// the new neighborhood byte for (x,y,z) itself. Port of parry update_neighbors_state.
// centerEmpty: did the center voxel just become empty?
static uint8_t b3VoxelsUpdateNeighborsState( b3Voxels* v, int x, int y, int z, bool centerEmpty )
{
	uint8_t keyData = 0;
	const int coord[3] = { x, y, z };
	for ( int k = 0; k < 3; ++k )
	{
		int c = coord[k];

		// Left neighbor at -1 along axis k. It sees `key` on its +1 side (bit k*2).
		int xl[3] = { x, y, z };
		xl[k] = c - 1;
		if ( b3VoxelsInBounds( v, xl[0], xl[1], xl[2] ) )
		{
			uint8_t* ls = &v->states[b3VoxelsIndex( v, xl[0], xl[1], xl[2] )];
			if ( !b3VoxelStateIsEmpty( *ls ) )
			{
				if ( centerEmpty )
				{
					*ls &= (uint8_t)~( 1 << ( k * 2 ) );
				}
				else
				{
					*ls |= (uint8_t)( 1 << ( k * 2 ) );
					keyData |= (uint8_t)( 1 << ( k * 2 + 1 ) ); // our -1 side is filled
				}
			}
		}

		// Right neighbor at +1 along axis k. It sees `key` on its -1 side (bit k*2+1).
		int xr[3] = { x, y, z };
		xr[k] = c + 1;
		if ( b3VoxelsInBounds( v, xr[0], xr[1], xr[2] ) )
		{
			uint8_t* rs = &v->states[b3VoxelsIndex( v, xr[0], xr[1], xr[2] )];
			if ( !b3VoxelStateIsEmpty( *rs ) )
			{
				if ( centerEmpty )
				{
					*rs &= (uint8_t)~( 1 << ( k * 2 + 1 ) );
				}
				else
				{
					*rs |= (uint8_t)( 1 << ( k * 2 + 1 ) );
					keyData |= (uint8_t)( 1 << ( k * 2 ) ); // our +1 side is filled
				}
			}
		}
	}
	return keyData;
}

b3VoxelState b3VoxelsSetVoxel( b3Voxels* v, int x, int y, int z, bool filled )
{
	if ( !b3VoxelsInBounds( v, x, y, z ) )
	{
		return B3_VOXEL_EMPTY_FACE_MASK;
	}
	int idx = b3VoxelsIndex( v, x, y, z );
	b3VoxelState prev = v->states[idx];
	bool wasEmpty = b3VoxelStateIsEmpty( prev );
	if ( wasEmpty == !filled )
	{
		return prev; // no flip
	}

	uint8_t newData = b3VoxelsUpdateNeighborsState( v, x, y, z, !filled );
	v->states[idx] = filled ? newData : B3_VOXEL_EMPTY_FACE_MASK;
	return prev;
}

b3Voxels* b3CreateVoxels( const b3VoxelsDef* def )
{
	if ( def == NULL || def->cx < 1 || def->cy < 1 || def->cz < 1 || def->voxelSize.x <= 0.0f ||
		 def->voxelSize.y <= 0.0f || def->voxelSize.z <= 0.0f )
	{
		B3_ASSERT( false );
		return NULL;
	}

	b3Voxels* v = (b3Voxels*)b3Alloc( sizeof( b3Voxels ) );
	if ( v == NULL )
	{
		return NULL;
	}

	int count = def->cx * def->cy * def->cz;
	v->cx = def->cx;
	v->cy = def->cy;
	v->cz = def->cz;
	v->voxelSize = def->voxelSize;
	v->origin = def->origin;
	v->states = (b3VoxelState*)b3Alloc( (size_t)count * sizeof( b3VoxelState ) );
	if ( v->states == NULL )
	{
		b3Free( v, sizeof( b3Voxels ) );
		return NULL;
	}

	for ( int i = 0; i < count; ++i )
	{
		v->states[i] = ( def->occupancy != NULL && def->occupancy[i] != 0 ) ? 0 : B3_VOXEL_EMPTY_FACE_MASK;
	}

	b3VoxelsRecomputeAllStates( v );

	// Cached local AABB = [origin, origin + voxelSize * dim].
	v->localAABB.lowerBound.x = v->origin.x;
	v->localAABB.lowerBound.y = v->origin.y;
	v->localAABB.lowerBound.z = v->origin.z;
	v->localAABB.upperBound.x = v->origin.x + (float)v->cx * v->voxelSize.x;
	v->localAABB.upperBound.y = v->origin.y + (float)v->cy * v->voxelSize.y;
	v->localAABB.upperBound.z = v->origin.z + (float)v->cz * v->voxelSize.z;

	return v;
}

int b3VoxelsFilledCount( const b3Voxels* v )
{
	int count = v->cx * v->cy * v->cz;
	int filled = 0;
	for ( int i = 0; i < count; ++i )
	{
		if ( !b3VoxelStateIsEmpty( v->states[i] ) )
		{
			filled += 1;
		}
	}
	return filled;
}

b3Vec3 b3VoxelsCentroid( const b3Voxels* v )
{
	b3Vec3 sum = b3Vec3_zero;
	int n = 0;
	for ( int z = 0; z < v->cz; ++z )
	{
		for ( int y = 0; y < v->cy; ++y )
		{
			for ( int x = 0; x < v->cx; ++x )
			{
				if ( !b3VoxelStateIsEmpty( v->states[b3VoxelsIndex( v, x, y, z )] ) )
				{
					sum = b3Add( sum, b3VoxelsVoxelCenter( v, x, y, z ) );
					n += 1;
				}
			}
		}
	}
	if ( n == 0 )
	{
		return b3AABB_Center( v->localAABB );
	}
	return b3MulSV( 1.0f / (float)n, sum );
}

b3MassData b3ComputeVoxelsMass( const b3Voxels* v, float density )
{
	b3MassData out = { 0 };

	// Each filled voxel is an identical uniform cuboid. Pass 1: mass + center of mass.
	b3Vec3 half = b3MulSV( 0.5f, v->voxelSize );
	float voxelVolume = v->voxelSize.x * v->voxelSize.y * v->voxelSize.z;
	float voxelMass = voxelVolume * density;

	b3Vec3 com = b3Vec3_zero;
	int n = 0;
	for ( int z = 0; z < v->cz; ++z )
	{
		for ( int y = 0; y < v->cy; ++y )
		{
			for ( int x = 0; x < v->cx; ++x )
			{
				if ( !b3VoxelStateIsEmpty( v->states[b3VoxelsIndex( v, x, y, z )] ) )
				{
					com = b3Add( com, b3VoxelsVoxelCenter( v, x, y, z ) );
					n += 1;
				}
			}
		}
	}

	if ( n == 0 )
	{
		return out;
	}

	com = b3MulSV( 1.0f / (float)n, com );

	// Inertia of one voxel about its own center.
	b3Matrix3 blockInertia = b3BoxInertia( voxelMass, b3Neg( half ), half );

	// Pass 2: sum shifted inertia tensors about the combined center of mass (parallel-axis).
	b3Matrix3 inertia = b3Mat3_zero;
	for ( int z = 0; z < v->cz; ++z )
	{
		for ( int y = 0; y < v->cy; ++y )
		{
			for ( int x = 0; x < v->cx; ++x )
			{
				if ( b3VoxelStateIsEmpty( v->states[b3VoxelsIndex( v, x, y, z )] ) )
				{
					continue;
				}
				b3Vec3 d = b3Sub( b3VoxelsVoxelCenter( v, x, y, z ), com );
				float dd = b3Dot( d, d );
				// shift = m * (|d|^2 * I - d (x) d)
				b3Matrix3 shift;
				shift.cx = (b3Vec3){ voxelMass * ( dd - d.x * d.x ), voxelMass * ( -d.x * d.y ), voxelMass * ( -d.x * d.z ) };
				shift.cy = (b3Vec3){ voxelMass * ( -d.y * d.x ), voxelMass * ( dd - d.y * d.y ), voxelMass * ( -d.y * d.z ) };
				shift.cz = (b3Vec3){ voxelMass * ( -d.z * d.x ), voxelMass * ( -d.z * d.y ), voxelMass * ( dd - d.z * d.z ) };
				inertia = b3AddMM( inertia, b3AddMM( blockInertia, shift ) );
			}
		}
	}

	out.mass = voxelMass * (float)n;
	out.center = com;
	out.inertia = inertia;
	return out;
}

// Sum of exposed face areas. Optionally weighted by |dot(faceNormal, planeNormal)| for projection.
static float b3VoxelsAreaImpl( const b3Voxels* v, const b3Vec3* planeNormal )
{
	// Face areas for the 6 faces (+x,-x,+y,-y,+z,-z).
	float areaX = v->voxelSize.y * v->voxelSize.z;
	float areaY = v->voxelSize.x * v->voxelSize.z;
	float areaZ = v->voxelSize.x * v->voxelSize.y;
	const b3Vec3 faceNormals[6] = {
		{ 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 },
	};
	const float faceAreas[6] = { areaX, areaX, areaY, areaY, areaZ, areaZ };

	float total = 0.0f;
	for ( int z = 0; z < v->cz; ++z )
	{
		for ( int y = 0; y < v->cy; ++y )
		{
			for ( int x = 0; x < v->cx; ++x )
			{
				b3VoxelState s = v->states[b3VoxelsIndex( v, x, y, z )];
				if ( b3VoxelStateIsEmpty( s ) )
				{
					continue;
				}
				uint8_t free = b3VoxelStateFreeFaces( s );
				for ( int f = 0; f < 6; ++f )
				{
					if ( free & ( 1 << f ) )
					{
						if ( planeNormal != NULL )
						{
							// Project: only the component facing the plane contributes.
							float w = b3Dot( faceNormals[f], *planeNormal );
							total += 0.5f * ( w > 0.0f ? w : -w ) * faceAreas[f];
						}
						else
						{
							total += faceAreas[f];
						}
					}
				}
			}
		}
	}
	return total;
}

float b3VoxelsSurfaceArea( const b3Voxels* v )
{
	return b3VoxelsAreaImpl( v, NULL );
}

float b3VoxelsProjectedArea( const b3Voxels* v, b3Vec3 planeNormal )
{
	return b3VoxelsAreaImpl( v, &planeNormal );
}

void b3DestroyVoxels( b3Voxels* v )
{
	if ( v == NULL )
	{
		return;
	}
	if ( v->states != NULL )
	{
		b3Free( v->states, (size_t)( v->cx * v->cy * v->cz ) * sizeof( b3VoxelState ) );
	}
	b3Free( v, sizeof( b3Voxels ) );
}
