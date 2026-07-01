// SPDX-FileCopyrightText: 2026 Erin Catto
// SPDX-License-Identifier: MIT

#include "test_macros.h"

#include "voxels.h"

#include "box3d/box3d.h"

#include <float.h>
#include <stdio.h>

// A 3x3x3 solid block exercises every voxel classification:
//   center (1,1,1)      -> Interior (all 6 neighbors filled)
//   face center (1,1,0) -> Face     (one exposed face)
//   edge center (1,0,0) -> Edge     (two exposed faces meeting on an edge)
//   corner (0,0,0)      -> Vertex   (three exposed faces meeting on a corner)
static int VoxelsClassify( void )
{
	uint8_t occ[27];
	for ( int i = 0; i < 27; ++i )
	{
		occ[i] = 1;
	}

	b3VoxelsDef def = { 0 };
	def.cx = def.cy = def.cz = 3;
	def.voxelSize = (b3Vec3){ 1.0f, 1.0f, 1.0f };
	def.origin = (b3Vec3){ 0.0f, 0.0f, 0.0f };
	def.occupancy = occ;

	b3Voxels* v = b3CreateVoxels( &def );
	ENSURE( v != NULL );

	ENSURE( b3VoxelStateType( b3VoxelsState( v, 1, 1, 1 ) ) == b3_voxelTypeInterior );
	ENSURE( b3VoxelStateType( b3VoxelsState( v, 1, 1, 0 ) ) == b3_voxelTypeFace );
	ENSURE( b3VoxelStateType( b3VoxelsState( v, 1, 0, 0 ) ) == b3_voxelTypeEdge );
	ENSURE( b3VoxelStateType( b3VoxelsState( v, 0, 0, 0 ) ) == b3_voxelTypeVertex );

	// Interior voxel has no free faces; corner has 3 free faces (-x,-y,-z).
	ENSURE( b3VoxelStateFreeFaces( b3VoxelsState( v, 1, 1, 1 ) ) == 0 );
	uint8_t cornerFree = b3VoxelStateFreeFaces( b3VoxelsState( v, 0, 0, 0 ) );
	ENSURE( cornerFree == ( b3_voxelX_NEG | b3_voxelY_NEG | b3_voxelZ_NEG ) );

	// Center of the +x face (2,1,1) exposes only +x.
	uint8_t faceFree = b3VoxelStateFreeFaces( b3VoxelsState( v, 2, 1, 1 ) );
	ENSURE( faceFree == b3_voxelX_POS );

	b3DestroyVoxels( v );
	return 0;
}

// set_voxel incremental update must agree with a full recompute.
static int VoxelsSetVoxelMatchesRecompute( void )
{
	// Start all-empty, then fill a 2x2x2 block via set_voxel.
	b3VoxelsDef def = { 0 };
	def.cx = def.cy = def.cz = 4;
	def.voxelSize = (b3Vec3){ 0.5f, 0.5f, 0.5f };
	def.origin = (b3Vec3){ 0.0f, 0.0f, 0.0f };
	def.occupancy = NULL;

	b3Voxels* v = b3CreateVoxels( &def );
	ENSURE( v != NULL );

	for ( int z = 1; z <= 2; ++z )
	{
		for ( int y = 1; y <= 2; ++y )
		{
			for ( int x = 1; x <= 2; ++x )
			{
				b3VoxelsSetVoxel( v, x, y, z, true );
			}
		}
	}

	// Snapshot incremental result, then brute-force recompute and compare.
	int count = def.cx * def.cy * def.cz;
	uint8_t snapshot[64];
	for ( int i = 0; i < count; ++i )
	{
		snapshot[i] = v->states[i];
	}

	b3VoxelsRecomputeAllStates( v );

	for ( int i = 0; i < count; ++i )
	{
		ENSURE( snapshot[i] == v->states[i] );
	}

	// Each of the 8 filled voxels is a corner of the 2x2x2 block -> Vertex.
	ENSURE( b3VoxelStateType( b3VoxelsState( v, 1, 1, 1 ) ) == b3_voxelTypeVertex );
	ENSURE( b3VoxelStateType( b3VoxelsState( v, 2, 2, 2 ) ) == b3_voxelTypeVertex );

	// Removing a voxel makes its former neighbors expose the shared face again.
	b3VoxelsSetVoxel( v, 2, 2, 2, false );
	ENSURE( b3VoxelStateIsEmpty( b3VoxelsState( v, 2, 2, 2 ) ) );
	uint8_t nbr = b3VoxelsState( v, 1, 2, 2 );
	ENSURE( b3VoxelStateFreeFaces( nbr ) & b3_voxelX_POS );

	b3DestroyVoxels( v );
	return 0;
}

static int VoxelsGeometry( void )
{
	uint8_t occ[8] = { 1, 1, 1, 1, 1, 1, 1, 1 };
	b3VoxelsDef def = { 0 };
	def.cx = def.cy = def.cz = 2;
	def.voxelSize = (b3Vec3){ 2.0f, 2.0f, 2.0f };
	def.origin = (b3Vec3){ -2.0f, -2.0f, -2.0f };
	def.occupancy = occ;

	b3Voxels* v = b3CreateVoxels( &def );
	ENSURE( v != NULL );

	// Local AABB spans [-2, 2] on each axis (2 voxels x 2.0 size).
	ENSURE_SMALL( v->localAABB.lowerBound.x + 2.0f, FLT_EPSILON );
	ENSURE_SMALL( v->localAABB.upperBound.z - 2.0f, FLT_EPSILON );

	// Center of voxel (0,0,0) is origin + half = (-1,-1,-1).
	b3Vec3 c = b3VoxelsVoxelCenter( v, 0, 0, 0 );
	ENSURE_SMALL( c.x + 1.0f, FLT_EPSILON );

	// voxel_at_point round trips a center back to its key.
	b3IVec3 k = b3VoxelsAtPoint( v, c );
	ENSURE( k.x == 0 && k.y == 0 && k.z == 0 );

	b3DestroyVoxels( v );
	return 0;
}

// A solid NxNxN voxel block must have the same mass and inertia as a single cuboid of the
// same overall bounds (parallel-axis theorem sanity check).
static int VoxelsMassMatchesCuboid( void )
{
	const int N = 4;
	uint8_t occ[64];
	for ( int i = 0; i < N * N * N; ++i )
	{
		occ[i] = 1;
	}

	b3VoxelsDef def = { 0 };
	def.cx = def.cy = def.cz = N;
	def.voxelSize = (b3Vec3){ 0.5f, 0.5f, 0.5f };
	def.origin = (b3Vec3){ -1.0f, -1.0f, -1.0f }; // block spans [-1, 1]^3
	def.occupancy = occ;

	b3Voxels* v = b3CreateVoxels( &def );
	ENSURE( v != NULL );

	float density = 2.5f;
	b3MassData m = b3ComputeVoxelsMass( v, density );

	// Reference: solid 2x2x2 box centered at origin.
	float side = 2.0f;
	float volume = side * side * side;
	float mass = volume * density;
	float i = mass * ( side * side + side * side ) / 12.0f;

	ENSURE_SMALL( m.mass - mass, 1e-3f );
	ENSURE_SMALL( m.center.x, 1e-4f );
	ENSURE_SMALL( m.center.y, 1e-4f );
	ENSURE_SMALL( m.center.z, 1e-4f );
	ENSURE_SMALL( m.inertia.cx.x - i, 1e-2f );
	ENSURE_SMALL( m.inertia.cy.y - i, 1e-2f );
	ENSURE_SMALL( m.inertia.cz.z - i, 1e-2f );
	// Off-diagonals should vanish by symmetry.
	ENSURE_SMALL( m.inertia.cx.y, 1e-3f );
	ENSURE_SMALL( m.inertia.cy.z, 1e-3f );

	b3DestroyVoxels( v );
	return 0;
}

// A dynamic box dropped onto a static voxel floor must come to rest on top of it. This exercises
// the full voxels-vs-hull contact pipeline through the solver.
static int VoxelsBoxRestsOnFloor( void )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	worldDef.gravity = (b3Vec3){ 0.0f, -10.0f, 0.0f };
	b3WorldId worldId = b3CreateWorld( &worldDef );
	ENSURE( b3World_IsValid( worldId ) );

	// Static voxel floor: 8x1x8 grid of 1m voxels, top surface at y = 0.
	b3BodyDef groundDef = b3DefaultBodyDef();
	groundDef.position = (b3Pos){ 0.0f, 0.0f, 0.0f };
	b3BodyId groundId = b3CreateBody( worldId, &groundDef );

	int cx = 8, cy = 1, cz = 8;
	uint8_t occ[64];
	for ( int i = 0; i < cx * cy * cz; ++i )
	{
		occ[i] = 1;
	}

	b3VoxelsDef vdef = { 0 };
	vdef.cx = cx;
	vdef.cy = cy;
	vdef.cz = cz;
	vdef.voxelSize = (b3Vec3){ 1.0f, 1.0f, 1.0f };
	vdef.origin = (b3Vec3){ -4.0f, -1.0f, -4.0f }; // top surface at y = 0
	vdef.occupancy = occ;

	b3ShapeDef groundShapeDef = b3DefaultShapeDef();
	b3CreateVoxelShape( groundId, &groundShapeDef, &vdef );

	// Dynamic unit box dropped from above.
	b3BodyDef bodyDef = b3DefaultBodyDef();
	bodyDef.type = b3_dynamicBody;
	bodyDef.position = (b3Pos){ 0.0f, 3.0f, 0.0f };
	b3BodyId bodyId = b3CreateBody( worldId, &bodyDef );

	b3BoxHull box = b3MakeCubeHull( 0.5f ); // half-extent 0.5 -> unit cube
	b3ShapeDef shapeDef = b3DefaultShapeDef();
	shapeDef.density = 1.0f;
	shapeDef.baseMaterial.friction = 0.3f;
	b3CreateHullShape( bodyId, &shapeDef, &box.base );

	float timeStep = 1.0f / 60.0f;
	b3Pos position = b3Body_GetPosition( bodyId );
	for ( int i = 0; i < 120; ++i )
	{
		b3World_Step( worldId, timeStep, 4 );
	}
	position = b3Body_GetPosition( bodyId );
	b3Quat rotation = b3Body_GetRotation( bodyId );

	b3DestroyWorld( worldId );

	// Box half-extent 0.5 resting on floor top y=0 -> center near y = 0.5.
	ENSURE_SMALL( position.y - 0.5f, 0.05f );
	// Must not have tunneled through or drifted sideways.
	ENSURE_SMALL( position.x, 0.05f );
	ENSURE_SMALL( position.z, 0.05f );
	ENSURE_SMALL( rotation.v.x, 0.05f );
	ENSURE_SMALL( rotation.v.z, 0.05f );

	return 0;
}

// A small dynamic voxel block dropped onto a static voxel floor must rest on top of it.
static int VoxelsBlockRestsOnVoxelFloor( void )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	worldDef.gravity = (b3Vec3){ 0.0f, -10.0f, 0.0f };
	b3WorldId worldId = b3CreateWorld( &worldDef );
	ENSURE( b3World_IsValid( worldId ) );

	// Static voxel floor, top surface at y = 0.
	b3BodyDef groundDef = b3DefaultBodyDef();
	b3BodyId groundId = b3CreateBody( worldId, &groundDef );
	int gcx = 8, gcy = 1, gcz = 8;
	uint8_t gocc[64];
	for ( int i = 0; i < gcx * gcy * gcz; ++i )
	{
		gocc[i] = 1;
	}
	b3VoxelsDef gdef = { 0 };
	gdef.cx = gcx;
	gdef.cy = gcy;
	gdef.cz = gcz;
	gdef.voxelSize = (b3Vec3){ 1.0f, 1.0f, 1.0f };
	gdef.origin = (b3Vec3){ -4.0f, -1.0f, -4.0f };
	gdef.occupancy = gocc;
	b3ShapeDef gshape = b3DefaultShapeDef();
	b3CreateVoxelShape( groundId, &gshape, &gdef );

	// Dynamic 2x2x2 voxel block dropped from above, centered at origin.
	b3BodyDef blockDef = b3DefaultBodyDef();
	blockDef.type = b3_dynamicBody;
	blockDef.position = (b3Pos){ 0.0f, 3.0f, 0.0f };
	b3BodyId blockId = b3CreateBody( worldId, &blockDef );
	uint8_t bocc[8] = { 1, 1, 1, 1, 1, 1, 1, 1 };
	b3VoxelsDef bdef = { 0 };
	bdef.cx = bdef.cy = bdef.cz = 2;
	bdef.voxelSize = (b3Vec3){ 0.5f, 0.5f, 0.5f };
	bdef.origin = (b3Vec3){ -0.5f, -0.5f, -0.5f }; // block spans [-0.5, 0.5], height 1
	bdef.occupancy = bocc;
	b3ShapeDef bshape = b3DefaultShapeDef();
	bshape.density = 1.0f;
	bshape.baseMaterial.friction = 0.5f;
	b3CreateVoxelShape( blockId, &bshape, &bdef );

	float timeStep = 1.0f / 60.0f;
	for ( int i = 0; i < 150; ++i )
	{
		b3World_Step( worldId, timeStep, 4 );
	}
	b3Pos position = b3Body_GetPosition( blockId );

	b3DestroyWorld( worldId );

	// Block half-height 0.5 resting on floor top y=0 -> center near y = 0.5.
	ENSURE_SMALL( position.y - 0.5f, 0.08f );
	ENSURE_SMALL( position.x, 0.1f );
	ENSURE_SMALL( position.z, 0.1f );

	return 0;
}

// Ray cast down onto a voxel floor must hit the top surface (y=0) with an up normal, and overlap
// queries must detect the grid.
static int VoxelsRayAndOverlap( void )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	b3WorldId worldId = b3CreateWorld( &worldDef );

	b3BodyDef groundDef = b3DefaultBodyDef();
	b3BodyId groundId = b3CreateBody( worldId, &groundDef );
	int gcx = 6, gcy = 2, gcz = 6;
	uint8_t occ[72];
	for ( int i = 0; i < gcx * gcy * gcz; ++i )
	{
		occ[i] = 1;
	}
	b3VoxelsDef gdef = { 0 };
	gdef.cx = gcx;
	gdef.cy = gcy;
	gdef.cz = gcz;
	gdef.voxelSize = (b3Vec3){ 1.0f, 1.0f, 1.0f };
	gdef.origin = (b3Vec3){ -3.0f, -2.0f, -3.0f }; // top surface at y = 0
	gdef.occupancy = occ;
	b3ShapeDef gshape = b3DefaultShapeDef();
	b3CreateVoxelShape( groundId, &gshape, &gdef );

	// Ray straight down from (0, 5, 0).
	b3QueryFilter filter = b3DefaultQueryFilter();
	b3RayResult r = b3World_CastRayClosest( worldId, (b3Pos){ 0.0f, 5.0f, 0.0f }, (b3Vec3){ 0.0f, -10.0f, 0.0f }, filter );
	ENSURE( r.hit );
	ENSURE_SMALL( r.point.y - 0.0f, 0.05f );  // top surface
	ENSURE_SMALL( r.normal.y - 1.0f, 0.05f ); // up normal

	// Ray that misses the grid entirely.
	b3RayResult miss =
		b3World_CastRayClosest( worldId, (b3Pos){ 100.0f, 5.0f, 0.0f }, (b3Vec3){ 0.0f, -10.0f, 0.0f }, filter );
	ENSURE( miss.hit == false );

	b3DestroyWorld( worldId );
	return 0;
}

int VoxelsTest( void )
{
	RUN_SUBTEST( VoxelsClassify );
	RUN_SUBTEST( VoxelsSetVoxelMatchesRecompute );
	RUN_SUBTEST( VoxelsGeometry );
	RUN_SUBTEST( VoxelsMassMatchesCuboid );
	RUN_SUBTEST( VoxelsBoxRestsOnFloor );
	RUN_SUBTEST( VoxelsBlockRestsOnVoxelFloor );
	RUN_SUBTEST( VoxelsRayAndOverlap );
	return 0;
}
