// SPDX-FileCopyrightText: 2026 Erin Catto
// SPDX-License-Identifier: MIT

#include "gfx/draw.h"
#include "imgui.h"
#include "sample.h"
#include "utils.h"

#include "box3d/box3d.h"

#include <cmath>
#include <cstdint>
#include <vector>

// The samples renderer draws every shape through the debug adapter, but the
// voxel grid is an opaque runtime type with no public read API, so the adapter
// skips it. These samples therefore keep the occupancy they hand to Box3D and
// draw the filled voxels themselves as unit cubes. Box3D still owns collision;
// this is purely presentation, so the two never disagree (the grids are static
// after creation, which is all the public API exposes anyway).
namespace
{

// sRGB 0xRRGGBB -> linear RGBA, matching the renderer's lit-shape color path so
// the voxel cubes shade like the debug-drawn primitives around them.
float SRGBToLinear( float c )
{
	return ( c <= 0.04045f ) ? ( c * ( 1.0f / 12.92f ) ) : powf( ( c + 0.055f ) * ( 1.0f / 1.055f ), 2.4f );
}

Vec4 LinearColor( uint32_t hex, float alpha = 1.0f )
{
	const float r = (float)( ( hex >> 16 ) & 0xFFu ) * ( 1.0f / 255.0f );
	const float g = (float)( ( hex >> 8 ) & 0xFFu ) * ( 1.0f / 255.0f );
	const float b = (float)( hex & 0xFFu ) * ( 1.0f / 255.0f );
	return MakeVec4( SRGBToLinear( r ), SRGBToLinear( g ), SRGBToLinear( b ), alpha );
}

Vec4 LerpColor( Vec4 a, Vec4 b, float t )
{
	return MakeVec4( a.x + ( b.x - a.x ) * t, a.y + ( b.y - a.y ) * t, a.z + ( b.z - a.z ) * t, a.w + ( b.w - a.w ) * t );
}

// A dense voxel occupancy grid plus the metadata needed to build a b3VoxelsDef
// and to render the filled cells. Index layout matches Box3D:
// index = x + y * cx + z * cx * cy.
struct VoxelModel
{
	int cx = 0;
	int cy = 0;
	int cz = 0;
	b3Vec3 voxelSize = { 1.0f, 1.0f, 1.0f };
	b3Vec3 origin = b3Vec3_zero;
	std::vector<uint8_t> occupancy;

	void Init( int nx, int ny, int nz, b3Vec3 size, b3Vec3 gridOrigin )
	{
		cx = nx;
		cy = ny;
		cz = nz;
		voxelSize = size;
		origin = gridOrigin;
		occupancy.assign( (size_t)cx * cy * cz, 0u );
	}

	int Index( int x, int y, int z ) const
	{
		return x + y * cx + z * cx * cy;
	}

	bool InBounds( int x, int y, int z ) const
	{
		return x >= 0 && x < cx && y >= 0 && y < cy && z >= 0 && z < cz;
	}

	void Set( int x, int y, int z, bool filled )
	{
		if ( InBounds( x, y, z ) )
		{
			occupancy[Index( x, y, z )] = filled ? 1u : 0u;
		}
	}

	bool Get( int x, int y, int z ) const
	{
		return InBounds( x, y, z ) && occupancy[Index( x, y, z )] != 0u;
	}

	// Local-space center of voxel (x, y, z).
	b3Vec3 Center( int x, int y, int z ) const
	{
		return { origin.x + ( (float)x + 0.5f ) * voxelSize.x, origin.y + ( (float)y + 0.5f ) * voxelSize.y,
				 origin.z + ( (float)z + 0.5f ) * voxelSize.z };
	}

	// A filled voxel is on the surface when at least one axis neighbor is empty
	// or off the grid. Interior cells are hidden behind their neighbors, so
	// skipping them keeps the cube instance count low without changing the look.
	bool IsSurface( int x, int y, int z ) const
	{
		if ( !Get( x, y, z ) )
		{
			return false;
		}
		return !Get( x + 1, y, z ) || !Get( x - 1, y, z ) || !Get( x, y + 1, z ) || !Get( x, y - 1, z ) ||
			   !Get( x, y, z + 1 ) || !Get( x, y, z - 1 );
	}

	b3ShapeId Attach( b3BodyId bodyId, const b3ShapeDef* shapeDef ) const
	{
		b3VoxelsDef def = {};
		def.cx = cx;
		def.cy = cy;
		def.cz = cz;
		def.voxelSize = voxelSize;
		def.origin = origin;
		def.occupancy = occupancy.data();
		return b3CreateVoxelShape( bodyId, shapeDef, &def );
	}

	// Draw the surface voxels of this grid at the given body transform. `low`
	// and `high` bracket a height gradient applied over the grid's Y extent.
	void Draw( b3WorldTransform bodyTransform, Vec4 low, Vec4 high ) const
	{
		const float invSpan = cy > 1 ? 1.0f / (float)( cy - 1 ) : 0.0f;
		for ( int z = 0; z < cz; ++z )
		{
			for ( int y = 0; y < cy; ++y )
			{
				for ( int x = 0; x < cx; ++x )
				{
					if ( !IsSurface( x, y, z ) )
					{
						continue;
					}

					b3Pos world = b3TransformWorldPoint( bodyTransform, Center( x, y, z ) );
					b3WorldTransform cube = { world, bodyTransform.q };
					Vec4 color = LerpColor( low, high, (float)y * invSpan );
					DrawCube( cube, voxelSize, color );
				}
			}
		}
	}
};

} // namespace

// A static voxel landscape built from a height function, with dynamic
// primitives (spheres, capsules, boxes) dropped onto it. Exercises the
// voxel-vs-primitive contact path against a large static grid.
class VoxelTerrain : public Sample
{
public:
	explicit VoxelTerrain( SampleContext* context )
		: Sample( context )
	{
		if ( context->restart == false )
		{
			m_camera->SetView( -45.0f, 30.0f, 45.0f, { 0.0f, 3.0f, 0.0f } );
		}

		const int nx = 48;
		const int ny = 20;
		const int nz = 48;
		const float s = 0.5f;
		const b3Vec3 size = { s, s, s };

		// Center the footprint on the body origin and rest the base on y = 0.
		const b3Vec3 gridOrigin = { -0.5f * nx * s, 0.0f, -0.5f * nz * s };
		m_terrain.Init( nx, ny, nz, size, gridOrigin );

		// Rolling hills: two sine ridges plus a broad bump in the middle.
		for ( int x = 0; x < nx; ++x )
		{
			for ( int z = 0; z < nz; ++z )
			{
				float fx = (float)x / (float)( nx - 1 );
				float fz = (float)z / (float)( nz - 1 );
				float h = 5.0f + 3.0f * sinf( fx * 6.0f ) * cosf( fz * 5.0f );

				// A central mound so dropped bodies have something to roll off.
				float dx = fx - 0.5f;
				float dz = fz - 0.5f;
				h += 5.0f * expf( -( dx * dx + dz * dz ) * 12.0f );

				int top = (int)h;
				if ( top < 1 )
				{
					top = 1;
				}
				if ( top > ny - 1 )
				{
					top = ny - 1;
				}
				for ( int y = 0; y <= top; ++y )
				{
					m_terrain.Set( x, y, z, true );
				}
			}
		}

		b3BodyDef bodyDef = b3DefaultBodyDef();
		m_terrainBody = b3CreateBody( m_worldId, &bodyDef );

		b3ShapeDef shapeDef = b3DefaultShapeDef();
		shapeDef.baseMaterial.friction = 0.6f;
		m_terrain.Attach( m_terrainBody, &shapeDef );

		Launch();
	}

	void Launch()
	{
		b3BodyDef bodyDef = b3DefaultBodyDef();
		bodyDef.type = b3_dynamicBody;
		b3ShapeDef shapeDef = b3DefaultShapeDef();
		shapeDef.baseMaterial.restitution = 0.1f;

		g_randomSeed = 42;
		for ( int i = 0; i < m_dropCount; ++i )
		{
			bodyDef.position = { RandomFloatRange( -8.0f, 8.0f ), 16.0f + 0.6f * i, RandomFloatRange( -8.0f, 8.0f ) };
			bodyDef.rotation = b3MakeQuatFromAxisAngle( { 0.0f, 1.0f, 0.0f }, RandomFloatRange( 0.0f, B3_PI ) );
			b3BodyId body = b3CreateBody( m_worldId, &bodyDef );

			int kind = i % 3;
			if ( kind == 0 )
			{
				b3Sphere sphere = { b3Vec3_zero, 0.6f };
				b3CreateSphereShape( body, &shapeDef, &sphere );
			}
			else if ( kind == 1 )
			{
				b3Capsule capsule = { { -0.5f, 0.0f, 0.0f }, { 0.5f, 0.0f, 0.0f }, 0.4f };
				b3CreateCapsuleShape( body, &shapeDef, &capsule );
			}
			else
			{
				b3BoxHull box = b3MakeBoxHull( 0.6f, 0.6f, 0.6f );
				b3CreateHullShape( body, &shapeDef, &box.base );
			}
		}
	}

	void Render() override
	{
		m_terrain.Draw( b3Body_GetTransform( m_terrainBody ), LinearColor( 0x3f6d3f ), LinearColor( 0xd8e6ee ) );
		Sample::Render();
	}

	bool DrawControls() override
	{
		ImGui::SliderInt( "Drop count", &m_dropCount, 1, 40 );
		if ( ImGui::Button( "Drop bodies" ) )
		{
			Launch();
		}
		return true;
	}

	static Sample* Create( SampleContext* context )
	{
		return new VoxelTerrain( context );
	}

	VoxelModel m_terrain;
	b3BodyId m_terrainBody;
	int m_dropCount = 20;
};

static int sampleVoxelTerrain = RegisterSample( "Voxels", "Terrain", VoxelTerrain::Create );

// Dynamic voxel chunks falling onto a static voxel floor and stacking on each
// other. Exercises voxel dynamic bodies and the voxels-vs-voxels contact path.
class VoxelChunks : public Sample
{
public:
	explicit VoxelChunks( SampleContext* context )
		: Sample( context )
	{
		if ( context->restart == false )
		{
			m_camera->SetView( -50.0f, 20.0f, 30.0f, { 0.0f, 4.0f, 0.0f } );
		}

		// Static voxel floor slab.
		const float s = 0.5f;
		const b3Vec3 size = { s, s, s };
		const int fx = 40;
		const int fz = 40;
		m_floor.Init( fx, 2, fz, size, { -0.5f * fx * s, -1.0f, -0.5f * fz * s } );
		for ( int x = 0; x < fx; ++x )
		{
			for ( int z = 0; z < fz; ++z )
			{
				m_floor.Set( x, 0, z, true );
				m_floor.Set( x, 1, z, true );
			}
		}

		b3BodyDef bodyDef = b3DefaultBodyDef();
		m_floorBody = b3CreateBody( m_worldId, &bodyDef );
		b3ShapeDef shapeDef = b3DefaultShapeDef();
		shapeDef.baseMaterial.friction = 0.7f;
		m_floor.Attach( m_floorBody, &shapeDef );

		// One reusable "plus"-shaped chunk template, centered on its body origin
		// so it tumbles about its middle.
		const int cn = 3;
		m_chunk.Init( cn, cn, cn, size, { -0.5f * cn * s, -0.5f * cn * s, -0.5f * cn * s } );
		for ( int x = 0; x < cn; ++x )
		{
			for ( int y = 0; y < cn; ++y )
			{
				for ( int z = 0; z < cn; ++z )
				{
					// Fill the three axis-aligned bars through the center.
					int mid = cn / 2;
					bool onX = ( y == mid && z == mid );
					bool onY = ( x == mid && z == mid );
					bool onZ = ( x == mid && y == mid );
					if ( onX || onY || onZ )
					{
						m_chunk.Set( x, y, z, true );
					}
				}
			}
		}

		Spawn();
	}

	void Spawn()
	{
		b3BodyDef bodyDef = b3DefaultBodyDef();
		bodyDef.type = b3_dynamicBody;
		b3ShapeDef shapeDef = b3DefaultShapeDef();
		shapeDef.baseMaterial.friction = 0.7f;

		g_randomSeed = 7;
		for ( int i = 0; i < m_chunkCount; ++i )
		{
			bodyDef.position = { RandomFloatRange( -4.0f, 4.0f ), 5.0f + 2.5f * i, RandomFloatRange( -4.0f, 4.0f ) };
			bodyDef.rotation = b3MakeQuatFromAxisAngle( b3Normalize( { RandomFloatRange( -1.0f, 1.0f ), 1.0f,
																	   RandomFloatRange( -1.0f, 1.0f ) } ),
														RandomFloatRange( 0.0f, B3_PI ) );
			b3BodyId body = b3CreateBody( m_worldId, &bodyDef );
			b3ShapeId shapeId = m_chunk.Attach( body, &shapeDef );
			(void)shapeId;
			m_chunkBodies.push_back( body );
		}
	}

	void Render() override
	{
		m_floor.Draw( b3Body_GetTransform( m_floorBody ), LinearColor( 0x4a4a55 ), LinearColor( 0x6d6d7a ) );

		Vec4 low = LinearColor( 0xd98032 );
		Vec4 high = LinearColor( 0xf2d06b );
		for ( b3BodyId body : m_chunkBodies )
		{
			m_chunk.Draw( b3Body_GetTransform( body ), low, high );
		}
		Sample::Render();
	}

	bool DrawControls() override
	{
		ImGui::SliderInt( "Chunk count", &m_chunkCount, 1, 12 );
		if ( ImGui::Button( "Spawn chunks" ) )
		{
			Spawn();
		}
		return true;
	}

	static Sample* Create( SampleContext* context )
	{
		return new VoxelChunks( context );
	}

	VoxelModel m_floor;
	VoxelModel m_chunk;
	b3BodyId m_floorBody;
	std::vector<b3BodyId> m_chunkBodies;
	int m_chunkCount = 6;
};

static int sampleVoxelChunks = RegisterSample( "Voxels", "Dynamic Chunks", VoxelChunks::Create );

// A ray sweeping across a static voxel sphere, cast every step through
// b3World_CastRayClosest. Exercises the voxel ray-cast query and visualizes the
// hit point and surface normal.
class VoxelRayCast : public Sample
{
public:
	explicit VoxelRayCast( SampleContext* context )
		: Sample( context )
	{
		if ( context->restart == false )
		{
			m_camera->SetView( -90.0f, 6.0f, 24.0f, { 0.0f, 4.0f, 0.0f } );
		}

		const int n = 24;
		const float s = 0.5f;
		const b3Vec3 size = { s, s, s };
		m_blob.Init( n, n, n, size, { -0.5f * n * s, 4.0f - 0.5f * n * s, -0.5f * n * s } );

		// Carve a solid sphere, then bore a couple of tunnels so the ray has
		// concavities to catch on.
		const float r = 0.5f * n;
		const float rr = r * r;
		for ( int x = 0; x < n; ++x )
		{
			for ( int y = 0; y < n; ++y )
			{
				for ( int z = 0; z < n; ++z )
				{
					float dx = (float)x - r + 0.5f;
					float dy = (float)y - r + 0.5f;
					float dz = (float)z - r + 0.5f;
					float d2 = dx * dx + dy * dy + dz * dz;
					bool inside = d2 <= rr;
					// Two cylindrical bores through the center.
					bool boreZ = ( dx * dx + dy * dy ) < 6.0f;
					bool boreX = ( dy * dy + dz * dz ) < 4.0f;
					m_blob.Set( x, y, z, inside && !boreZ && !boreX );
				}
			}
		}

		b3BodyDef bodyDef = b3DefaultBodyDef();
		m_blobBody = b3CreateBody( m_worldId, &bodyDef );
		b3ShapeDef shapeDef = b3DefaultShapeDef();
		m_blob.Attach( m_blobBody, &shapeDef );
	}

	void Render() override
	{
		m_blob.Draw( b3Body_GetTransform( m_blobBody ), LinearColor( 0x6a5acd ), LinearColor( 0xb0a8ee ) );

		// Sweep the ray origin in a circle around the blob, aiming at its center.
		m_angle += 0.01f;
		float radius = 12.0f;
		b3Pos center = { 0.0f, 4.0f, 0.0f };
		b3Pos origin = { center.x + radius * cosf( m_angle ), 4.0f + 3.0f * sinf( m_angle * 0.7f ), center.z + radius * sinf( m_angle ) };
		b3Vec3 translation = { center.x - origin.x, center.y - origin.y, center.z - origin.z };

		b3RayResult result = b3World_CastRayClosest( m_worldId, origin, translation, b3DefaultQueryFilter() );

		if ( result.hit )
		{
			b3Pos hit = result.point;
			DrawLine( origin, hit, LinearColor( 0xffd400 ) );
			DrawPoint( hit, 8.0f, LinearColor( 0xff3030 ) );
			b3Pos normalEnd = { hit.x + result.normal.x, hit.y + result.normal.y, hit.z + result.normal.z };
			DrawArrow( hit, normalEnd, LinearColor( 0x30ff30 ) );
		}
		else
		{
			b3Pos end = { origin.x + translation.x, origin.y + translation.y, origin.z + translation.z };
			DrawLine( origin, end, LinearColor( 0x808080 ) );
		}

		DrawTextLine( result.hit ? "ray: hit" : "ray: miss" );
		Sample::Render();
	}

	static Sample* Create( SampleContext* context )
	{
		return new VoxelRayCast( context );
	}

	VoxelModel m_blob;
	b3BodyId m_blobBody;
	float m_angle = 0.0f;
};

static int sampleVoxelRayCast = RegisterSample( "Voxels", "Ray Cast", VoxelRayCast::Create );
