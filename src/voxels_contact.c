// SPDX-FileCopyrightText: 2026 Erin Catto
// SPDX-License-Identifier: MIT
//
// Voxels vs convex-shape contact manifolds — port of parry's
// `contact_manifolds_voxels_shape.rs`. See G:/parry/VOXEL_COLLISIONS_GUIDE.md sections 6, 7, 10.
//
// The core trick: each exposed voxel is grown into a "canonical cuboid" that extends outward along
// its free (air-facing) axes. A standard convex-convex manifold is computed against that grown
// cuboid, then contacts are filtered back to the real single voxel with a two-stage test. This is
// what eliminates the internal-edges problem (objects snagging on seams between adjacent voxels).

#include "voxels_contact.h"

#include "contact.h"
#include "core.h"
#include "manifold.h"
#include "physics_world.h"
#include "shape.h"

#include "box3d/collision.h"

#include <float.h>
#include <string.h>

// Build the canonical voxel range for an exposed voxel. Port of CanonicalVoxelShape::from_voxel.
//
// The cuboid is pushed outward to the domain limit along each axis direction whose face is FILLED
// (not free) — i.e. it grows INTO the solid neighborhood, leaving the exposed faces bounded at the
// real voxel. For example a flat-floor voxel with only its +Y face exposed grows into a wide, deep
// slab whose top face coincides with the voxel top, so a body resting on it gets a clean +Y contact
// with no internal edges between adjacent top voxels. (Parry: `if !free_faces.contains(axis) { grow }`.)
b3CanonicalVoxel b3CanonicalFromVoxel( const b3Voxels* v, b3IVec3 key, uint8_t freeFaces )
{
	// Offset the domain by 1 so we can expand past the last voxel (matches parry).
	b3IVec3 mins = { -1, -1, -1 };
	b3IVec3 maxs = { v->cx, v->cy, v->cz };

	b3CanonicalVoxel c;
	c.low = key;
	c.high = key;

	if ( !( freeFaces & b3_voxelX_POS ) ) c.high.x = maxs.x;
	if ( !( freeFaces & b3_voxelX_NEG ) ) c.low.x = mins.x;
	if ( !( freeFaces & b3_voxelY_POS ) ) c.high.y = maxs.y;
	if ( !( freeFaces & b3_voxelY_NEG ) ) c.low.y = mins.y;
	if ( !( freeFaces & b3_voxelZ_POS ) ) c.high.z = maxs.z;
	if ( !( freeFaces & b3_voxelZ_NEG ) ) c.low.z = mins.z;

	return c;
}

// Materialize the canonical cuboid as a box hull. Expanded axes are clamped to `domain`
// (the other shape's AABB grown by a margin) so the "infinite" cuboid stays finite. Port of
// CanonicalVoxelShape::cuboid. Outputs the local-space center and half extents.
void b3CanonicalCuboid( const b3Voxels* v, b3IVec3 key, b3CanonicalVoxel canon, b3AABB domain, b3Vec3* outCenter,
						b3Vec3* outHalf )
{
	b3Vec3 radius = b3MulSV( 0.5f, v->voxelSize );
	b3Vec3 cmin = b3VoxelsCenterKey( v, canon.low );
	b3Vec3 cmax = b3VoxelsCenterKey( v, canon.high );

	// Clamp expanded sides to the other shape's loosened AABB.
	if ( canon.low.x != key.x ) cmin.x = b3MaxFloat( cmin.x, domain.lowerBound.x );
	if ( canon.low.y != key.y ) cmin.y = b3MaxFloat( cmin.y, domain.lowerBound.y );
	if ( canon.low.z != key.z ) cmin.z = b3MaxFloat( cmin.z, domain.lowerBound.z );
	if ( canon.high.x != key.x ) cmax.x = b3MinFloat( cmax.x, domain.upperBound.x );
	if ( canon.high.y != key.y ) cmax.y = b3MinFloat( cmax.y, domain.upperBound.y );
	if ( canon.high.z != key.z ) cmax.z = b3MinFloat( cmax.z, domain.upperBound.z );

	outHalf->x = 0.5f * ( cmax.x - cmin.x ) + radius.x;
	outHalf->y = 0.5f * ( cmax.y - cmin.y ) + radius.y;
	outHalf->z = 0.5f * ( cmax.z - cmin.z ) + radius.z;
	*outCenter = b3Lerp( cmin, cmax, 0.5f );
}

// Collide a canonical box hull (axis-aligned at the origin of the box frame) against shapeB.
// The caller has already moved shapeB into that box frame via transformBtoA. Fills geomManifold in
// the box frame. Returns point count.
static int b3CollideCanonicalWithShape( b3LocalManifold* geomManifold, int capacity, b3Vec3 half,
										const b3Shape* shapeB, b3Transform transformBtoA, b3ContactCache* cache )
{
	b3BoxHull box = b3MakeCollisionBoxHull( half.x, half.y, half.z );

	// transformBtoA maps shapeB into the box frame (identity pose, origin-centered).
	switch ( shapeB->type )
	{
		case b3_sphereShape:
			b3CollideHullAndSphere( geomManifold, capacity, &box.base, &shapeB->sphere, transformBtoA, &cache->simplexCache );
			break;
		case b3_capsuleShape:
			b3CollideHullAndCapsule( geomManifold, capacity, &box.base, &shapeB->capsule, transformBtoA, &cache->simplexCache );
			break;
		case b3_hullShape:
			b3CollideHulls( geomManifold, capacity, &box.base, shapeB->hull, transformBtoA, &cache->satCache );
			break;
		default:
			B3_ASSERT( false );
			return 0;
	}
	return geomManifold->pointCount;
}

// A single voxel cuboid at `center` with half extents = voxelSize/2, used for the membership test.
// Returns true if `point` (voxel-local) lies inside the inflated real voxel.
static bool b3PointInVoxel( b3Vec3 point, b3Vec3 center, b3Vec3 radius )
{
	// Parry inflates the test cuboid by an absolute 1e-2. Cap it at a fraction of the voxel radius
	// so the test does not degenerate to always-true for very small voxels (radius < ~0.04).
	float minR = b3MinFloat( radius.x, b3MinFloat( radius.y, radius.z ) );
	float eps = b3MinFloat( 1e-2f, 0.25f * minR );
	b3Vec3 d = b3Abs( b3Sub( point, center ) );
	return d.x <= radius.x + eps && d.y <= radius.y + eps && d.z <= radius.z + eps;
}

// Support point of an axis-aligned voxel cuboid (center, half extents) along direction d.
static b3Vec3 b3VoxelSupport( b3Vec3 center, b3Vec3 radius, b3Vec3 d )
{
	b3Vec3 s;
	s.x = center.x + ( d.x >= 0.0f ? radius.x : -radius.x );
	s.y = center.y + ( d.y >= 0.0f ? radius.y : -radius.y );
	s.z = center.z + ( d.z >= 0.0f ? radius.z : -radius.z );
	return s;
}

// Support point of shapeB, expressed in the voxel-local frame (frame A), along direction dirA.
// transformBtoA maps shapeB-local space into frame A. The proxy and R = mat(transformBtoA.q) are
// loop-invariant across all voxels/points, so the caller precomputes them once and passes them in.
static b3Vec3 b3SupportShapeBinA( const b3ShapeProxy* proxy, b3Matrix3 R, b3Transform transformBtoA, b3Vec3 dirA )
{
	float best = -FLT_MAX;
	b3Vec3 bestP = b3Vec3_zero;
	for ( int i = 0; i < proxy->count; ++i )
	{
		b3Vec3 pA = b3Add( b3MulMV( R, proxy->points[i] ), transformBtoA.p );
		float dot = b3Dot( pA, dirA );
		if ( dot > best )
		{
			best = dot;
			bestP = pA;
		}
	}
	// Account for the proxy's external radius (sphere/capsule).
	if ( proxy->radius > 0.0f )
	{
		bestP = b3MulAdd( bestP, proxy->radius, b3Normalize( dirA ) );
	}
	return bestP;
}

bool b3ComputeVoxelsManifolds( b3World* world, int workerIndex, b3Contact* contact, const b3Shape* shapeA,
							   b3WorldTransform xfA, const b3Shape* shapeB, b3WorldTransform xfB, bool isFast, b3Arena arena )
{
	B3_ASSERT( shapeA->type == b3_voxelShape );
	B3_UNUSED( isFast );
	B3_UNUSED( workerIndex );

	const b3Voxels* v = shapeA->voxels;
	float voxelHalfMin = 0.5f * b3MinFloat( v->voxelSize.x, b3MinFloat( v->voxelSize.y, v->voxelSize.z ) );

	// shapeB's AABB in the voxel-local frame, expanded by the speculative margin so near-touching
	// (not yet overlapping) contacts are still generated. Without this, fast bodies that stop just
	// short of a voxel face get no contact and tunnel on the next step.
	b3Transform transformBtoA = b3InvMulWorldTransforms( xfA, xfB );
	float margin = B3_MAX_AABB_MARGIN + B3_SPECULATIVE_DISTANCE;
	b3AABB aabbB_A = b3ComputeShapeAABB( shapeB, transformBtoA );
	aabbB_A.lowerBound = b3Sub( aabbB_A.lowerBound, (b3Vec3){ margin, margin, margin } );
	aabbB_A.upperBound = b3Add( aabbB_A.upperBound, (b3Vec3){ margin, margin, margin } );

	// Domain for clamping "infinite" canonical axes: shapeB's AABB grown by 10x voxel half-size.
	float grow = voxelHalfMin * 10.0f;
	b3AABB domain2_1 = { b3Sub( aabbB_A.lowerBound, (b3Vec3){ grow, grow, grow } ),
						 b3Add( aabbB_A.upperBound, (b3Vec3){ grow, grow, grow } ) };

	// Broad-phase: voxel range intersecting shapeB's AABB, clamped to the grid.
	b3VoxelRange range = b3VoxelsRangeIntersectingAABB( v, aabbB_A );
	int x0 = range.mins.x < 0 ? 0 : range.mins.x;
	int y0 = range.mins.y < 0 ? 0 : range.mins.y;
	int z0 = range.mins.z < 0 ? 0 : range.mins.z;
	int x1 = range.maxs.x > v->cx ? v->cx : range.maxs.x;
	int y1 = range.maxs.y > v->cy ? v->cy : range.maxs.y;
	int z1 = range.maxs.z > v->cz ? v->cz : range.maxs.z;

	b3Vec3 radius = b3MulSV( 0.5f, v->voxelSize );

	// Scratch: gather manifolds (up to one per voxel).
	int maxVoxels = ( x1 - x0 ) * ( y1 - y0 ) * ( z1 - z0 );
	if ( maxVoxels <= 0 )
	{
		if ( contact->manifoldCount > 0 )
		{
			b3FreeManifolds( world, contact->manifolds, contact->manifoldCount );
			contact->manifolds = NULL;
			contact->manifoldCount = 0;
		}
		return false;
	}

	// Per-voxel geometry manifold scratch.
	int pointCapacity = 32;
	b3LocalManifoldPoint* pointBuffer = b3Bump( &arena, pointCapacity * sizeof( b3LocalManifoldPoint ) );

	// Collected finished manifolds (frame A local, converted to world at the end).
	b3LocalManifold* collected = b3Bump( &arena, maxVoxels * sizeof( b3LocalManifold ) );
	b3LocalManifoldPoint* collectedPoints = b3Bump( &arena, maxVoxels * B3_MAX_MANIFOLD_POINTS * sizeof( b3LocalManifoldPoint ) );
	int collectedCount = 0;

	b3ContactCache localCache = { 0 };

	// Loop-invariant inputs to the stage-1 penetration filter (b3SupportShapeBinA), hoisted out of
	// the per-voxel / per-point inner loops.
	b3ShapeProxy proxyB = b3MakeShapeProxy( shapeB );
	b3Matrix3 matrixBtoA = b3MakeMatrixFromQuat( transformBtoA.q );

	// Memoize the box-vs-shapeB geometry manifold across voxels that share the same canonical cuboid.
	// A canonical cuboid's integer range depends only on a voxel's free faces and its position along
	// those free axes, so across a flat contact region every exposed voxel yields the identical cuboid
	// — and thus the identical geometry manifold. Only the per-voxel filtering below differs. Computing
	// the box-vs-shapeB SAT/clip once per distinct cuboid instead of once per voxel is the dominant
	// narrow-phase saving (a flat resting box touches hundreds of coplanar voxels sharing one cuboid).
	b3CanonicalVoxel prevCanon = { 0 };
	bool havePrevCanon = false;
	b3LocalManifold geom = { 0 };
	geom.points = pointBuffer;
	int count = 0;

	for ( int z = z0; z < z1; ++z )
	{
		for ( int y = y0; y < y1; ++y )
		{
			for ( int x = x0; x < x1; ++x )
			{
				b3VoxelState state = v->states[b3VoxelsIndex( v, x, y, z )];
				if ( b3VoxelStateIsEmpty( state ) )
				{
					continue;
				}
				b3VoxelType type = b3VoxelStateType( state );
				if ( type == b3_voxelTypeInterior )
				{
					continue; // interior voxels are fully surrounded (parry skips them)
				}

				b3IVec3 key = { x, y, z };
				uint8_t freeFaces = b3VoxelStateFreeFaces( state );
				b3CanonicalVoxel canon = b3CanonicalFromVoxel( v, key, freeFaces );

				b3Vec3 canonCenter, canonHalf;
				b3CanonicalCuboid( v, key, canon, domain2_1, &canonCenter, &canonHalf );

				// Recompute the geometry manifold only when the canonical cuboid changes (see above).
				// Equal integer canon ranges give bit-identical canonCenter/canonHalf/transformBtoBox,
				// so the cached geom is exactly what a recompute would produce.
				if ( !havePrevCanon || memcmp( &canon, &prevCanon, sizeof( canon ) ) != 0 )
				{
					// Move shapeB into the canonical box frame (box is at canonCenter, identity rotation).
					b3Transform boxToLocal = { canonCenter, b3Quat_identity };
					b3Transform transformBtoBox = b3InvMulTransforms( boxToLocal, transformBtoA );

					geom.pointCount = 0;
					count = b3CollideCanonicalWithShape( &geom, pointCapacity, canonHalf, shapeB, transformBtoBox, &localCache );
					prevCanon = canon;
					havePrevCanon = true;
				}
				if ( count == 0 )
				{
					continue;
				}

				// The geometry manifold is in the box's frame. Translate points/normal into the
				// voxel-local frame (frame A) by adding canonCenter; the normal is unchanged.
				b3Vec3 voxelCenter = b3VoxelsVoxelCenter( v, x, y, z );

				// Two-stage filter (guide section 10). We test points against the real single voxel.
				b3LocalManifold* out = collected + collectedCount;
				out->points = collectedPoints + collectedCount * B3_MAX_MANIFOLD_POINTS;
				out->normal = geom.normal;
				out->pointCount = 0;

				// The geometry normal is in the box frame; it equals the voxel-local normal (box has
				// identity rotation). Used by the stage-1 penetration sanity check.
				b3Vec3 normalA = geom.normal;

				for ( int i = 0; i < count && out->pointCount < B3_MAX_MANIFOLD_POINTS; ++i )
				{
					b3LocalManifoldPoint p = geom.points[i];
					// Point is in box frame -> voxel-local frame.
					b3Vec3 pLocal = b3Add( p.point, canonCenter );

					// Stage 1 (guide section 10): penetration sanity check against the REAL single
					// voxel, not the grown canonical cuboid. This discards phantom contacts the
					// canonical shape invents by extending into empty space (e.g. a thin floor voxel
					// whose column grows above the real surface).
					if ( p.separation < 0.0f )
					{
						b3Vec3 sp1 = b3VoxelSupport( voxelCenter, radius, b3Neg( normalA ) );
						b3Vec3 sp2 = b3SupportShapeBinA( &proxyB, matrixBtoA, transformBtoA, normalA );
						float testDist = b3Dot( b3Sub( sp2, sp1 ), b3Neg( normalA ) );
						if ( testDist >= p.separation )
						{
							continue; // spurious penetration invented by the grown cuboid
						}
					}

					// Stage 2: membership test. Keep only contacts within the inflated real voxel.
					if ( !b3PointInVoxel( pLocal, voxelCenter, radius ) )
					{
						continue;
					}

					b3LocalManifoldPoint* dst = out->points + out->pointCount;
					*dst = p;
					dst->point = pLocal;
					dst->triangleIndex = b3VoxelsIndex( v, x, y, z ); // stable per-voxel id for warm start
					out->pointCount += 1;
				}

				if ( out->pointCount > 0 )
				{
					collectedCount += 1;
				}
			}
		}
	}

	if ( collectedCount == 0 )
	{
		if ( contact->manifoldCount > 0 )
		{
			b3FreeManifolds( world, contact->manifolds, contact->manifoldCount );
			contact->manifolds = NULL;
			contact->manifoldCount = 0;
		}
		return false;
	}

	// Preserve previous manifolds for impulse warm-starting (match by featureId + voxel id).
	int oldManifoldCount = contact->manifoldCount;
	b3Manifold* oldManifolds = NULL;
	if ( oldManifoldCount > 0 )
	{
		oldManifolds = b3Bump( &arena, oldManifoldCount * sizeof( b3Manifold ) );
		memcpy( oldManifolds, contact->manifolds, oldManifoldCount * sizeof( b3Manifold ) );
	}

	if ( oldManifoldCount != collectedCount )
	{
		b3FreeManifolds( world, contact->manifolds, contact->manifoldCount );
		contact->manifolds = b3AllocateManifolds( world, collectedCount );
		contact->manifoldCount = (uint16_t)collectedCount;
	}
	else
	{
		memset( contact->manifolds, 0, contact->manifoldCount * sizeof( b3Manifold ) );
	}

	b3Matrix3 matrixA = b3MakeMatrixFromQuat( xfA.q );
	b3Vec3 offsetB = b3SubPos( xfA.p, xfB.p );

	for ( int i = 0; i < collectedCount; ++i )
	{
		b3LocalManifold* src = collected + i;
		b3Manifold* dst = contact->manifolds + i;
		dst->pointCount = src->pointCount;
		dst->normal = b3MulMV( matrixA, src->normal );

		for ( int j = 0; j < src->pointCount; ++j )
		{
			const b3LocalManifoldPoint* sp = src->points + j;
			b3ManifoldPoint* tp = dst->points + j;
			tp->anchorA = b3MulMV( matrixA, sp->point );
			tp->anchorB = b3Add( tp->anchorA, offsetB );
			tp->separation = sp->separation;
			tp->featureId = b3MakeFeatureId( sp->pair );
			tp->triangleIndex = sp->triangleIndex;
			tp->normalVelocity = 0.0f;
			tp->normalImpulse = 0.0f;
			tp->totalNormalImpulse = 0.0f;
			tp->persisted = false;

			// Warm start: match against old manifolds by (featureId, voxel id).
			for ( int m = 0; m < oldManifoldCount && !tp->persisted; ++m )
			{
				b3Manifold* om = oldManifolds + m;
				for ( int k = 0; k < om->pointCount; ++k )
				{
					if ( om->points[k].featureId == tp->featureId && om->points[k].triangleIndex == tp->triangleIndex )
					{
						tp->normalImpulse = om->points[k].normalImpulse;
						tp->persisted = true;
						om->points[k].triangleIndex = B3_NULL_INDEX;
						break;
					}
				}
			}
		}
	}

	// Friction / restitution: voxels use the single base material of shapeA.
	const b3SurfaceMaterial* materialA = b3GetShapeMaterials( shapeA );
	const b3SurfaceMaterial* materialB = b3GetShapeMaterials( shapeB );
	contact->friction = world->frictionCallback( materialA->friction, materialA->userMaterialId, materialB->friction,
												 materialB->userMaterialId );
	contact->restitution = world->restitutionCallback( materialA->restitution, materialA->userMaterialId, materialB->restitution,
													   materialB->userMaterialId );

	float radiusB = 0.0f;
	if ( shapeB->type == b3_sphereShape )
	{
		radiusB = shapeB->sphere.radius;
	}
	else if ( shapeB->type == b3_capsuleShape )
	{
		radiusB = shapeB->capsule.radius;
	}
	else if ( shapeB->type == b3_hullShape )
	{
		radiusB = shapeB->hull->innerRadius;
	}
	contact->rollingResistance = materialB->rollingResistance * radiusB;

	b3Vec3 tangentVelocityA = b3RotateVector( xfA.q, materialA->tangentVelocity );
	b3Vec3 tangentVelocityB = b3RotateVector( xfB.q, materialB->tangentVelocity );
	contact->tangentVelocity = b3Sub( tangentVelocityA, tangentVelocityB );

	return true;
}
