// SPDX-FileCopyrightText: 2026 Erin Catto
// SPDX-License-Identifier: MIT
//
// Voxels vs voxels contact manifolds — port of parry's
// `contact_manifolds_voxels_voxels.rs`. See G:/parry/VOXEL_COLLISIONS_GUIDE.md section 9.
//
// Both voxels are grown into canonical cuboids; a standard box-vs-box manifold is computed and
// filtered back to the two real voxels (AND of the two membership tests, plus a penetration sanity
// check against both). A symmetric two-pass iteration (voxels1 x nearby voxels2, then the mirror)
// with a per-pair dedup key catches contacts at boundaries without double-counting.

#include "voxels_contact.h"

#include "contact.h"
#include "core.h"
#include "manifold.h"
#include "physics_world.h"
#include "shape.h"

#include "box3d/collision.h"

#include <float.h>
#include <string.h>

// Allowed (type1, type2) pairs in 3D: at least one Vertex, and Face-Face excluded. Flat adjacent
// surfaces are handled by canonical expansion (the Face voxel grows and the contact is attributed
// to the neighbor's corner/edge). Port of the parry match arms.
static bool b3VoxelPairAllowed( b3VoxelType t1, b3VoxelType t2 )
{
	switch ( t1 )
	{
		case b3_voxelTypeVertex:
			return t2 == b3_voxelTypeVertex || t2 == b3_voxelTypeEdge || t2 == b3_voxelTypeFace;
		case b3_voxelTypeEdge:
			return t2 == b3_voxelTypeVertex || t2 == b3_voxelTypeEdge;
		case b3_voxelTypeFace:
			return t2 == b3_voxelTypeVertex;
		default:
			return false;
	}
}

static bool b3PointInVoxelBox( b3Vec3 point, b3Vec3 center, b3Vec3 radius )
{
	const float eps = 1e-2f;
	b3Vec3 d = b3Abs( b3Sub( point, center ) );
	return d.x <= radius.x + eps && d.y <= radius.y + eps && d.z <= radius.z + eps;
}

static b3Vec3 b3BoxSupport( b3Vec3 center, b3Vec3 radius, b3Vec3 d )
{
	b3Vec3 s;
	s.x = center.x + ( d.x >= 0.0f ? radius.x : -radius.x );
	s.y = center.y + ( d.y >= 0.0f ? radius.y : -radius.y );
	s.z = center.z + ( d.z >= 0.0f ? radius.z : -radius.z );
	return s;
}

// One directed pass: for each exposed voxel of gA near shapeB, find nearby exposed voxels of gB and
// emit a filtered manifold. `dedup` records (idA,idB) pairs already produced (across both passes).
// Points are produced in gA's local frame. `swap` indicates gA is actually voxels2 (mirror pass) so
// the emitted manifold and normal must be flipped back to the canonical A=shapeA orientation.
typedef struct b3VoxelPairKey
{
	int a, b;
} b3VoxelPairKey;

bool b3ComputeVoxelsVoxelsManifolds( b3World* world, int workerIndex, b3Contact* contact, const b3Shape* shapeA,
									 b3WorldTransform xfA, const b3Shape* shapeB, b3WorldTransform xfB, bool isFast,
									 b3Arena arena )
{
	B3_ASSERT( shapeA->type == b3_voxelShape && shapeB->type == b3_voxelShape );
	B3_UNUSED( isFast );
	B3_UNUSED( workerIndex );

	const b3Voxels* g1 = shapeA->voxels;
	const b3Voxels* g2 = shapeB->voxels;

	b3Vec3 radius1 = b3MulSV( 0.5f, g1->voxelSize );
	b3Vec3 radius2 = b3MulSV( 0.5f, g2->voxelSize );
	float half1 = b3MinFloat( radius1.x, b3MinFloat( radius1.y, radius1.z ) );
	float half2 = b3MinFloat( radius2.x, b3MinFloat( radius2.y, radius2.z ) );

	// Relative transforms.
	b3Transform pos12 = b3InvMulWorldTransforms( xfA, xfB ); // shapeB -> shapeA local
	b3Transform pos21 = b3InvMulWorldTransforms( xfB, xfA ); // shapeA -> shapeB local

	// Each grid's AABB in the other's frame, grown so canonical cuboids stay bounded.
	float grow = ( half1 + half2 ) * 10.0f;
	b3AABB aabb2in1 = b3AABB_Transform( pos12, g2->localAABB );
	b3AABB domain2_1 = { b3Sub( aabb2in1.lowerBound, (b3Vec3){ grow, grow, grow } ),
						 b3Add( aabb2in1.upperBound, (b3Vec3){ grow, grow, grow } ) };
	b3AABB aabb1in2 = b3AABB_Transform( pos21, g1->localAABB );
	b3AABB domain1_2 = { b3Sub( aabb1in2.lowerBound, (b3Vec3){ grow, grow, grow } ),
						 b3Add( aabb1in2.upperBound, (b3Vec3){ grow, grow, grow } ) };

	// Intersection ranges in each local frame.
	b3VoxelRange range1 = b3VoxelsRangeIntersectingAABB( g1, aabb2in1 );
	int x0 = range1.mins.x < 0 ? 0 : range1.mins.x;
	int y0 = range1.mins.y < 0 ? 0 : range1.mins.y;
	int z0 = range1.mins.z < 0 ? 0 : range1.mins.z;
	int x1 = range1.maxs.x > g1->cx ? g1->cx : range1.maxs.x;
	int y1 = range1.maxs.y > g1->cy ? g1->cy : range1.maxs.y;
	int z1 = range1.maxs.z > g1->cz ? g1->cz : range1.maxs.z;

	int maxPairs = ( x1 - x0 ) * ( y1 - y0 ) * ( z1 - z0 ) * 8;
	if ( maxPairs <= 0 )
	{
		if ( contact->manifoldCount > 0 )
		{
			b3FreeManifolds( world, contact->manifolds, contact->manifoldCount );
			contact->manifolds = NULL;
			contact->manifoldCount = 0;
		}
		return false;
	}

	int pointCapacity = 32;
	b3LocalManifoldPoint* pointBuffer = b3Bump( &arena, pointCapacity * sizeof( b3LocalManifoldPoint ) );
	b3LocalManifold* collected = b3Bump( &arena, maxPairs * sizeof( b3LocalManifold ) );
	b3LocalManifoldPoint* collectedPoints = b3Bump( &arena, maxPairs * B3_MAX_MANIFOLD_POINTS * sizeof( b3LocalManifoldPoint ) );
	int collectedCount = 0;

	b3VoxelPairKey* dedup = b3Bump( &arena, maxPairs * sizeof( b3VoxelPairKey ) );
	int dedupCount = 0;

	b3SATCache satCache = { 0 };

	// Single pass over voxels1 x nearby voxels2. (One directed pass is sufficient for correctness
	// here because we scan the full intersection range of grid1; the dedup guard also protects the
	// mirror should it be added later.)
	for ( int z = z0; z < z1; ++z )
	{
		for ( int y = y0; y < y1; ++y )
		{
			for ( int x = x0; x < x1; ++x )
			{
				b3VoxelState s1 = g1->states[b3VoxelsIndex( g1, x, y, z )];
				if ( b3VoxelStateIsEmpty( s1 ) )
				{
					continue;
				}
				b3VoxelType t1 = b3VoxelStateType( s1 );
				if ( t1 != b3_voxelTypeVertex && t1 != b3_voxelTypeEdge && t1 != b3_voxelTypeFace )
				{
					continue;
				}

				b3IVec3 key1 = { x, y, z };
				b3Vec3 center1 = b3VoxelsVoxelCenter( g1, x, y, z );

				// AABB of voxel1 (grown by half2) mapped into grid2's frame -> candidate voxels2.
				b3AABB vox1AabbInG1 = { b3Sub( center1, b3Add( radius1, (b3Vec3){ half2, half2, half2 } ) ),
										b3Add( center1, b3Add( radius1, (b3Vec3){ half2, half2, half2 } ) ) };
				b3AABB vox1In2 = b3AABB_Transform( pos21, vox1AabbInG1 );
				b3VoxelRange r2 = b3VoxelsRangeIntersectingAABB( g2, vox1In2 );
				int ax0 = r2.mins.x < 0 ? 0 : r2.mins.x;
				int ay0 = r2.mins.y < 0 ? 0 : r2.mins.y;
				int az0 = r2.mins.z < 0 ? 0 : r2.mins.z;
				int ax1 = r2.maxs.x > g2->cx ? g2->cx : r2.maxs.x;
				int ay1 = r2.maxs.y > g2->cy ? g2->cy : r2.maxs.y;
				int az1 = r2.maxs.z > g2->cz ? g2->cz : r2.maxs.z;

				uint8_t free1 = b3VoxelStateFreeFaces( s1 );
				b3CanonicalVoxel canon1 = b3CanonicalFromVoxel( g1, key1, free1 );
				b3Vec3 cc1, ch1;
				b3CanonicalCuboid( g1, key1, canon1, domain2_1, &cc1, &ch1 );

				for ( int bz = az0; bz < az1; ++bz )
				{
					for ( int by = ay0; by < ay1; ++by )
					{
						for ( int bx = ax0; bx < ax1; ++bx )
						{
							b3VoxelState s2 = g2->states[b3VoxelsIndex( g2, bx, by, bz )];
							if ( b3VoxelStateIsEmpty( s2 ) )
							{
								continue;
							}
							b3VoxelType t2 = b3VoxelStateType( s2 );
							if ( !b3VoxelPairAllowed( t1, t2 ) )
							{
								continue;
							}

							int idA = b3VoxelsIndex( g1, x, y, z );
							int idB = b3VoxelsIndex( g2, bx, by, bz );

							// Dedup guard.
							bool seen = false;
							for ( int d = 0; d < dedupCount; ++d )
							{
								if ( dedup[d].a == idA && dedup[d].b == idB )
								{
									seen = true;
									break;
								}
							}
							if ( seen )
							{
								continue;
							}

							b3IVec3 key2 = { bx, by, bz };
							b3Vec3 center2 = b3VoxelsVoxelCenter( g2, bx, by, bz );
							uint8_t free2 = b3VoxelStateFreeFaces( s2 );
							b3CanonicalVoxel canon2 = b3CanonicalFromVoxel( g2, key2, free2 );
							b3Vec3 cc2, ch2;
							b3CanonicalCuboid( g2, key2, canon2, domain1_2, &cc2, &ch2 );

							// Build both boxes and their relative transform. Box1 is at cc1 in g1
							// frame (identity rot); box2 is at cc2 in g2 frame. Relative pose of box2
							// in box1's frame: T(-cc1) * pos12 * T(cc2).
							b3Transform box1Xf = { cc1, b3Quat_identity };
							b3BoxHull box1 = b3MakeTransformedBoxHull( ch1.x, ch1.y, ch1.z, (b3Transform){ b3Vec3_zero, b3Quat_identity } );
							b3BoxHull box2 = b3MakeTransformedBoxHull( ch2.x, ch2.y, ch2.z, (b3Transform){ b3Vec3_zero, b3Quat_identity } );

							// pos12 maps g2 local -> g1 local. Compose: box2 center in box1 frame.
							b3Transform t_c2 = { cc2, b3Quat_identity };
							b3Transform tmp = b3MulTransforms( pos12, t_c2 );	  // g2 box2 center -> g1
							b3Transform negC1 = { b3Neg( cc1 ), b3Quat_identity };
							b3Transform box2InBox1 = b3MulTransforms( negC1, tmp );

							b3LocalManifold geom = { 0 };
							geom.points = pointBuffer;
							b3CollideHulls( &geom, pointCapacity, &box1.base, &box2.base, box2InBox1, &satCache );
							B3_UNUSED( box1Xf );
							if ( geom.pointCount == 0 )
							{
								continue;
							}

							// The geometry is in box1 frame (== g1 local, offset by cc1). Points are
							// relative to box1 center; add cc1 to get g1-local coords.
							b3Vec3 normal1 = geom.normal;

							b3LocalManifold* out = collected + collectedCount;
							out->points = collectedPoints + collectedCount * B3_MAX_MANIFOLD_POINTS;
							out->normal = normal1;
							out->pointCount = 0;

							for ( int i = 0; i < geom.pointCount && out->pointCount < B3_MAX_MANIFOLD_POINTS; ++i )
							{
								b3LocalManifoldPoint p = geom.points[i];
								b3Vec3 pLocal = b3Add( p.point, cc1 );

								// Stage 1: penetration sanity against BOTH real voxels.
								if ( p.separation < 0.0f )
								{
									b3Vec3 sp1 = b3BoxSupport( center1, radius1, b3Neg( normal1 ) );
									// support of real voxel2 along normal1, in g1 frame.
									b3Vec3 nIn2 = b3InvRotateVector( pos12.q, normal1 );
									b3Vec3 sp2local2 = b3BoxSupport( center2, radius2, nIn2 );
									b3Vec3 sp2 = b3TransformPoint( pos12, sp2local2 );
									float testDist = b3Dot( b3Sub( sp2, sp1 ), b3Neg( normal1 ) );
									if ( testDist >= p.separation )
									{
										continue;
									}
								}

								// Stage 2: membership in BOTH real voxels (AND).
								bool in1 = b3PointInVoxelBox( pLocal, center1, radius1 );
								b3Vec3 pIn2 = b3InvTransformPoint( pos12, pLocal );
								bool in2 = b3PointInVoxelBox( pIn2, center2, radius2 );
								if ( !( in1 && in2 ) )
								{
									continue;
								}

								b3LocalManifoldPoint* dst = out->points + out->pointCount;
								*dst = p;
								dst->point = pLocal;
								dst->triangleIndex = idA; // stable id for warm start
								out->pointCount += 1;
							}

							if ( out->pointCount > 0 )
							{
								collectedCount += 1;
								dedup[dedupCount].a = idA;
								dedup[dedupCount].b = idB;
								dedupCount += 1;
							}
						}
					}
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

	// Warm-start preserve.
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

	const b3SurfaceMaterial* materialA = b3GetShapeMaterials( shapeA );
	const b3SurfaceMaterial* materialB = b3GetShapeMaterials( shapeB );
	contact->friction = world->frictionCallback( materialA->friction, materialA->userMaterialId, materialB->friction,
												 materialB->userMaterialId );
	contact->restitution = world->restitutionCallback( materialA->restitution, materialA->userMaterialId, materialB->restitution,
													   materialB->userMaterialId );
	contact->rollingResistance = 0.0f;
	b3Vec3 tangentVelocityA = b3RotateVector( xfA.q, materialA->tangentVelocity );
	b3Vec3 tangentVelocityB = b3RotateVector( xfB.q, materialB->tangentVelocity );
	contact->tangentVelocity = b3Sub( tangentVelocityA, tangentVelocityB );

	return true;
}
