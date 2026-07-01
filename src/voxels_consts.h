// SPDX-FileCopyrightText: 2025 Erin Catto
// SPDX-License-Identifier: MIT
//
// Voxel collision tables — ported verbatim from parry's `src/shape/voxels/voxels_consts.rs`.
// These encode the internal-edges fix: each voxel's 6-neighborhood byte indexes into three
// lookup tables that classify the voxel and describe its exposed convex features per octant.
// DO NOT edit by hand — the bit layout must match parry exactly or projections invert.
// See G:/parry/VOXEL_COLLISIONS_GUIDE.md section 2.

#pragma once

#include <stdint.h>

// 3D neighbor byte layout: 0b_00_zz_yy_xx, two bits per axis.
//   axis k (x=0,y=1,z=2): rightmost bit (1<<k*2)   = +1 neighbor filled
//                          leftmost bit  (1<<k*2+1) = -1 neighbor filled
// Special sentinels:
#define B3_VOXEL_INTERIOR_FACE_MASK 0b00111111u // all 6 neighbors filled
#define B3_VOXEL_EMPTY_FACE_MASK    0b01000000u // the cell itself is empty (bit 6, outside the 6-window)

// AxisMask bitflag over the 6 neighbor faces (matches the byte's low 6 bits).
typedef enum b3VoxelAxisMask
{
	b3_voxelX_POS = 1 << 0,
	b3_voxelX_NEG = 1 << 1,
	b3_voxelY_POS = 1 << 2,
	b3_voxelY_NEG = 1 << 3,
	b3_voxelZ_POS = 1 << 4,
	b3_voxelZ_NEG = 1 << 5,
	b3_voxelAllFaces = 0b00111111,
} b3VoxelAxisMask;

// Voxel classification by neighborhood (parry VoxelType).
typedef enum b3VoxelType
{
	b3_voxelTypeEmpty = 0,
	b3_voxelTypeVertex,
	b3_voxelTypeEdge, // 3D only
	b3_voxelTypeFace,
	b3_voxelTypeInterior,
} b3VoxelType;

// OctantPattern: 3-bit value packed per octant in FACES_TO_OCTANT_MASKS.
// Tells what feature a voxel exposes in a given octant direction.
typedef enum b3VoxelOctantPattern
{
	b3_voxelOctantInterior = 0,
	b3_voxelOctantVertex = 1,
	b3_voxelOctantEdgeX = 2,
	b3_voxelOctantEdgeY = 3,
	b3_voxelOctantEdgeZ = 4,
	b3_voxelOctantFaceX = 5,
	b3_voxelOctantFaceY = 6,
	b3_voxelOctantFaceZ = 7,
} b3VoxelOctantPattern;

// FACES_TO_VOXEL_TYPES[i] classifies a voxel from its 6-neighborhood byte i.
// 64 real cases + Interior (0x3F) + Empty (0x40) = 65 entries.
// Values are b3VoxelType: Vertex=1, Edge=2, Face=3, Interior=4, Empty=0.
static const uint8_t b3_facesToVoxelTypes[65] = {
	1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2, 2, 2, 2, 3, // 0..15
	1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2, 2, 2, 2, 3, // 16..31
	1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2, 2, 2, 2, 3, // 32..47
	2, 2, 2, 3, 2, 2, 2, 3, 2, 2, 2, 3, 3, 3, 3, 4, // 48..63
	0,												// 64 = Empty
};

// FACES_TO_FEATURE_MASKS[i]: a u16 whose meaning depends on the voxel type:
//   Vertex -> bit v => AABB corner v is convex
//   Edge   -> bit e => AABB edge e (EDGES_VERTEX_IDS) is convex
//   Face   -> bit f => AABB face f (FACES_VERTEX_IDS) is exposed
static const uint16_t b3_facesToFeatureMasks[65] = {
	0b11111111,		0b10011001,		0b1100110,		0b1010101,		0b110011,		0b10001,		0b100010,		0b10001,
	0b11001100,		0b10001000,		0b1000100,		0b1000100,		0b10101010,		0b10001000,		0b100010,		0b110000,
	0b1111,			0b1001,			0b110,			0b101,			0b11,			0b1,			0b10,			0b1,
	0b1100,			0b1000,			0b100,			0b100,			0b1010,			0b1000,			0b10,			0b100000,
	0b11110000,		0b10010000,		0b1100000,		0b1010000,		0b110000,		0b10000,		0b100000,		0b10000,
	0b11000000,		0b10000000,		0b1000000,		0b1000000,		0b10100000,		0b10000000,		0b100000,		0b10000,
	0b111100000000, 0b100100000000, 0b11000000000, 0b1100,			0b1100000000,	0b100000000,	0b1000000000,	0b1000,
	0b110000000000, 0b100000000000, 0b10000000000, 0b100,			0b11,			0b10,			0b1,			0b1111111111111111,
	0,
};

// FACES_TO_OCTANT_MASKS[i]: 8 octants x 3 bits = 24 bits. Each field is a VoxelOctantPattern
// telling what feature the voxel exposes approached from that octant. Priority vertex>edge>face
// is baked in — this is what prevents two adjacent voxels from both claiming a shared edge.
static const uint32_t b3_facesToOctantMasks[65] = {
	0b1001001001001001001001,	 0b1010010001001010010001,	 0b10001001010010001001010, 0b10010010010010010010010,
	0b11011001001011011001001,	 0b11111010001011111010001,	 0b111011001010111011001010, 0b111111010010111111010010,
	0b1001011011001001011011,	 0b1010111011001010111011,	 0b10001011111010001011111, 0b10010111111010010111111,
	0b11011011011011011011011,	 0b11111111011011111111011,	 0b111011011111111011011111, 0b111111111111111111111111,
	0b100100100100001001001001,	 0b100110110100001010010001, 0b110100100110010001001010, 0b110110110110010010010010,
	0b101101100100011011001001,	 0b101000110100011111010001, 0b101100110111011001010,	  0b110110111111010010,
	0b100100101101001001011011,	 0b100110000101001010111011, 0b110100101000010001011111, 0b110110000000010010111111,
	0b101101101101011011011011,	 0b101000000101011111111011, 0b101101000111011011111,	  0b111111111111,
	0b1001001001100100100100,	 0b1010010001100110110100,	  0b10001001010110100100110, 0b10010010010110110110110,
	0b11011001001101101100100,	 0b11111010001101000110100,	  0b111011001010000101100110, 0b111111010010000000110110,
	0b1001011011100100101101,	 0b1010111011100110000101,	  0b10001011111110100101000, 0b10010111111110110000000,
	0b11011011011101101101101,	 0b11111111011101000000101,	  0b111011011111000101101000, 0b111111111111000000000000,
	0b100100100100100100100100,	 0b100110110100100110110100,	 0b110100100110110100100110, 0b110110110110110110110110,
	0b101101100100101101100100,	 0b101000110100101000110100,	 0b101100110000101100110,	  0b110110000000110110,
	0b100100101101100100101101,	 0b100110000101100110000101,	 0b110100101000110100101000, 0b110110000000110110000000,
	0b101101101101101101101101,	 0b101000000101101000000101,	 0b101101000000101101000,	  0b0,
	0,
};

// AABB vertex/edge/face topology, matching parry's conventions (so the octant remap and feature
// masks line up). Vertices are indexed as a unit cube: bit0=x, bit1=y, bit2=z (8 corners).
// Face order: 0=+X 1=-X 2=+Y 3=-Y 4=+Z 5=-Z (matches VoxelAxisMask bit order).
static const uint8_t b3_aabbFacesVertexIds[6][4] = {
	{ 1, 2, 6, 5 }, // +X
	{ 0, 3, 7, 4 }, // -X
	{ 2, 3, 7, 6 }, // +Y
	{ 1, 0, 4, 5 }, // -Y
	{ 4, 5, 6, 7 }, // +Z
	{ 0, 1, 2, 3 }, // -Z
};

// Edge order and which axis each edge runs along (parry EDGE_AXIS).
static const uint8_t b3_aabbEdgesVertexIds[12][2] = {
	{ 0, 1 }, { 1, 2 }, { 3, 2 }, { 0, 3 }, { 4, 5 }, { 5, 6 },
	{ 7, 6 }, { 4, 7 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 },
};
// EDGE_AXIS[eid]: which VoxelOctantPattern (EDGE_X/Y/Z) edge eid exposes.
static const uint8_t b3_aabbEdgeAxis[12] = {
	b3_voxelOctantEdgeX, b3_voxelOctantEdgeY, b3_voxelOctantEdgeX, b3_voxelOctantEdgeY,
	b3_voxelOctantEdgeX, b3_voxelOctantEdgeY, b3_voxelOctantEdgeX, b3_voxelOctantEdgeY,
	b3_voxelOctantEdgeZ, b3_voxelOctantEdgeZ, b3_voxelOctantEdgeZ, b3_voxelOctantEdgeZ,
};
// FACE_NORMALS[fid]: which VoxelOctantPattern (FACE_X/Y/Z) face fid exposes.
static const uint8_t b3_aabbFaceNormal[6] = {
	b3_voxelOctantFaceX, b3_voxelOctantFaceX, b3_voxelOctantFaceY,
	b3_voxelOctantFaceY, b3_voxelOctantFaceZ, b3_voxelOctantFaceZ,
};

// Remap from an easy octant key (bit0=x>=0, bit1=y>=0, bit2=z>=0) to the AABB vertex key used
// by the tables above. Must match parry's AABB_OCTANT_KEYS.
static const uint8_t b3_aabbOctantKeys[8] = { 0, 1, 3, 2, 4, 5, 7, 6 };
