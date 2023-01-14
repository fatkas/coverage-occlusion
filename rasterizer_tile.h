#pragma once

#include <assert.h>
#include "rasterizer_math.h"

#include <tinystl/allocator.h>
#include <tinystl/vector.h>

namespace stl = tinystl;

// 112bytes, 28 bytes per priangle
struct ALIGN16 Triangle
{
    // x1, x2
    // dx1, dx2, dx3
    // y123
	vec4_t x0, x1, x2;
    vec4_t y0, y1, y2;
    uint32_t flag;
    uint32_t mask;
};
static_assert(sizeof(Triangle) == 112, "this is important");

// 52 bytes per 4 triangles, 13 bytes per triangle
struct TrianglePacked
{
    // 24 bytes
    uint16_t x0[4];
    uint16_t x1[4];
    uint16_t x2[4];

    // 24 bytes
    uint16_t y0[4];
    uint16_t y1[4];
    uint16_t y2[4];

    uint32_t flag : 28;
    uint32_t mask : 4;
};
static_assert(sizeof(TrianglePacked) == 52, "this is important");

struct ALIGN16 Tile
{
    static constexpr int    g_tile_height = 32;
    static constexpr int    g_tile_width = 128;

    vec4i_t                 m_frame_buffer[g_tile_height];
    uint64_t                m_mask = 0;
    uint32_t                m_triangles_drawn_total = 0;
    uint32_t                m_triangles_drawn_occluder_total = 0;
    uint32_t                m_triangles_drawn_occludee_total = 0;
    uint32_t                m_triangles_skipped = 0;
    int                     m_x = 0;
    int                     m_y = 0;

    Tile(int x, int y);

    void clear()
    {
        memset(m_frame_buffer, 0, sizeof(m_frame_buffer));
        m_mask = 0;

        m_triangles_drawn_total = 0;
        m_triangles_drawn_occluder_total = 0;
        m_triangles_drawn_occludee_total = 0;
        m_triangles_skipped = 0;
    }
};
