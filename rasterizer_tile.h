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
    uint32_t* flag;
    uint32_t z;
    uint32_t mask;
};

struct ALIGN16 Tile
{
    static constexpr int    g_tile_height = 32;
    static constexpr int    g_tile_width = 128;

    vec4i_t                 m_frame_buffer[g_tile_height];
    stl::vector<vec4i_t>    m_shifts;
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
