#pragma once

#include <assert.h>
#include "rasterizer_math.h"

#include <tinystl/allocator.h>
#include <tinystl/vector.h>

namespace stl = tinystl;

// 112bytes, 28 bytes per priangle
struct _MM_ALIGN16 Triangle
{
    // x1, x2
    // dx1, dx2, dx3
    // y123
	__m128 x0, x1, x2;
	__m128 y0, y1, y2;
    int* flag;
    uint32_t z;
    uint32_t mask;
};

// 96 bytes
struct TrianglePacked
{
    uint16_t x0[4]; // 10:6
    uint16_t x1[4]; // 10:6

    uint16_t y0[4]; // 10:6
    uint16_t y1[4]; // 10:6
    uint16_t y2[4]; // 10:6

    int32_t dx0[4]; // 15:16
    int32_t dx1[4]; // 15:16
    int32_t dx2[4]; // 15:16
};

struct _MM_ALIGN16 Tile
{
    static constexpr int        g_tile_height = 32;
    static constexpr int        g_tile_width = 128;
    static constexpr int        g_max_triangles = 64 * 1024;

    __m128i                 m_frame_buffer[g_tile_height];
    stl::vector<uint64_t>   m_triangles;
    uint32_t                m_triangle_count = 0;
    stl::vector<__m128i>    m_shifts;
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
        m_triangle_count = 0;

        m_mask = 0;

        m_triangles_drawn_total = 0;
        m_triangles_drawn_occluder_total = 0;
        m_triangles_drawn_occludee_total = 0;
        m_triangles_skipped = 0;
    }
};