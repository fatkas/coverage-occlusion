
#include "rasterizer.h"
#include "rasterizer_math.h"

__forceinline bool Rasterizer::draw_scanlines(Tile& tile, int& xs1_global, int& xs2_global, int y1, int y2, int xa1, int xa2, const RESTRICT vec4i_t* masks, RESTRICT uint32_t* flag)
{
    uint32_t mask = 0;
    vec4i_t full_span = m_full_span;
    int xs1 = xs1_global, xs2 = xs2_global, bits = g_fixed_point_bits;
    for (int scanline = y1; scanline < y2; ++scanline)
    {
        int xb = xs1 >> bits;
        int xe = xs2 >> bits;
        xs1 += xa1;
        xs2 += xa2;

        assert(scanline >= 0);
        assert(scanline < Tile::g_tile_height);

        vec4i_t span = VecIntXor(masks[xb], masks[xe]);
        if (flag == nullptr)
        {
            tile.m_frame_buffer[scanline] = VecIntOr(tile.m_frame_buffer[scanline], span);
            uint32_t bit = VecIntMask(VecIntCmpEqual(VecIntAnd(tile.m_frame_buffer[scanline], full_span), full_span)) == 65535;
            mask |= bit << scanline;
        }
        else
        {
            if (VecIntMask(VecIntCmpEqual(VecIntAnd(tile.m_frame_buffer[scanline], span), span)) != 65535)
            {
                *flag = 1;
                return true;
            }
        }
    }
    xs1_global = xs1;
    xs2_global = xs2;
    tile.m_mask |= mask;
    return tile.m_mask == ~0u;
}

__forceinline void Rasterizer::draw_4triangles(Tile& tile, const TriangleType& RESTRICT tri, uint32_t* RESTRICT flag, const vec4i_t* RESTRICT masks)
{
    vec4_t inv_quantizer = m_inv_fixed_point;
#if USE_PACKED_TRIANGLES
    vec4i_t x0x1 = VecIntLoad(tri.x0);
    vec4i_t x2y0 = VecIntLoad(tri.x2);
    vec4i_t y1y2 = VecIntLoad(tri.y1);

    vec4i_t lo = VecIntUnpackLo(x0x1, VecIntZero());
    vec4i_t hi = VecIntUnpackHi(x0x1, VecIntZero());
    vec4_t x0 = VecInt2Float(lo);
    vec4_t x1 = VecInt2Float(hi);
    lo = VecIntUnpackLo(x2y0, VecIntZero());
    hi = VecIntUnpackHi(x2y0, VecIntZero());
    vec4_t x2 = VecInt2Float(lo);
    vec4_t y0 = VecMul(VecInt2Float(hi), inv_quantizer);
    lo = VecIntUnpackLo(y1y2, VecIntZero());
    hi = VecIntUnpackHi(y1y2, VecIntZero());
    vec4_t y1 = VecMul(VecInt2Float(lo), inv_quantizer);
    vec4_t y2 = VecMul(VecInt2Float(hi), inv_quantizer);

    vec4_t vx0 = x0, vx1 = x1, vx2 = x2;
    vec4_t vy0 = y0, vy1 = y1, vy2 = y2;
#else
    vec4_t vx0 = tri.x0, vx1 = tri.x1, vx2 = tri.x2;
    vec4_t vy0 = VecMul(tri.y0, inv_quantizer), vy1 = VecMul(tri.y1, inv_quantizer), vy2 = VecMul(tri.y2, inv_quantizer);
#endif

    ALIGN16 int iy0[4], iy1[4], iy2[4], ix0[4], ix1[4], ix2[4], dx1[4], dx2[4], dx3[4];

    vec4_t vdx1 = VecMul(VecSub(vx2, vx0), VecRcp(VecSub(vy2, vy0)));
    vec4_t vdx2 = VecMul(VecSub(vx1, vx0), VecRcp(VecSub(vy1, vy0)));
    vec4_t vdx3 = VecMul(VecSub(vx2, vx1), VecRcp(VecSub(vy2, vy1)));
    VecIntStore(dx1, VecFloat2Int(vdx1));
    VecIntStore(dx2, VecFloat2Int(vdx2));
    VecIntStore(dx3, VecFloat2Int(vdx3));

    // adjust y to be inside tile bounds
    vy0 = VecSub(vy0, tile.m_y);
    vy1 = VecSub(vy1, tile.m_y);
    vy2 = VecSub(vy2, tile.m_y);
    vec4_t dy0 = VecMax(VecZero(), -vy0);
    vec4_t dy1 = VecMax(VecZero(), -vy1);
    vy0 = VecMin(VecMax(vy0, VecZero()), m_tile_height_v);
    vy1 = VecMin(VecMax(vy1, VecZero()), m_tile_height_v);
    vy2 = VecMin(VecMax(vy2, VecZero()), m_tile_height_v);

    VecIntStore(iy0, VecFloat2Int(vy0));
    VecIntStore(iy1, VecFloat2Int(vy1));
    VecIntStore(iy2, VecFloat2Int(vy2));
    VecIntStore(ix0, VecFloat2Int(VecMad(vdx1, dy0, vx0)));
    VecIntStore(ix1, VecFloat2Int(VecMad(vdx2, dy0, vx0)));
    VecIntStore(ix2, VecFloat2Int(VecMad(vdx3, dy1, vx1)));

    bool skip_full = m_skip_full;
    for (size_t i = 0; i < 4; ++i)
    {
        assert(iy0[i] <= 32);
        assert(iy2[i] <= 32);
        uint32_t span = (0xffffffff << iy0[i]) & (0xffffffff >> (Tile::g_tile_height - iy2[i]));
        if (skip_full && (tile.m_mask & span) == span)
        {
#if USE_STATS
            ++tile.m_triangles_skipped;
#endif
            continue;
        }

#if USE_STATS
        tile.m_triangles_drawn_total++;
        if (flag == nullptr)
            tile.m_triangles_drawn_occluder_total++;
        else
            tile.m_triangles_drawn_occludee_total++;
#endif

        int xs1 = ix0[i], xs2 = ix1[i], xs3 = ix2[i];
        if (draw_scanlines(tile, xs1, xs2, iy0[i], iy1[i], dx1[i], dx2[i], masks, flag))
            return;
        if (draw_scanlines(tile, xs1, xs3, iy1[i], iy2[i], dx1[i], dx3[i], masks, flag))
            return;
    }
}

void Rasterizer::draw_triangles(uint32_t tile_index)
{
    auto & tile = m_tiles[tile_index];
    auto & tile_indices = m_data.data.tiles[tile_index];

    if (tile_indices.triangle_index_count == 0)
        return;

    assert(tile_indices.triangle_index_count <= m_data.data.triangle_count);
    const SortKey* tri = tile_indices.triangle_index_data, *tri_end = tri + tile_indices.triangle_index_count;
    const TriangleType* tri_data = m_data.data.triangle_data;
    uint32_t** flags = m_flags.data();
    vec4i_t* masks = m_masks.data() + tile.m_x*g_max_masks_per_tile;
    bool skip_full = m_skip_full;
    while (tri != tri_end)
    {
        assert(tri->index < m_data.data.triangle_count);
        auto & key = *tri++;
        auto & triangle = tri_data[key.index];

        uint32_t* flag = key.flag ? flags[key.flag] : nullptr;
        if (key.flag && *flag)
            continue;
        draw_4triangles(tile, triangle, flag, masks);
        if (skip_full && tile.m_mask == ~0u)
        {
            break;
        }
    }
#if USE_STATS
    if (m_mt)
    {
        atomic_add(m_triangles_drawn_total, tile.m_triangles_drawn_total);
        atomic_add(m_triangles_drawn_occluder_total, tile.m_triangles_drawn_occluder_total);
        atomic_add(m_triangles_drawn_occludee_total, tile.m_triangles_drawn_occludee_total);
        atomic_add(m_triangles_skipped, tile.m_triangles_skipped);
    }
    else
    {
        m_triangles_drawn_total += tile.m_triangles_drawn_total;
        m_triangles_drawn_occluder_total += tile.m_triangles_drawn_occluder_total;
        m_triangles_drawn_occludee_total += tile.m_triangles_drawn_occludee_total;
        m_triangles_skipped += tile.m_triangles_skipped;
    }
#endif
}

void Rasterizer::draw_triangles()
{
    for (uint32_t i = 0; i < g_width*g_height; ++i)
        draw_triangles(i);
}
