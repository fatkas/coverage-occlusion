
#include "rasterizer.h"

template <bool is_occluder> __forceinline
bool Rasterizer::draw_scanlines(Tile& tile, int& xs1, int& xs2, int y1, int y2, int xa1, int xa2, const vec4i_t* masks, uint32_t* flag)
{
    assert((is_occluder && !flag) || (flag && !is_occluder));
    for (int scanline = y1; scanline < y2; ++scanline)
    {
        int xb = xs1 >> 16;
        int xe = xs2 >> 16;
        xs1 += xa1;
        xs2 += xa2;

        assert(scanline >= 0);
        assert(scanline < Tile::g_tile_height);

        vec4i_t span = VecIntXor(masks[xb], masks[xe]);
        if (is_occluder)
        {
            tile.m_frame_buffer[scanline] = VecIntOr(tile.m_frame_buffer[scanline], span);
            uint64_t bit = VecIntMask(VecIntCmpEqual(VecIntAnd(tile.m_frame_buffer[scanline], m_full_span), m_full_span)) == 65535;
            tile.m_mask |= bit << scanline;
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
    return false;
}

template < bool is_occluder >
__forceinline void Rasterizer::draw_4triangles(Tile& tile, const Triangle& tri)
{
    vec4_t vx0 = tri.x0, vx1 = tri.x1, vx2 = tri.x2;
    vec4_t vy0 = tri.y0, vy1 = tri.y1, vy2 = tri.y2;

    ALIGN16 int iy0[4], iy1[4], iy2[4], ix0[4], ix1[4], ix2[4], dx1[4], dx2[4], dx3[4];

    vec4_t vdx1 = VecMul(VecSub(vx2, vx0), VecRcp(VecSub(vy2, vy0)));
    vec4_t vdx2 = VecMul(VecSub(vx1, vx0), VecRcp(VecSub(vy1, vy0)));
    vec4_t vdx3 = VecMul(VecSub(vx2, vx1), VecRcp(VecSub(vy2, vy1)));
    VecIntStore(dx1, VecFloat2Int(vdx1));
    VecIntStore(dx2, VecFloat2Int(vdx2));
    VecIntStore(dx3, VecFloat2Int(vdx3));

    // adjust y to be inside tile bounds
    vy0 = VecSub(vy0, Vector4(tile.m_y*Tile::g_tile_height));
    vy1 = VecSub(vy1, Vector4(tile.m_y*Tile::g_tile_height));
    vy2 = VecSub(vy2, Vector4(tile.m_y*Tile::g_tile_height));
    vec4_t dy0 = VecMax(VecZero(), -vy0);
    vec4_t dy1 = VecMax(VecZero(), -vy1);
    vy0 = VecMin(VecMax(vy0, VecZero()), Vector4(Tile::g_tile_height));
    vy1 = VecMin(VecMax(vy1, VecZero()), Vector4(Tile::g_tile_height));
    vy2 = VecMin(VecMax(vy2, VecZero()), Vector4(Tile::g_tile_height));

    VecIntStore(iy0, VecFloat2Int(vy0));
    VecIntStore(iy1, VecFloat2Int(vy1));
    VecIntStore(iy2, VecFloat2Int(vy2));
    VecIntStore(ix0, VecFloat2Int(VecMad(vdx1, dy0, vx0)));
    VecIntStore(ix1, VecFloat2Int(VecMad(vdx2, dy0, vx0)));
    VecIntStore(ix2, VecFloat2Int(VecMad(vdx3, dy1, vx1)));

    for (size_t i = 0, mask = 1; i < 4; ++i, mask <<= 1)
    {
        if ((tri.mask & mask ) == 0)
            continue;

        assert(iy0[i] <= 32);
        assert(iy2[i] <= 32);
        uint32_t span0 = 0xffffffff << iy0[i];
        uint32_t span1 = 0xffffffff >> (32-iy2[i]);
        uint32_t span = span0 & span1;
        if (m_skip_full && (tile.m_mask & span) == span)
        {
            ++tile.m_triangles_skipped;
            continue;
        }

        tile.m_triangles_drawn_total++;
        if (is_occluder)
            tile.m_triangles_drawn_occluder_total++;
        else
            tile.m_triangles_drawn_occludee_total++;

        int xs1 = ix0[i], xs2 = ix1[i], xs3 = ix2[i];
        if (draw_scanlines<is_occluder>(tile, xs1, xs2, iy0[i], iy1[i], dx1[i], dx2[i], tile.m_shifts.data(), tri.flag))
            return;
        if (draw_scanlines<is_occluder>(tile, xs1, xs3, iy1[i], iy2[i], dx1[i], dx3[i], tile.m_shifts.data(), tri.flag))
            return;
    }
}

void Rasterizer::draw_triangles(uint32_t tile_index)
{
    auto & tile = m_tiles[tile_index];
    auto & tile_indices = m_data.data.tiles[tile_index];

    assert(tile_indices.triangle_index_count <= m_data.data.triangle_count);
    const SortKey* tri = tile_indices.triangle_index_data, *tri_end = tri + tile_indices.triangle_index_count;
    const Triangle* tri_data = m_data.data.triangle_data;
    while (tri != tri_end)
    {
        assert(tri->index < m_data.data.triangle_count);
        auto & key = *tri++;
        auto & triangle = tri_data[key.index];
        if (triangle.flag)
        {
            if (m_skip_full && *triangle.flag)
                tile.m_triangles_skipped += bx::uint32_cntbits(triangle.mask);
            else
                draw_4triangles<false>(tile, triangle);
        }
        else
        {
            draw_4triangles<true>(tile, triangle);
            if (m_skip_full && tile.m_mask == ~0u)
            {
                while (tri != tri_end)
                {
                    auto & lkey = *tri++;
                    tile.m_triangles_skipped += bx::uint32_cntbits(tri_data[lkey.index].mask);
                }
                break;
            }
        }
    }
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
}

void Rasterizer::draw_triangles()
{
    for (uint32_t i = 0; i < g_width*g_height; ++i)
        draw_triangles(i);
}
