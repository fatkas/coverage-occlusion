
#include "rasterizer.h"
#include <stdio.h>

#include <algorithm>

inline static void sort(vec4_t& RESTRICT A0, vec4_t& RESTRICT A1, vec4_t& RESTRICT B0, vec4_t& RESTRICT B1)
{
    vec4_t mask = VecCmpLe(B0, B1);
    vec4_t sx = VecAdd(A0, A1);
    A0 = VecOr(VecAnd(mask, A0), VecAndNot(mask, A1));
    A1 = VecSub(sx, A0);

    vec4_t sy = VecAdd(B0, B1);
    B0 = VecOr(VecAnd(mask, B0), VecAndNot(mask, B1));
    B1 = VecSub(sy, B0);
}

inline static void ExtractMatrix(const Matrix& RESTRICT m, vec4_t* RESTRICT matrix)
{
    #define EXTRACT(line) matrix[line*4+0] = VecShuffle( m.r[line], m.r[line], VecShuffleMask(0, 0, 0, 0));  \
                        matrix[line*4+1] = VecShuffle( m.r[line], m.r[line], VecShuffleMask(1, 1, 1, 1)); \
                        matrix[line*4+2] = VecShuffle( m.r[line], m.r[line], VecShuffleMask(2, 2, 2, 2)); \
                        matrix[line*4+3] = VecShuffle( m.r[line], m.r[line], VecShuffleMask(3, 3, 3, 3));
    EXTRACT(0);
    EXTRACT(1);
    EXTRACT(2);
    EXTRACT(3);

    #undef EXTRACT
}

template< int v, bool use_plane >
__forceinline static vec4_t intersectLineZ(vec4_t a, vec4_t b, vec4_t plane)
{
    //  t = (a.x - a.w) / (b.w - a.w - b.x + a.x);
    vec4_t bz = VecShuffle(b, b, VecShuffleMask(v, v, v, v));
    vec4_t az = VecShuffle(a, a, VecShuffleMask(v, v, v, v));
    vec4_t bw = VecMul(plane, VecShuffle(b, b, VecShuffleMask(3, 3, 3, 3)));
    vec4_t aw = VecMul(plane, VecShuffle(a, a, VecShuffleMask(3, 3, 3, 3)));
    return VecAdd(b, VecMul(VecSub(a, b), use_plane ? VecMul(VecSub(bw, bz), VecRcp(VecAdd(VecSub(az, bz), VecSub(bw, aw)))) : VecMul(bz, VecRcp(VecSub(bz, az)))));
}

template< int vertex_component, bool cmp_func >
__forceinline static uint32_t clip_triangle(vec4_t* RESTRICT input, uint32_t& vertex_count, const uint16_t* RESTRICT indices, uint32_t index_count,
                                            const uint8_t* RESTRICT normal_masks, uint32_t normal_mask, uint16_t* RESTRICT output, vec4_t plane)
{
    uint32_t output_indices = 0;
    for (uint32_t i = 0; i < index_count / 3; ++i)
    {
#if USE_NORMAL_MASKS
        if (normal_masks && (normal_masks[i/3] & normal_mask) == 0)
        {
            indices += 3;
            continue;
        }
#endif
        uint16_t i0 = *indices++;
        uint16_t i1 = *indices++;
        uint16_t i2 = *indices++;
        vec4_t v0 = input[i0];
        vec4_t v1 = input[i1];
        vec4_t v2 = input[i2];
        vec4_t tmp = VecShuffle(v0, v1, VecShuffleMask(vertex_component, 0, vertex_component, 0));
        vec4_t w = VecShuffle(tmp, v2, VecShuffleMask(0, 2, vertex_component, vertex_component));

        int mask = cmp_func ? VecMask(VecCmpGt(w, VecMul(VecShuffle(VecShuffle(v0, v1, VecShuffleMask(3, 3, 3, 3)), v2, VecShuffleMask(0, 2, 3, 3)), plane))) & 7
                            : VecMask(VecCmpLt(w, VecZero())) & 7;

        switch (mask)
        {
        case 7:
            break;
        case 0:
            output[output_indices + 0] = i0;
            output[output_indices + 1] = i1;
            output[output_indices + 2] = i2;
            output_indices += 3;
            break;
        case 1:
        case 2:
        case 4:
            {
                #define SHUFFLE(type,d0,d1,d2,s0,s1,s2) {type tv = s0; d0 = s1; d1 = s2; d2 = tv; }
                if (mask & 1)
                {
                    SHUFFLE(vec4_t, v0, v1, v2, v0, v1, v2);
                    SHUFFLE(uint16_t, i0, i1, i2, i0, i1, i2);
                }
                else if (mask & 2)
                {
                    SHUFFLE(vec4_t, v0, v2, v1, v0, v2, v1);
                    SHUFFLE(uint16_t, i0, i2, i1, i0, i2, i1);
                }
                #undef SHUFFLE

                uint16_t ba = vertex_count;
                input[vertex_count++] = intersectLineZ<vertex_component, cmp_func>(v2, v0, plane);
                uint16_t bb = vertex_count;
                input[vertex_count++] = intersectLineZ<vertex_component, cmp_func>(v2, v1, plane);

                output[output_indices + 0] = i0;
                output[output_indices + 1] = i1;
                output[output_indices + 2] = bb;
                output[output_indices + 3] = bb;
                output[output_indices + 4] = ba;
                output[output_indices + 5] = i0;
                output_indices += 6;
            }
            break;
        default:
            {
                vec4_t in = (mask & 1) == 0 ? v0 : ((mask & 4) == 0 ? v2 : v1);
                if (mask & 1)
                {
                    output[output_indices + 0] = vertex_count;
                    input[vertex_count++] = intersectLineZ< vertex_component, cmp_func >(v0, in, plane);
                }
                else
                    output[output_indices + 0] = i0;
                if (mask & 2)
                {
                    output[output_indices + 1] = vertex_count;
                    input[vertex_count++] = intersectLineZ< vertex_component, cmp_func >(v1, in, plane);
                }
                else
                    output[output_indices + 1] = i1;
                if (mask & 4)
                {
                    output[output_indices + 2] = vertex_count;
                    input[vertex_count++] = intersectLineZ< vertex_component, cmp_func >(v2, in, plane);
                }
                else
                    output[output_indices + 2] = i2;
                output_indices += 3;
            }
            break;
        }
    }
    return output_indices;
}

__forceinline static uint32_t clip_triangles(vec4_t* RESTRICT vertices, uint32_t& vertex_count, const uint16_t* RESTRICT indices, uint32_t index_count,
                                             const uint8_t* RESTRICT normal_masks, uint32_t normal_mask, uint16_t* RESTRICT output_indices)
{
    vec4_t g_total_width_v = Vector4(Rasterizer::g_total_width);
    vec4_t g_total_height_v = Vector4(Rasterizer::g_total_height);

    int count = 4;
    uint16_t input_array[1024], output_array[1024];
    count = clip_triangle<0, true>(vertices, vertex_count, indices, index_count, normal_masks, normal_mask, input_array, g_total_width_v); // x > 1280
    count = clip_triangle<0, false>(vertices, vertex_count, input_array, count, nullptr, 0, output_array, VecZero()); // x < 0
    count = clip_triangle<1, false>(vertices, vertex_count, output_array, count, nullptr, 0, input_array, VecZero()); // y < 0
    count = clip_triangle<1, true>(vertices, vertex_count, input_array, count, nullptr, 0, output_array, g_total_height_v); // y > 720
    count = clip_triangle<2, false>(vertices, vertex_count, output_array, count, nullptr, 0, output_indices, VecZero()); // z < 0
    return count;
}

__forceinline void Rasterizer::push_4triangles(TrianagleData& RESTRICT data, uint32_t flag, int* RESTRICT bounds_array,
                                               vec2_t RESTRICT v0[4], vec2_t RESTRICT v1[4], vec2_t RESTRICT v2[4], uint32_t w)
{
    Triangle t;
    vec4_t p0_0 = VecLoad(v0 + 0); // x0_0 y0_0 x0_1 y0_1
    vec4_t p0_1 = VecLoad(v0 + 2); // x0_2 y0_2 x0_3 y0_3
    t.x0 = VecShuffle(p0_0, p0_1, VecShuffleMask(0, 2, 0, 2));
    t.y0 = VecShuffle(p0_0, p0_1, VecShuffleMask(1, 3, 1, 3));
    p0_0 = VecLoad(v1 + 0); // x1_0 y1_0 x1_1 y1_1
    p0_1 = VecLoad(v1 + 2); // x1_2 y1_2 x1_3 y1_3
    t.x1 = VecShuffle(p0_0, p0_1, VecShuffleMask(0, 2, 0, 2));
    t.y1 = VecShuffle(p0_0, p0_1, VecShuffleMask(1, 3, 1, 3));
    p0_0 = VecLoad(v2 + 0); // x2_0 y2_0 x2_1 y2_1
    p0_1 = VecLoad(v2 + 2); // x2_2 y2_2 x2_3 y2_3
    t.x2 = VecShuffle(p0_0, p0_1, VecShuffleMask(0, 2, 0, 2));
    t.y2 = VecShuffle(p0_0, p0_1, VecShuffleMask(1, 3, 1, 3));

    sort(t.x0, t.x1, t.y0, t.y1);
    sort(t.x1, t.x2, t.y1, t.y2);
    sort(t.x0, t.x1, t.y0, t.y1);

    assert(data.triangle_count < data.triangles.size());
    TriangleType & tri = data.triangle_data[data.triangle_count];

#if USE_PACKED_TRIANGLES
    vec4i_t x0x1 = VecIntPack16(VecFloat2Int(t.x0), VecFloat2Int(t.x1));
    vec4i_t x2y0 = VecIntPack16(VecFloat2Int(t.x2), VecFloat2Int(t.y0));
    vec4i_t y1y2 = VecIntPack16(VecFloat2Int(t.y1), VecFloat2Int(t.y2));

    VecIntStore(tri.x0, x0x1);
    VecIntStore(tri.x2, x2y0);
    VecIntStore(tri.y1, y1y2);
#else
    tri = t;
#endif

    SortKey key;
    key.z = w;
    key.index = data.triangle_count++;
    key.flag = flag;
    for (int yy = bounds_array[1]; yy < bounds_array[3]; ++yy)
        for (int xx = bounds_array[0]; xx < bounds_array[2]; ++xx)
        {
            uint32_t tile_index = xx + yy*g_width;
            assert(xx < g_width);
            assert(yy < g_height);
            auto & tile = data.tiles[tile_index];
            assert(tile.triangle_index_count < tile.triangle_indices.size());
            tile.triangle_index_data[tile.triangle_index_count++] = key;
        }
}

// https://guru.multimedia.cx/avoiding-branchesifconditionals/
__forceinline int32_t max4(int32_t a, int32_t b, int32_t c, int32_t d)
{
    int32_t t0 = a-((a-b)&((a-b)>>31));
    int32_t t1 = c-((c-d)&((c-d)>>31));
    return t0-((t0-t1)&((t0-t1)>>31));
}

__forceinline int32_t min4(int32_t a, int32_t b, int32_t c, int32_t d)
{
    int32_t t0 = b+((a-b)&((a-b)>>31));
    int32_t t1 = d+((c-d)&((c-d)>>31));
    return t1+((t0-t1)&((t0-t1)>>31));
}

__forceinline void Rasterizer::push_triangle_batched(TrianagleData& RESTRICT data, uint32_t flag,
                                                     const vec2_t* RESTRICT src, const uint32_t* RESTRICT w, const uint16_t* RESTRICT indices, int count,
                                                     const uint8_t* RESTRICT normal_masks, uint32_t normal_mask,
                                                     int* RESTRICT bounds_array, bool select_tiles)
{
    ALIGN16 vec2_t v0[5];
    ALIGN16 vec2_t v1[5];
    ALIGN16 vec2_t v2[5];
    vec4_t v_min = m_total_size, v_max = VecZero();
    vec4_t inv_quantizer = m_inv_fixed_point;

    uint32_t tris = 0, max_group_w = 65535*64, group_w = flag ? max_group_w : 0;
    for (int i = 0; i < count; i += 3)
    {
        if (normal_masks && (normal_masks[i/3] & normal_mask) == 0)
        {
#if USE_STATS
            if (!m_mt)
                m_triangles_backface++;
#endif
            indices += 3;
            continue;
        }

        int i0 = *indices++;
        int i1 = *indices++;
        int i2 = *indices++;
        vec4_t p0 = VecLoadU(&src[i0].x);
        vec4_t p1 = VecLoadU(&src[i1].x);
        vec4_t p2 = VecLoadU(&src[i2].x);
        if (select_tiles)
        {
            v_min = VecMin(VecMin(v_min, p0), VecMin(p1, p2));
            v_max = VecMax(VecMax(v_max, p0), VecMax(p1, p2));
        }
        VecStoreU(&v0[tris].x, p0);
        VecStoreU(&v1[tris].x, p1);
        VecStoreU(&v2[tris].x, p2);
        tris++;

        if (flag)
            group_w = min4(w[i0], w[i1], w[i2], group_w); // std::min(std::min(w[i0], w[i1]), std::min(group_w, w[i2]));
        else
            group_w = max4(w[i0], w[i1], w[i2], group_w);

        if (tris == 4)
        {
            if (select_tiles)
                VecIntStore(bounds_array, VecFloat2Int(get_tile_bounds(VecMul(VecMoveLH(v_min, v_max), inv_quantizer))));
            push_4triangles(data, flag, bounds_array, v0, v1, v2, group_w);
            tris = 0;
            group_w = flag ? max_group_w : 0;
            v_min = m_total_size, v_max = VecZero();
        }
    }
    if (tris)
    {
        for (int i = tris; i < 4; ++i)
            v0[i] = v0[tris-1], v1[i] = v1[tris - 1], v2[i] = v2[tris-1];
        if (select_tiles)
            VecIntStore(bounds_array, VecFloat2Int(get_tile_bounds(VecMul(VecMoveLH(v_min, v_max), inv_quantizer))));
        push_4triangles(data, flag, bounds_array, v0, v1, v2, group_w);
    }
}

bool Rasterizer::occlude_object(const vec4_t* RESTRICT m, vec4_t v_min, vec4_t v_max, int* RESTRICT bounds_array)
{
    vec4_t g_total_width_v = Vector4(g_total_width);
    vec4_t g_total_height_v = Vector4(g_total_height);

    vec4_t pt[4];
    vec4_t vTmp = VecUnpackLo(v_min, v_max);                // x, X, y, Y
    pt[0] = VecShuffle(v_min, v_max, VecShuffleMask(0, 0, 0, 0)); // xxXX
    pt[1] = VecShuffle(vTmp, vTmp, VecShuffleMask(2, 3, 2, 3)); // yYyY
    pt[2] = VecShuffle(v_min, v_min, VecShuffleMask(2, 2, 2, 2)); // zzzz
    pt[3] = VecShuffle(v_max, v_max, VecShuffleMask(2, 2, 2, 2)); // ZZZZ

    vec4_t xxxx0 = VecMad(m[8], pt[2], VecMad(m[4], pt[1], VecMad(m[0], pt[0], m[12])));
    vec4_t yyyy0 = VecMad(m[9], pt[2], VecMad(m[5], pt[1], VecMad(m[1], pt[0], m[13])));
    vec4_t zzzz0 = VecMad(m[10], pt[2], VecMad(m[6], pt[1], VecMad(m[2], pt[0], m[14])));
    vec4_t wwww0 = VecMad(m[11], pt[2], VecMad(m[7], pt[1], VecMad(m[3], pt[0], m[15])));

    vec4_t xxxx1 = VecMad(m[8], pt[3], VecMad(m[4], pt[1], VecMad(m[0], pt[0], m[12])));
    vec4_t yyyy1 = VecMad(m[9], pt[3], VecMad(m[5], pt[1], VecMad(m[1], pt[0], m[13])));
    vec4_t zzzz1 = VecMad(m[10], pt[3], VecMad(m[6], pt[1], VecMad(m[2], pt[0], m[14])));
    vec4_t wwww1 = VecMad(m[11], pt[3], VecMad(m[7], pt[1], VecMad(m[3], pt[0], m[15])));

    vec4_t v_mask00 = VecAnd(VecCmpGt(xxxx0, VecZero()), VecCmpGt(xxxx1, VecZero()));
    vec4_t v_mask01 = VecAnd(VecCmpGt(yyyy0, VecZero()), VecCmpGt(yyyy1, VecZero()));
    vec4_t v_mask10 = VecAnd(VecCmpLt(xxxx0, VecMul(wwww0, g_total_width_v)), VecCmpLt(xxxx1, VecMul(wwww1, g_total_width_v)));
    vec4_t v_mask11 = VecAnd(VecCmpLt(yyyy0, VecMul(wwww0, g_total_height_v)), VecCmpLt(yyyy1, VecMul(wwww1, g_total_height_v)));

    vec4_t v_mask0 = VecAnd(v_mask00, v_mask10);
    vec4_t v_mask1 = VecAnd(v_mask01, v_mask11);
    int mask = VecMask(VecAnd(v_mask0, v_mask1));

    bool intersect_near = VecMask(VecAnd(VecCmpGt(zzzz0, VecZero()), VecCmpGt(zzzz1, VecZero()))) != 15;

    vec4_t x_min, x_max, y_min, y_max;
    if (intersect_near == false)
    {
        vec4_t x0 = VecMul(xxxx0, VecRcp(wwww0));
        vec4_t y0 = VecMul(yyyy0, VecRcp(wwww0));
        vec4_t x1 = VecMul(xxxx1, VecRcp(wwww1));
        vec4_t y1 = VecMul(yyyy1, VecRcp(wwww1));

        x_min = VecMin(x0, x1);
        x_max = VecMax(x0, x1);
        y_min = VecMin(y0, y1);
        y_max = VecMax(y0, y1);
    }
    else
    {
#define INTERSECT_EDGE(a,b,c) VecMad(VecSub(a, b), t, b)
        vec4_t xxxx0_1 = VecShuffle(xxxx0, xxxx0, VecShuffleMask(1, 3, 0, 2));
        vec4_t yyyy0_1 = VecShuffle(yyyy0, yyyy0, VecShuffleMask(1, 3, 0, 2));
        vec4_t zzzz0_1 = VecShuffle(zzzz0, zzzz0, VecShuffleMask(1, 3, 0, 2));
        vec4_t wwww0_1 = VecShuffle(wwww0, wwww0, VecShuffleMask(1, 3, 0, 2));

        vec4_t xxxx1_1 = VecShuffle(xxxx1, xxxx1, VecShuffleMask(1, 3, 0, 2));
        vec4_t yyyy1_1 = VecShuffle(yyyy1, yyyy1, VecShuffleMask(1, 3, 0, 2));
        vec4_t zzzz1_1 = VecShuffle(zzzz1, zzzz1, VecShuffleMask(1, 3, 0, 2));
        vec4_t wwww1_1 = VecShuffle(wwww1, wwww1, VecShuffleMask(1, 3, 0, 2));

        vec4_t t = VecMul(zzzz1, VecRcp(VecSub( zzzz1, zzzz0)));
        vec4_t new_xxxx0 = INTERSECT_EDGE(xxxx0, xxxx1, t);
        vec4_t new_yyyy0 = INTERSECT_EDGE(yyyy0, yyyy1, t);
        vec4_t new_wwww0 = INTERSECT_EDGE(wwww0, wwww1, t);

        t = VecMul(zzzz0_1, VecRcp(VecSub( zzzz0_1, zzzz0)));
        vec4_t new_xxxx1 = INTERSECT_EDGE(xxxx0, xxxx0_1, t );
        vec4_t new_yyyy1 = INTERSECT_EDGE(yyyy0, yyyy0_1, t );
        vec4_t new_wwww1 = INTERSECT_EDGE(wwww0, wwww0_1, t );

        t = VecMul(zzzz1_1, VecRcp(VecSub(zzzz1_1, zzzz0)));
        vec4_t new_xxxx2 = INTERSECT_EDGE(xxxx1, xxxx1_1, t);
        vec4_t new_yyyy2 = INTERSECT_EDGE(yyyy1, yyyy1_1, t);
        vec4_t new_wwww2 = INTERSECT_EDGE(wwww1, wwww1_1, t);

        vec4_t x0 = VecMul(xxxx0, VecRcp(wwww0));
        vec4_t y0 = VecMul(yyyy0, VecRcp(wwww0));
        vec4_t x1 = VecMul(xxxx1, VecRcp(wwww1));
        vec4_t y1 = VecMul(yyyy1, VecRcp(wwww1));
        vec4_t x2 = VecMul(new_xxxx0, VecRcp(new_wwww0));
        vec4_t y2 = VecMul(new_yyyy0, VecRcp(new_wwww0));
        vec4_t x3 = VecMul(new_xxxx1, VecRcp(new_wwww1));
        vec4_t y3 = VecMul(new_yyyy1, VecRcp(new_wwww1));
        vec4_t x4 = VecMul(new_xxxx2, VecRcp(new_wwww2));
        vec4_t y4 = VecMul(new_yyyy2, VecRcp(new_wwww2));

        x_min = VecMin(VecMin(x0, x1), VecMin(x2, VecMin(x3, x4)));
        x_max = VecMax(VecMax(x0, x1), VecMax(x2, VecMax(x3, x4)));
        y_min = VecMin(VecMin(y0, y1), VecMin(y2, VecMin(y3, y4)));
        y_max = VecMax(VecMax(y0, y1), VecMax(y2, VecMax(y3, y4)));
#undef INTERSECT_EDGE
    }

    vec4_t min_0 = VecMin(VecShuffle(x_min, y_min, VecShuffleMask(0, 1, 0, 1)), VecShuffle(x_min, y_min, VecShuffleMask(2, 3, 2, 3)));
    vec4_t max_0 = VecMax(VecShuffle(x_max, y_max, VecShuffleMask(0, 1, 0, 1)), VecShuffle(x_max, y_max, VecShuffleMask(2, 3, 2, 3)));

    vec4_t min_1 = VecMin(VecShuffle(min_0, min_0, VecShuffleMask(0, 2, 0, 0)), VecShuffle(min_0, min_0, VecShuffleMask(1, 3, 0, 0)));
    vec4_t max_1 = VecMax(VecShuffle(max_0, max_0, VecShuffleMask(0, 2, 0, 0)), VecShuffle(max_0, max_0, VecShuffleMask(1, 3, 0, 0)));

    VecIntStore(bounds_array, VecFloat2Int(get_tile_bounds(VecMoveLH(min_1, max_1))));
    return mask == 15 && intersect_near == false;
}

__forceinline vec4_t Rasterizer::get_tile_bounds(vec4_t minmax)
{
    vec4_t tile_bounds = VecMad(minmax, m_tile_size, m_almost_one);
    return VecMax(VecMin(tile_bounds, m_tile_bounds), VecZero());
}

void Rasterizer::flush_thread_data(ThreadData& RESTRICT thread_data)
{
    assert(m_mt);

    // copy back triangles & indices
    uint32_t triangle_offset = 0;
    if (thread_data.data.triangle_count)
    {
        triangle_offset = atomic_add(m_data.data.triangle_count, thread_data.data.triangle_count);
        assert(triangle_offset + thread_data.data.triangle_count <= m_data.data.triangles.size());
        for (uint32_t i = 0; i < thread_data.data.triangle_count; ++i)
            m_data.data.triangle_data[triangle_offset + i] = thread_data.data.triangle_data[i];
    }
    for (uint32_t tile_index = 0; tile_index < g_width*g_height; ++tile_index)
    {
        auto & tile = thread_data.data.tiles[tile_index];
        if (tile.triangle_index_count == 0)
            continue;
        uint32_t offset = atomic_add(m_data.data.tiles[tile_index].triangle_index_count, tile.triangle_index_count);
        assert(offset + tile.triangle_index_count <= m_data.data.tiles[tile_index].triangle_indices.size());
        for (uint32_t i = 0; i < tile.triangle_index_count; ++i)
        {
            auto & key = tile.triangle_index_data[i];
            key.index += triangle_offset;
            m_data.data.tiles[tile_index].triangle_index_data[offset + i] = key;
        }
    }
    thread_data.clear();
}

__forceinline static void TransformVertices(const vec4_t* RESTRICT m, const vec4_t* RESTRICT src, uint32_t count, vec2_t* RESTRICT dst, uint32_t* RESTRICT depths)
{
    vec4_t w_quantizer = Vector4(64.f);
    vec4_t xy_quantizer = Vector4(float(1<<Rasterizer::g_fixed_point_bits));
    for (uint32_t i = 0; i < count; i += 4, src += 4, dst += 4, depths += 4)
    {
        vec4_t src0 = src[0];
        vec4_t src1 = src[1];
        vec4_t src2 = src[2];
        vec4_t src3 = src[3];

        vec4_t tmp0 = VecUnpackLo(src0, src1); // x0 x1 y0 y1
        vec4_t tmp1 = VecUnpackHi(src0, src1); // z0 z1 w0 w1
        vec4_t tmp2 = VecUnpackLo(src2, src3); // x2 x3 y2 y3
        vec4_t tmp3 = VecUnpackHi(src2, src3); // z2 z3 w2 w3

        vec4_t xxxx = VecMoveLH(tmp0, tmp2); // x0 x1 x2 x3
        vec4_t yyyy = VecMoveHL(tmp0, tmp2); // y0 y1 y2 y3
        vec4_t zzzz = VecMoveLH(tmp1, tmp3); // z0 z1 z2 z3

        vec4_t tw = VecMad(m[11], zzzz, VecMad(m[7], yyyy, VecMad(m[3], xxxx, m[15])));
        vec4_t inv_tw = VecMul(VecRcp(tw), xy_quantizer);

        vec4_t tx = VecMul(inv_tw, VecMad(m[ 8], zzzz, VecMad(m[4], yyyy, VecMad(m[0], xxxx, m[12]))));
        vec4_t ty = VecMul(inv_tw, VecMad(m[ 9], zzzz, VecMad(m[5], yyyy, VecMad(m[1], xxxx, m[13]))));

        VecStore(dst, VecUnpackLo(tx, ty));
        VecStore(dst + 2, VecUnpackHi(tx, ty));
        VecIntStore(depths, VecFloat2Int(VecMul(tw, w_quantizer)));
    }
}

void Rasterizer::push_objects(const Object* RESTRICT objects, uint32_t object_count, uint32_t thread_index)
{
    ThreadData & thread_data = m_mt ? m_thread_data[thread_index] : m_data;
    if (m_mt)
        thread_data.clear();

    uint32_t flags_start = 0, **flag_data = m_flags.data();
    if (m_mt)
    {
        flags_start = atomic_add(m_flag_count, object_count);
        assert(flags_start + object_count <= m_flags.size());
    }
    else
    {
        assert(m_flag_count + object_count <= m_flags.size());
        flags_start = m_flag_count, m_flag_count += object_count;
    }
    assert(flags_start + object_count <= m_flags.size());

    ALIGN16 int bounds_array[4] = { 0 };
    vec4_t matrix[ 16 ];

    uint16_t output_indices[1024];

#if USE_STATS
    uint32_t triangles_total = 0;
    uint32_t triangles_occluder_total = 0;
    uint32_t triangles_occludee_total = 0;
    uint32_t triangles_offscreen = 0;
#endif
    for (uint32_t idx = 0; idx < object_count; ++idx)
    {
        auto & obj = objects[idx];

        ExtractMatrix(obj.transform * m_transform, matrix);

        uint32_t normal_mask = 0;
#if USE_NORMAL_MASKS
        for (uint32_t n = 0; n < obj.normal_count; ++n)
        {
            vec4_t normal = Vector3TransformNormal(obj.transform, obj.normals[n]);
            normal_mask |= (Vector3Dot(normal, m_camera_direction) < 0 ? 1 : 0) << n;
        }
#endif

#if USE_STATS
        triangles_total += obj.index_count / 3;
        if (obj.visibility)
            triangles_occludee_total += obj.index_count / 3;
        else
            triangles_occluder_total += obj.index_count / 3;
#endif
        if (obj.visibility)
            *obj.visibility = 0;

        uint32_t flag = obj.visibility ? flags_start + idx : 0;
        flag_data[flag] = (uint32_t*)obj.visibility;

        bool inside = occlude_object(matrix, obj.bound_min, obj.bound_max, bounds_array);
        if (bounds_array[0] == bounds_array[2] || bounds_array[1] == bounds_array[3])
        {
#if USE_STATS
            triangles_offscreen += obj.index_count / 3;
#endif
            continue;
        }

        uint32_t triangles_4count = inside ? obj.index_count / 12 : obj.index_count / 6;
        if (m_mt && thread_data.data.triangle_count + triangles_4count > thread_data.data.triangles.size())
            flush_thread_data(thread_data);
        else
            assert(thread_data.data.triangle_count + triangles_4count <= thread_data.data.triangles.size());

        vec2_t* positions = thread_data.positions.data();
        uint32_t* depths = thread_data.depths.data();

        bool select_tiles = !(bounds_array[0] + 2 > bounds_array[2] && bounds_array[1] + 2 > bounds_array[3]);
        if (inside)
        {
            TransformVertices(matrix, obj.vertices, obj.vertex_count, positions, depths);
            push_triangle_batched(thread_data.data, flag, positions, depths, obj.indices, obj.index_count, obj.normal_masks, normal_mask, bounds_array, select_tiles);
        }
        else
        {
            vec4_t* transformed_vertices = thread_data.vertices.data();

            size_t aligned_count = (obj.vertex_count + 3) & ~3;
            assert(aligned_count < thread_data.vertices.size());
            for (size_t i = 0; i < aligned_count; i += 4)
                Vector3TransformCoord4Homogeneous(matrix, obj.vertices + i, transformed_vertices + i);

            uint32_t vertex_count = obj.vertex_count;
            uint32_t clipped_indices = clip_triangles(transformed_vertices, vertex_count, obj.indices, obj.index_count, nullptr, 0, output_indices);
            if (clipped_indices == 0)
            {
#if USE_STATS
                triangles_offscreen += obj.index_count / 3;
#endif
                continue;
            }

            vec4_t w_quantizer = Vector4(64.f);
            vec4_t xy_quantizer = Vector4(float(1<<Rasterizer::g_fixed_point_bits));

            aligned_count = vertex_count & ~3;
            for (size_t i = 0; i < aligned_count; i += 4)
            {
                vec4_t src0 = transformed_vertices[i + 0];
                vec4_t src1 = transformed_vertices[i + 1];
                vec4_t src2 = transformed_vertices[i + 2];
                vec4_t src3 = transformed_vertices[i + 3];

                vec4_t tmp0 = VecUnpackLo(src0, src1); // x0 x1 y0 y1
                vec4_t tmp1 = VecUnpackHi(src0, src1); // z0 z1 w0 w1
                vec4_t tmp2 = VecUnpackLo(src2, src3); // x2 x3 y2 y3
                vec4_t tmp3 = VecUnpackHi(src2, src3); // z2 z3 w2 w3

                vec4_t tw = VecMoveHL(tmp1, tmp3);
                vec4_t inv_tw = VecMul(VecRcp(tw), xy_quantizer);

                vec4_t xxxx = VecMul(inv_tw, VecMoveLH(tmp0, tmp2)); // x0 x1 x2 x3
                vec4_t yyyy = VecMul(inv_tw, VecMoveHL(tmp0, tmp2)); // y0 y1 y2 y3

                VecStore(positions + i, VecUnpackLo(xxxx, yyyy));
                VecStore(positions + i + 2, VecUnpackHi(xxxx, yyyy));
                VecIntStore(depths + i, VecFloat2Int(VecMul(tw, w_quantizer)));
            }

            push_triangle_batched(thread_data.data, flag, positions, depths, output_indices, clipped_indices, nullptr, 0, bounds_array, select_tiles);
        }
    }

    if (m_mt)
    {
#if USE_STATS
        atomic_add(m_triangles_total, triangles_total);
        atomic_add(m_triangles_offscreen, triangles_offscreen);
        atomic_add(m_triangles_occludee_total, triangles_occludee_total);
        atomic_add(m_triangles_occluder_total, triangles_occluder_total);
#endif

        flush_thread_data(thread_data);
    }
    else
    {
#if USE_STATS
        m_triangles_total += triangles_total;
        m_triangles_offscreen += triangles_offscreen;
        m_triangles_occludee_total += triangles_occludee_total;
        m_triangles_occluder_total += triangles_occluder_total;
#endif
    }
}

static uint32_t InitThreadData(Rasterizer::ThreadData& data, uint32_t max_4triangles, uint32_t max_triangle_indices)
{
    uint32_t total_mem = 0;

    constexpr uint32_t max_vertices = 1024;
    data.vertices.resize(max_vertices);
    data.positions.resize(max_vertices);
    data.depths.resize(max_vertices);
    total_mem += max_vertices * (sizeof(vec4_t) + sizeof(vec2_t) + sizeof(uint32_t));

    data.data.triangles.resize(max_4triangles);
    data.data.triangle_data = data.data.triangles.data();
    data.data.triangle_count = 0;
    total_mem += max_4triangles * sizeof(Rasterizer::TriangleType);
    for (uint32_t i = 0; i < Rasterizer::g_width*Rasterizer::g_height; ++i)
    {
        auto & tile = data.data.tiles[i];
        tile.triangle_indices.resize(max_triangle_indices);
        tile.triangle_index_data = tile.triangle_indices.data();
        tile.triangle_index_count = 0;
        total_mem += max_triangle_indices * sizeof(Rasterizer::SortKey);
    }

    return total_mem;
}

inline static int Shift( int val )
{
    if( val > 31 )
        return 0;
    if( val < 0 )
        return 0xffffffff;
    return 0xffffffff >> val;
}

void Rasterizer::Init(uint32_t num_threads)
{
    uint32_t total_mem = 0, thread_mem = 0, main_data = 0;

    m_tiles.reserve(g_width*g_height);
    for (int j = 0; j < g_height; ++j)
        for (int i = 0; i < g_width; ++i)
        {
            m_tiles.push_back(Tile(i, j));
        }

    m_tile_height_v = Vector4(Tile::g_tile_height);

    m_masks.resize(g_width*g_max_masks_per_tile);
    for (int tile = 0; tile < g_width; ++tile)
    {
        for (int i = 0; i < tile*128; ++i)
            m_masks[tile*g_max_masks_per_tile + i] = Vector4Int(0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff);
        for (int i = (tile+1)*128; i < g_max_masks_per_tile; ++i)
            m_masks[tile*g_max_masks_per_tile + i] = VecIntZero();
        for (int i = 0; i < 128; ++i)
            m_masks[tile*g_max_masks_per_tile + tile*128 + i] = Vector4Int(Shift(i - 96), Shift(i - 64), Shift(i - 32), Shift(i));
    }

    m_full_span = Vector4Int(0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff);

    constexpr uint32_t max_4triangles = 260 * 1024, max4_triangles_thread = 1024;
    constexpr uint32_t max_triangle_indices = 16*1024, max_triangle_indices_thread = max4_triangles_thread;
    main_data = InitThreadData(m_data, max_4triangles, max_triangle_indices);

    assert(num_threads >= 1);
    m_thread_data.resize(num_threads);
    for (auto & th : m_thread_data)
    {
        thread_mem = InitThreadData(th, max4_triangles_thread, max_triangle_indices_thread) + sizeof(SortKey)*max_triangle_indices;
        th.sort.resize(max_triangle_indices);

        total_mem += thread_mem;
    }
    total_mem += main_data;

    printf("total mem kb %d (%d %d*%d)\n", total_mem >> 10, main_data >> 10, thread_mem >> 10, num_threads);

    m_tile_size = Vector4(1.f / (float)Tile::g_tile_width, 1.f / (float)Tile::g_tile_height, 1.f / (float)Tile::g_tile_width, 1.f / (float)Tile::g_tile_height);
    m_almost_one = Vector4(0.0f, 0.0f, 0.999f, 0.999f);
    m_tile_bounds = Vector4((float)Rasterizer::g_width, (float)Rasterizer::g_height, (float)Rasterizer::g_width, (float)Rasterizer::g_height);
    m_total_size = m_tile_bounds * Vector4((float)(Tile::g_tile_width<<g_fixed_point_bits),(float)(Tile::g_tile_height<<g_fixed_point_bits), 0.f, 0.f);

    constexpr uint32_t max_flags = 128*1024 + 1;
    m_flags.resize(max_flags);

    m_fixed_point = Vector4(1<<g_fixed_point_bits);
    m_inv_fixed_point = VecRcp(m_fixed_point);
}

void Rasterizer::begin(const Matrix& m, vec4_t cam)
{
    m_transform = m;
    m_camera_direction = cam;

    m_data.clear();
    m_flag_count = 1;

    for (auto & t : m_tiles)
        t.clear();

    m_triangles_total = 0;
    m_triangles_occluder_total = 0;
    m_triangles_occludee_total = 0;
    m_triangles_drawn_total = 0;
    m_triangles_drawn_occluder_total = 0;
    m_triangles_drawn_occludee_total = 0;
    m_triangles_skipped = 0;
    m_triangles_offscreen = 0;
    m_triangles_backface = 0;
}

void Rasterizer::ThreadData::clear()
{
    for (uint32_t i = 0; i < g_width*g_height; ++i)
        data.tiles[i].triangle_index_count = 0;
    data.triangle_count = 0;
}
