
#include "rasterizer.h"
#include <stdio.h>

#include <algorithm>

inline static void sort(vec4_t &A0, vec4_t &A1, vec4_t &B0, vec4_t &B1 )
{
    vec4_t mask = VecCmpLe(B0, B1);
    vec4_t sx = VecAdd(A0, A1);
    A0 = VecOr(VecAnd(mask, A0), VecAndNot(mask, A1));
    A1 = VecSub(sx, A0);

    vec4_t sy = VecAdd(B0, B1);
    B0 = VecOr(VecAnd(mask, B0), VecAndNot(mask, B1));
    B1 = VecSub(sy, B0);
}

inline static void ExtractMatrix(const Matrix& m, vec4_t* matrix)
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
__forceinline static int clip_triangle(const vec4_t* input, int count, vec4_t* output, vec4_t plane)
{
    int vertices = 0;
    for ( int i = 0; i < count; ++i )
    {
        vec4_t a0 = input[i*3 + 0];
        vec4_t b0 = input[i*3 + 1];
        vec4_t c0 = input[i*3 + 2];
        vec4_t tmp = VecShuffle(a0, b0, VecShuffleMask( vertex_component, 0, vertex_component, 0));
        vec4_t w = VecShuffle(tmp, c0, VecShuffleMask( 0, 2, vertex_component, vertex_component));

        int mask = cmp_func ? VecMask(VecCmpGt(w, VecMul(VecShuffle(VecShuffle(a0, b0, VecShuffleMask(3, 3, 3, 3)), c0, VecShuffleMask(0, 2, 3, 3)), plane))) & 7 : VecMask(VecCmpLt(w, VecZero())) & 7;
        switch ( mask )
        {
        case 7:
            break;
        case 0:
            output[vertices++] = a0;
            output[vertices++] = b0;
            output[vertices++] = c0;
            break;
        case 1:
        case 2:
        case 4:
            {
                if ( mask & 1 ) { vec4_t t = a0; a0 = b0; b0 = c0; c0 = t; }
                else if ( mask & 2 ) { vec4_t t = a0; a0 = c0; c0 = b0; b0 = t; }

                vec4_t ba = intersectLineZ<vertex_component, cmp_func>( c0, a0, plane );
                vec4_t bb = intersectLineZ<vertex_component, cmp_func>( c0, b0, plane );
                output[vertices++] = a0;
                output[vertices++] = b0;
                output[vertices++] = bb;
                output[vertices++] = bb;
                output[vertices++] = ba;
                output[vertices++] = a0;
            }
            break;
        default:
            {
                vec4_t in = ( mask & 1 ) == 0 ? a0 : ( ( mask & 4 ) == 0 ? c0 : b0 );
                output[vertices++] = mask & 1 ? intersectLineZ< vertex_component, cmp_func >( a0, in, plane ) : a0;
                output[vertices++] = mask & 2 ? intersectLineZ< vertex_component, cmp_func >( b0, in, plane ) : b0;
                output[vertices++] = mask & 4 ? intersectLineZ< vertex_component, cmp_func >( c0, in, plane ) : c0;
            }
            break;
        }
    }
    assert( vertices < 196 );
    return vertices / 3;
}

__forceinline int static clip_triangle(vec4_t* v, vec4_t* dst)
{
    vec4_t g_total_width_v = Vector4(Rasterizer::g_total_width);
    vec4_t g_total_height_v = Vector4(Rasterizer::g_total_height);

    int count = 4;
    vec4_t input_array[196], output_array[196];
    count = clip_triangle< 1, false >(v, count, input_array, VecZero()); // y < 0
    count = clip_triangle< 0, false >(input_array, count, output_array, VecZero()); // x < 0
    count = clip_triangle< 0, true >(output_array, count, input_array, g_total_width_v); // x > 1280
    count = clip_triangle< 2, false >(input_array, count, output_array, VecZero()); // z < 0
    count = clip_triangle< 1, true >(output_array, count, dst, g_total_height_v); // y > 720
    return count;
}

template <bool select_tiles>
__forceinline void Rasterizer::push_4triangles(TrianagleData& data, uint32_t flag, int* bounds_array,
                                               const vec4_t* x, const vec4_t* y, const vec4_t* w)
{
    const vec4_t local_fixed_point = Vector4(1<<g_fixed_point_bits);
    vec4_t x0 = VecMul(x[0], local_fixed_point);
    vec4_t x1 = VecMul(x[1], local_fixed_point);
    vec4_t x2 = VecMul(x[2], local_fixed_point);
    vec4_t y0 = y[0];
    vec4_t y1 = y[1];
    vec4_t y2 = y[2];

    uint32_t mask = VecMask(VecCmpLt(VecSub(VecMul(VecSub(x1, x0), VecSub(y2, y0)), VecMul(VecSub(x2, x0), VecSub(y1, y0))), VecZero()));
    if (mask == 0)
        return;

    sort(x0, x1, y0, y1);
    sort(x1, x2, y1, y2);
    sort(x0, x1, y0, y1);

    assert(data.triangle_count < data.triangles.size());
    TriangleType & t = data.triangle_data[data.triangle_count];

    t.mask = mask;
    t.flag = flag;
#if USE_PACKED_TRIANGLES
    vec4i_t x0x1 = VecIntPack16(VecFloat2Int(x0), VecFloat2Int(x1));
    vec4i_t x2y0 = VecIntPack16(VecFloat2Int(x2), VecFloat2Int(VecMul(y0, local_fixed_point)));
    vec4i_t y1y2 = VecIntPack16(VecFloat2Int(VecMul(y1, local_fixed_point)), VecFloat2Int(VecMul(y2, local_fixed_point)));


    VecIntStore(t.x0, x0x1);
    VecIntStore(t.x2, x2y0);
    VecIntStore(t.y1, y1y2);
#else
    t.x0 = x0;
    t.x1 = x1;
    t.x2 = x2;
    t.y0 = y0;
    t.y1 = y1;
    t.y2 = y2;
#endif

    uint32_t index = data.triangle_count++;

    constexpr float z_quantizator = 64.f;
    ALIGN16 int zz[4];
    if (flag)
    {
        vec4_t ww = VecMul(VecMin(w[0], VecMin(w[1], w[2])), Vector4(z_quantizator));
        vec4_t www = VecMin(ww, VecShuffle(ww, ww, VecShuffleMask(2, 3, 2, 3)));
        VecIntStore(zz, VecFloat2Int(VecMin(www, VecShuffle(www, www, VecShuffleMask(1, 1, 1, 1)))));
    }
    else
    {
        vec4_t ww = VecMul(VecMax(w[0], VecMax(w[1], w[2])), Vector4(z_quantizator));
        vec4_t www = VecMax(ww, VecShuffle(ww, ww, VecShuffleMask(2, 3, 2, 3)));
        VecIntStore(zz, VecFloat2Int(VecMax(www, VecShuffle(www, www, VecShuffleMask(1, 1, 1, 1)))));
    }
    assert(zz[0] >= 0 && zz[0] <= 65535*z_quantizator);
    uint32_t z = zz[0];

    ALIGN16 int transformed_bounds[4];
    if ( select_tiles )
    {
        vec4_t x_min = VecMin(x[0], VecMin(x[1], x[2]));
        vec4_t x_max = VecMax(x[0], VecMax(x[1], x[2]));
        vec4_t y_min = VecMin(y[0], VecMin(y[1], y[2]));
        vec4_t y_max = VecMax(y[0], VecMax(y[1], y[2]));

        vec4_t min_0 = VecMin(VecShuffle(x_min, y_min, VecShuffleMask(0, 1, 0, 1)), VecShuffle(x_min, y_min, VecShuffleMask(2, 3, 2, 3)));
        vec4_t max_0 = VecMax(VecShuffle(x_max, y_max, VecShuffleMask(0, 1, 0, 1)), VecShuffle(x_max, y_max, VecShuffleMask(2, 3, 2, 3)));

        vec4_t min_1 = VecMin(VecShuffle(min_0, min_0, VecShuffleMask(0, 2, 0, 0)), VecShuffle(min_0, min_0, VecShuffleMask(1, 3, 0, 0)));
        vec4_t max_1 = VecMax(VecShuffle(max_0, max_0, VecShuffleMask(0, 2, 0, 0)), VecShuffle(max_0, max_0, VecShuffleMask(1, 3, 0, 0)));

        vec4_t tile_bound = get_tile_bounds(min_1, max_1);
        VecIntStore(transformed_bounds, VecFloat2Int( tile_bound ));
        assert(transformed_bounds[0] >= bounds_array[0]);
        assert(transformed_bounds[1] >= bounds_array[1]);
        assert(transformed_bounds[2] <= bounds_array[2]);
        assert(transformed_bounds[3] <= bounds_array[3]);

        bounds_array = transformed_bounds;
    }
    for (int yy = bounds_array[1]; yy < bounds_array[3]; ++yy)
        for (int xx = bounds_array[0]; xx < bounds_array[2]; ++xx)
        {
            uint32_t tile_index = xx + yy*g_width;
            assert(xx < g_width);
            assert(yy < g_height);
            auto & tile = data.tiles[tile_index];
            assert(tile.triangle_index_count < tile.triangle_indices.size());
            tile.triangle_index_data[tile.triangle_index_count++] = {z, index};
        }
}

template <bool use_indices>
__forceinline void load_4vertices(vec4_t& x, vec4_t& y, vec4_t& w, const vec4_t* src, const unsigned short* indices, uint32_t base_index)
{
#define IDX(num)(use_indices ? indices[base_index + num] : base_index + num)
    vec4_t v0_0 = src[IDX(0)];
    vec4_t v0_1 = src[IDX(3)];
    vec4_t v0_2 = src[IDX(6)];
    vec4_t v0_3 = src[IDX(9)];
    vec4_t tmp0 = VecUnpackLo(v0_0, v0_1); // x0_0 x1_0 y0_0 y1_0
    vec4_t tmp1 = VecUnpackLo(v0_2, v0_3); // x2_0 x3_0 y2_0 y3_0
    vec4_t tmp0_0 = VecUnpackHi(v0_0, v0_1); // z0_0 z1_0 w0_0 w1_0
    vec4_t tmp1_0 = VecUnpackHi(v0_2, v0_3); // z2_0 z3_0 w2_0 w3_0
    w = VecMoveHL(tmp0_0, tmp1_0); // w0_0 w1_0 w2_0 w3_0
    x = VecMul(VecMoveLH(tmp0, tmp1), VecRcp(w));
    y = VecMul(VecMoveHL(tmp0, tmp1), VecRcp(w));
#undef IDX
}

template <bool select_tiles, bool use_indices>
__forceinline void Rasterizer::push_triangle_batched(TrianagleData& data,uint32_t flag, const vec4_t* src, int count, const unsigned short* indices, int* bounds_array)
{
    assert(( (count / 3) & 3 ) == 0);
    for ( int i = 0; i < count; i += 12 )
    {
        vec4_t x[3], y[3], w[3];
        load_4vertices<use_indices>(x[0], y[0], w[0], src, indices, i + 0);
        load_4vertices<use_indices>(x[1], y[1], w[1], src, indices, i + 1);
        load_4vertices<use_indices>(x[2], y[2], w[2], src, indices, i + 2);
        push_4triangles<select_tiles>(data, flag, bounds_array, x, y, w);
    }
}

template < bool select_tiles >
void Rasterizer::push_object_clipped(ThreadData& thread_data, const vec4_t* matrix, const uint16_t* indices, int index_count,
                                     const vec4_t* vertices, int vertex_count, int* bounds_array, uint32_t flag)
{
    vec4_t* transformed_vertices = thread_data.vertices.data();

    size_t aligned_count = (vertex_count + 3 ) & ~3;
    assert(aligned_count < thread_data.vertices.size());
    for (size_t i = 0; i < aligned_count; i += 4)
        Vector3TransformCoord4Homogeneous(matrix, vertices + i, transformed_vertices + i);

    assert(index_count >= 12);

    constexpr uint32_t max_triangles_in_object = 1024;
    vec4_t clipped_triangles[max_triangles_in_object*3];

    uint32_t clipped_triangle_count = 0;
    for (int i = 0; i < index_count; i += 12)
    {
        assert(indices[ i + 0 ] < vertex_count);
        assert(indices[ i + 1 ] < vertex_count);
        assert(indices[ i + 2 ] < vertex_count);

        vec4_t v[12];

        v[0] = transformed_vertices[indices[i + 0]];
        v[1] = transformed_vertices[indices[i + 1]];
        v[2] = transformed_vertices[indices[i + 2]];

        v[3] = transformed_vertices[indices[i + 3]];
        v[4] = transformed_vertices[indices[i + 4]];
        v[5] = transformed_vertices[indices[i + 5]];

        v[6] = transformed_vertices[indices[i + 6]];
        v[7] = transformed_vertices[indices[i + 7]];
        v[8] = transformed_vertices[indices[i + 8]];

        v[9] = transformed_vertices[indices[i + 9 ]];
        v[10] = transformed_vertices[indices[i + 10 ]];
        v[11] = transformed_vertices[indices[i + 11 ]];

        clipped_triangle_count += clip_triangle(v, clipped_triangles + clipped_triangle_count * 3);
    }
    assert(clipped_triangle_count < max_triangles_in_object);

    int tris_to_pad = ( ( clipped_triangle_count + 3 ) & ~3 ) - clipped_triangle_count;
    vec4_t& v0 = clipped_triangles[ clipped_triangle_count * 3 - 3 ];
    vec4_t& v1 = clipped_triangles[ clipped_triangle_count * 3 - 2 ];
    vec4_t& v2 = clipped_triangles[ clipped_triangle_count * 3 - 1 ];
    for ( int i = 0; i < tris_to_pad; ++i, ++clipped_triangle_count )
    {
        clipped_triangles[ clipped_triangle_count*3 + 0 ] = v0;
        clipped_triangles[ clipped_triangle_count*3 + 1 ] = v1;
        clipped_triangles[ clipped_triangle_count*3 + 2 ] = v2;
    }
    push_triangle_batched<select_tiles, false>(thread_data.data, flag, clipped_triangles, clipped_triangle_count*3, 0, bounds_array);
}

bool Rasterizer::occlude_object(const vec4_t* m, vec4_t v_min, vec4_t v_max, int* bounds_array)
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

    VecIntStore(bounds_array, VecFloat2Int(get_tile_bounds(min_1, max_1)));
    return mask == 15 && intersect_near == false;
}

__forceinline vec4_t Rasterizer::get_tile_bounds( vec4_t v_min, vec4_t v_max )
{
    vec4_t minmax = VecMoveLH(v_min, v_max); // xyXY
    vec4_t tile_bounds = VecMad(minmax, m_tile_size, m_almost_one);
    return VecMax(VecMin(tile_bounds, m_tile_bounds), VecZero());
}

void Rasterizer::flush_thread_data(ThreadData& thread_data)
{
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
            m_data.data.tiles[tile_index].triangle_index_data[offset + i] = {tile.triangle_index_data[i].z, tile.triangle_index_data[i].index + triangle_offset};
    }
    thread_data.clear();
}

void Rasterizer::push_objects(const Object* objects, uint32_t object_count, uint32_t thread_index)
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

    uint32_t triangles_total = 0;
    uint32_t triangles_occluder_total = 0;
    uint32_t triangles_occludee_total = 0;
    uint32_t triangles_offscreen = 0;
    for (uint32_t idx = 0; idx < object_count; ++idx)
    {
        auto & obj = objects[idx];

        ExtractMatrix(obj.transform * m_transform, matrix);

        triangles_total += obj.index_count / 3;
        if (obj.visibility)
            *obj.visibility = 0, triangles_occludee_total += obj.index_count / 3;
        else
            triangles_occluder_total += obj.index_count / 3;

        uint32_t flag = obj.visibility ? flags_start + idx : 0;
        flag_data[flag] = (uint32_t*)obj.visibility;

        bool inside = occlude_object(matrix, obj.bound_min, obj.bound_max, bounds_array);
        if ( bounds_array[0] == bounds_array[2] || bounds_array[1] == bounds_array[3] )
        {
            triangles_offscreen += obj.index_count / 3;
            continue;
        }

        uint32_t triangles_4count = inside ? obj.index_count / 12 : obj.index_count / 6;
        if (m_mt && thread_data.data.triangle_count + triangles_4count > thread_data.data.triangles.size())
            flush_thread_data(thread_data);
        else
            assert(thread_data.data.triangle_count + triangles_4count <= thread_data.data.triangles.size());

        if (inside)
        {
            vec4_t* transformed_vertices = thread_data.vertices.data();

            size_t aligned_count = (obj.vertex_count + 3) & ~3;
            assert(aligned_count < thread_data.vertices.size());
            for (size_t i = 0; i < aligned_count; i += 4)
                Vector3TransformCoord4Homogeneous( matrix, obj.vertices + i, transformed_vertices + i );

            if ( bounds_array[0] + 2 > bounds_array[2] && bounds_array[1] + 2 > bounds_array[3] )
                push_triangle_batched<false, true>(thread_data.data, flag, transformed_vertices, obj.index_count, obj.indices, bounds_array );
            else
                push_triangle_batched<true, true>(thread_data.data, flag, transformed_vertices, obj.index_count, obj.indices, bounds_array );
        }
        else
        {
            if ( bounds_array[0] + 2 > bounds_array[2] && bounds_array[1] + 2 > bounds_array[3] )
                push_object_clipped<false>(thread_data, matrix, obj.indices, obj.index_count, obj.vertices, obj.vertex_count, bounds_array, flag);
            else
                push_object_clipped<true>(thread_data, matrix, obj.indices, obj.index_count, obj.vertices, obj.vertex_count, bounds_array, flag);
        }
    }

    if (m_mt)
    {
        atomic_add(m_triangles_total, triangles_total);
        atomic_add(m_triangles_offscreen, triangles_offscreen);
        atomic_add(m_triangles_occludee_total, triangles_occludee_total);
        atomic_add(m_triangles_occluder_total, triangles_occluder_total);

        flush_thread_data(thread_data);
    }
    else
    {
        m_triangles_total += triangles_total;
        m_triangles_offscreen += triangles_offscreen;
        m_triangles_occludee_total += triangles_occludee_total;
        m_triangles_occluder_total += triangles_occluder_total;
    }
}

static uint32_t InitThreadData(Rasterizer::ThreadData& data, uint32_t max_4triangles, uint32_t max_triangle_indices)
{
    uint32_t total_mem = 0;

    constexpr uint32_t max_vertices = 1024;
    data.vertices.resize(max_vertices);
    total_mem += max_vertices * sizeof(vec4_t);


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
            total_mem += sizeof(Tile) + m_tiles.back().m_shifts.size()*sizeof(vec4i_t);
        }

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

    constexpr uint32_t max_4triangles = 384 * 1024, max4_triangles_thread = 1024;
    constexpr uint32_t max_triangle_indices = 48*1024, max_triangle_indices_thread = 1024;
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
    m_almost_one = Vector4(0.f, 0.f, 0.9999f, 0.9999f);
    m_tile_bounds = Vector4((float)Rasterizer::g_width, (float)Rasterizer::g_height, (float)Rasterizer::g_width, (float)Rasterizer::g_height);

    constexpr uint32_t max_flags = 128*1024 + 1;
    m_flags.resize(max_flags);
}

void Rasterizer::begin(const Matrix& m)
{
    m_transform = m;

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
}

void Rasterizer::ThreadData::clear()
{
    for (uint32_t i = 0; i < g_width*g_height; ++i)
        data.tiles[i].triangle_index_count = 0;
    data.triangle_count = 0;
}
