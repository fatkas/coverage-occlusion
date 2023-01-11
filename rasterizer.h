#pragma once

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "rasterizer_math.h"

#include "rasterizer_tile.h"

#include <bx/uint32_t.h>

struct _MM_ALIGN16 Rasterizer
{
	static constexpr int g_width = 6;
    static constexpr int g_height = 12;
	static constexpr int g_total_width = g_width * Tile::g_tile_width;
	static constexpr int g_total_height = g_height * Tile::g_tile_height;

    static constexpr int g_max_4triangles = 512 * 1024;

private:

    stl::vector<Tile>       m_tiles;
	Matrix			        m_transform;
    stl::vector<Triangle>   m_triangles;
    uint32_t                m_triangle_count = 0;
    stl::vector<uint64_t>   m_sort;
    __m128i                 m_full_span;
    bool                    m_mt = false;

    vec4_t                  g_tile_size = Vector4(1.f / (float)Tile::g_tile_width, 1.f / (float)Tile::g_tile_height, 1.f / (float)Tile::g_tile_width, 1.f / (float)Tile::g_tile_height);
    vec4_t                  g_almost_one = Vector4(0.f, 0.f, 0.9999f, 0.9999f);
    vec4_t                  g_tile_bounds = Vector4((float)Rasterizer::g_width, (float)Rasterizer::g_height, (float)Rasterizer::g_width, (float)Rasterizer::g_height);

	inline void ExtractMatrix( const Matrix& m, vec4_t* matrix )
	{
		#define EXTRACT(line) matrix[line*4+0] = VecShuffle( m.r[line], m.r[line], VecShuffleMask( 0, 0, 0, 0 ) );  \
							matrix[line*4+1] = VecShuffle( m.r[line], m.r[line], VecShuffleMask( 1, 1, 1, 1 ) ); \
							matrix[line*4+2] = VecShuffle( m.r[line], m.r[line], VecShuffleMask( 2, 2, 2, 2 ) ); \
							matrix[line*4+3] = VecShuffle( m.r[line], m.r[line], VecShuffleMask( 3, 3, 3, 3 ) );
		EXTRACT( 0 );
		EXTRACT( 1 );
		EXTRACT( 2 );
		EXTRACT( 3 );

		#undef EXTRACT
	}

    inline void sort( __m128 &A0, __m128 &A1, __m128 &B0, __m128 &B1 )
    {
        __m128 mask = _mm_cmple_ps( B0, B1 );
        __m128 sx = _mm_add_ps( A0, A1 );
        A0 = _mm_or_ps( _mm_and_ps( mask, A0 ), _mm_andnot_ps( mask, A1 ) );
        A1 = _mm_sub_ps( sx, A0 );

        __m128 sy = _mm_add_ps( B0, B1 );
        B0 = _mm_or_ps( _mm_and_ps( mask, B0 ), _mm_andnot_ps( mask, B1 ) );
        B1 = _mm_sub_ps(sy, B0);
    }

    template <bool select_tiles> inline void push_4triangles(int* flag, int* bounds_array,
                                                             vec4_t x0, vec4_t y0, vec4_t w0,
                                                             vec4_t x1, vec4_t y1, vec4_t w1,
                                                             vec4_t x2, vec4_t y2, vec4_t w2)
    {
        Triangle t;

        t.flag = flag;

        const vec4_t local_fixed_point = Vector4(65536.f);
        t.x0 = VecMul(x0, local_fixed_point);
        t.y0 = y0;
        t.x1 = VecMul(x1, local_fixed_point);
        t.y1 = y1;
        t.x2 = VecMul(x2, local_fixed_point);
        t.y2 = y2;

        t.mask = _mm_movemask_ps(_mm_cmplt_ps(VecSub(VecMul(VecSub(t.x1, t.x0), VecSub(t.y2, t.y0)), VecMul(VecSub(t.x2, t.x0), VecSub(t.y1, t.y0))), VecZero()));
        if (t.mask == 0)
            return;

        assert(m_triangle_count < g_max_4triangles);

        _MM_ALIGN16 int zz[4];
        if (flag)
        {
            vec4_t w = VecMul(VecMin(w0, VecMin(w1, w2)), Vector4(256.0f));
            vec4_t ww = VecMin(VecShuffle(w, w, VecShuffleMask(0, 1, 0, 1)), VecShuffle(w, w, VecShuffleMask(2, 3, 2, 3)));
            VecIntStore(zz, VecFloat2Int(VecMin(VecShuffle(ww, ww, VecShuffleMask(0, 0, 0, 0)), VecShuffle(w, w, VecShuffleMask(1, 1, 1, 1)))));
        }
        else
        {
            vec4_t w = VecMul(VecMax(w0, VecMax(w1, w2)), Vector4(256.0f));
            vec4_t ww = VecMax(VecShuffle(w, w, VecShuffleMask(0, 1, 0, 1)), VecShuffle(w, w, VecShuffleMask(2, 3, 2, 3)));
            VecIntStore(zz, VecFloat2Int(VecMax(VecShuffle(ww, ww, VecShuffleMask(0, 0, 0, 0)), VecShuffle(w, w, VecShuffleMask(1, 1, 1, 1)))));
        }
        assert(zz[0] >= 0 && zz[0] <= 65535*256);
        t.z = zz[0];

        sort(t.x0, t.x1, t.y0, t.y1);
        sort(t.x1, t.x2, t.y1, t.y2);
        sort(t.x0, t.x1, t.y0, t.y1);

        uint64_t index = 0, tri_index = 0;
        if (m_mt)
            index = __atomic_add_fetch(&m_triangle_count, 1, __ATOMIC_SEQ_CST) - 1;
        else
            index = m_triangle_count++;
        m_triangles[index] = t;


        if ( select_tiles )
        {
            vec4_t x_min = VecMin(x0, VecMin(x1, x2));
            vec4_t x_max = VecMax(x0, VecMax(x1, x2));
            vec4_t y_min = VecMin(t.y0, VecMin(t.y1, t.y2));
            vec4_t y_max = VecMax(t.y0, VecMax(t.y1, t.y2));

            vec4_t min_0 = VecMin(VecShuffle(x_min, y_min, VecShuffleMask(0, 1, 0, 1)), VecShuffle(x_min, y_min, VecShuffleMask(2, 3, 2, 3)));
            vec4_t max_0 = VecMax(VecShuffle(x_max, y_max, VecShuffleMask(0, 1, 0, 1)), VecShuffle(x_max, y_max, VecShuffleMask(2, 3, 2, 3)));

            vec4_t min_1 = VecMin(VecShuffle(min_0, min_0, VecShuffleMask(0, 2, 0, 0)), VecShuffle(min_0, min_0, VecShuffleMask(1, 3, 0, 0)));
            vec4_t max_1 = VecMax(VecShuffle(max_0, max_0, VecShuffleMask(0, 2, 0, 0)), VecShuffle(max_0, max_0, VecShuffleMask(1, 3, 0, 0)));

            vec4_t tile_bound = get_tile_bounds(min_1, max_1);

            _MM_ALIGN16 int transformed_bounds[4];
            VecIntStore(transformed_bounds, VecFloat2Int( tile_bound ));
            assert(transformed_bounds[0] >= bounds_array[0]);
            assert(transformed_bounds[1] >= bounds_array[1]);
            assert(transformed_bounds[2] <= bounds_array[2]);
            assert(transformed_bounds[3] <= bounds_array[3]);

            for (int y = transformed_bounds[1]; y < transformed_bounds[3]; ++y)
                for (int x = transformed_bounds[0]; x < transformed_bounds[2]; ++x)
                {
                    assert(x < g_width);
                    assert(y < g_height);
                    auto & tt = m_tiles[x + y*g_width];
                    assert(tt.m_triangle_count < Tile::g_max_triangles);
                    if (m_mt)
                        tri_index = __atomic_add_fetch(&tt.m_triangle_count, 1, __ATOMIC_SEQ_CST) - 1;
                    else
                        tri_index = tt.m_triangle_count++;
                    assert(tri_index < Tile::g_max_triangles);
                    tt.m_triangles[tri_index] = (index<<32)|t.z;
                }
        }
        else
        {
            assert(bounds_array[0] < g_width);
            assert(bounds_array[1] < g_height);
            auto & tt = m_tiles[bounds_array[0] + bounds_array[1]*g_width];
            assert(tt.m_triangle_count < Tile::g_max_triangles);
            if (m_mt)
                tri_index = __atomic_add_fetch(&tt.m_triangle_count, 1, __ATOMIC_SEQ_CST) - 1;
            else
                tri_index = tt.m_triangle_count++;
            assert(tri_index < Tile::g_max_triangles);
            tt.m_triangles[tri_index] = (index<<32)|t.z;
        }
    }

	template < bool select_tiles, bool use_indices > inline void push_triangle_batched(int* flag, const vec4_t* src, int count, const unsigned short* indices, int* bounds_array)
	{
		assert(( (count / 3) & 3 ) == 0);

		for ( int i = 0; i < count; i += 12 )
		{
			#define IDX(num)( use_indices ? indices[ i + num ] : i + num )
            vec4_t v0_0 = src[IDX(0)];
            vec4_t v0_1 = src[IDX(3)];
            vec4_t v0_2 = src[IDX(6)];
            vec4_t v0_3 = src[IDX(9)];
            vec4_t tmp0 = _mm_unpacklo_ps(v0_0, v0_1); // x0_0 x1_0 y0_0 y1_0
            vec4_t tmp1 = _mm_unpacklo_ps(v0_2, v0_3); // x2_0 x3_0 y2_0 y3_0
            vec4_t tmp0_0 = _mm_unpackhi_ps(v0_0, v0_1); // z0_0 z1_0 w0_0 w1_0
            vec4_t tmp1_0 = _mm_unpackhi_ps(v0_2, v0_3); // z2_0 z3_0 w2_0 w3_0
            vec4_t w0 = _mm_movehl_ps(tmp1_0, tmp0_0); // w0_0 w1_0 w2_0 w3_0
            vec4_t x0 = VecMul(_mm_movelh_ps(tmp0, tmp1), VecRcp(w0));
            vec4_t y0 = VecMul(_mm_movehl_ps(tmp1, tmp0), VecRcp(w0));
			
            vec4_t v1_0 = src[IDX(1)];
            vec4_t v1_1 = src[IDX(4)];
            vec4_t v1_2 = src[IDX(7)];
            vec4_t v1_3 = src[IDX(10)];
            vec4_t tmp2 = _mm_unpacklo_ps( v1_0, v1_1 ); // x0_1 x1_1 y0_1 y1_1
            vec4_t tmp3 = _mm_unpacklo_ps( v1_2, v1_3 ); // x2_1 x3_1 y2_1 y3_1
            vec4_t tmp2_0 = _mm_unpackhi_ps(v1_0, v1_1);
            vec4_t tmp3_0 = _mm_unpackhi_ps(v1_2, v1_3);
            vec4_t w1 = _mm_movehl_ps(tmp3_0, tmp2_0); // w0_1 w1_1 w2_1 w3_1
            vec4_t x1 = VecMul(_mm_movelh_ps(tmp2, tmp3), VecRcp(w1));
            vec4_t y1 = VecMul(_mm_movehl_ps(tmp3, tmp2), VecRcp(w1));

            vec4_t v2_0 = src[IDX(2)];
            vec4_t v2_1 = src[IDX(5)];
            vec4_t v2_2 = src[IDX(8)];
            vec4_t v2_3 = src[IDX(11)];
            vec4_t tmp4 = _mm_unpacklo_ps(v2_0, v2_1); // x0_2 x1_2 y0_2 y1_2
            vec4_t tmp5 = _mm_unpacklo_ps(v2_2, v2_3); // x2_2 x3_2 y2_2 y3_2
            vec4_t tmp4_0 = _mm_unpackhi_ps(v2_0, v2_1);
            vec4_t tmp5_0 = _mm_unpackhi_ps(v2_2, v2_3);
            vec4_t w2 = _mm_movehl_ps(tmp5_0, tmp4_0); // w0_2 w1_2 w2_2 w3_2
            vec4_t x2 = VecMul(_mm_movelh_ps(tmp4, tmp5), VecRcp(w2));
            vec4_t y2 = VecMul(_mm_movehl_ps(tmp5, tmp4), VecRcp(w2));
			#undef IDX

            push_4triangles<select_tiles>(flag, bounds_array, x0, y0, w0, x1, y1, w1, x2, y2, w2);
		}
	}

	bool occlude_object(const __m128* m, const __m128& v_min, const __m128& v_max, int* bounds_array)
	{
        vec4_t g_total_width_v = Vector4(g_total_width);
        vec4_t g_total_height_v = Vector4(g_total_height);

        vec4_t pt[ 4 ];
        vec4_t vTmp = _mm_unpacklo_ps(v_min, v_max);				// x, X, y, Y
		pt[0] = VecShuffle(v_min, v_max, _MM_SHUFFLE( 0, 0, 0, 0)); // xxXX
		pt[1] = VecShuffle(vTmp, vTmp, _MM_SHUFFLE( 2, 3, 2, 3)); // yYyY
		pt[2] = VecShuffle(v_min, v_min, _MM_SHUFFLE( 2, 2, 2, 2)); // zzzz
		pt[3] = VecShuffle(v_max, v_max, _MM_SHUFFLE( 2, 2, 2, 2)); // ZZZZ

        vec4_t xxxx0 = VecMad(m[8], pt[2], VecMad(m[4], pt[1], VecMad(m[0], pt[0], m[12])));
        vec4_t yyyy0 = VecMad(m[9], pt[2], VecMad(m[5], pt[1], VecMad(m[1], pt[0], m[13])));
        vec4_t zzzz0 = VecMad(m[10], pt[2], VecMad(m[6], pt[1], VecMad(m[2], pt[0], m[14])));
        vec4_t wwww0 = VecMad(m[11], pt[2], VecMad(m[7], pt[1], VecMad(m[3], pt[0], m[15])));

        vec4_t xxxx1 = VecMad(m[8], pt[3], VecMad(m[4], pt[1], VecMad(m[0], pt[0], m[12])));
        vec4_t yyyy1 = VecMad(m[9], pt[3], VecMad(m[5], pt[1], VecMad(m[1], pt[0], m[13])));
        vec4_t zzzz1 = VecMad(m[10], pt[3], VecMad(m[6], pt[1], VecMad(m[2], pt[0], m[14])));
        vec4_t wwww1 = VecMad(m[11], pt[3], VecMad(m[7], pt[1], VecMad(m[3], pt[0], m[15])));

        vec4_t v_mask00 = _mm_and_ps(_mm_cmpgt_ps(xxxx0, _mm_setzero_ps()), _mm_cmpgt_ps(xxxx1, _mm_setzero_ps()));
        vec4_t v_mask01 = _mm_and_ps(_mm_cmpgt_ps(yyyy0, _mm_setzero_ps()), _mm_cmpgt_ps(yyyy1, _mm_setzero_ps()));
        vec4_t v_mask10 = _mm_and_ps(_mm_cmplt_ps(xxxx0, _mm_mul_ps(wwww0, g_total_width_v)), _mm_cmplt_ps(xxxx1, _mm_mul_ps(wwww1, g_total_width_v )));
        vec4_t v_mask11 = _mm_and_ps(_mm_cmplt_ps(yyyy0, _mm_mul_ps(wwww0, g_total_height_v)), _mm_cmplt_ps(yyyy1, _mm_mul_ps(wwww1, g_total_height_v)));

        vec4_t v_mask0 = _mm_and_ps(v_mask00, v_mask10);
        vec4_t v_mask1 = _mm_and_ps(v_mask01, v_mask11);
		int mask = _mm_movemask_ps(_mm_and_ps(v_mask0, v_mask1));

		bool intersect_near = _mm_movemask_ps(_mm_and_ps(_mm_cmpgt_ps(zzzz0, _mm_setzero_ps()), _mm_cmpgt_ps(zzzz1, _mm_setzero_ps()))) != 15;

        vec4_t x_min, x_max, y_min, y_max;
		if (intersect_near == false)
		{
            vec4_t x0 = VecMul(xxxx0, VecRcp(wwww0));
            vec4_t y0 = VecMul(yyyy0, VecRcp(wwww0));
            vec4_t x1 = VecMul(xxxx1, VecRcp(wwww1));
            vec4_t y1 = VecMul(yyyy1, VecRcp(wwww1));

			x_min = _mm_min_ps(x0, x1);
			x_max = _mm_max_ps(x0, x1);
			y_min = _mm_min_ps(y0, y1);
			y_max = _mm_max_ps(y0, y1);
        }
		else
		{
#define INTERSECT_EDGE(a,b,c) _mm_add_ps(b, _mm_mul_ps(_mm_sub_ps(a, b), t))
            vec4_t xxxx0_1 = VecShuffle(xxxx0, xxxx0, VecShuffleMask(1, 3, 0, 2));
            vec4_t yyyy0_1 = VecShuffle(yyyy0, yyyy0, VecShuffleMask(1, 3, 0, 2));
            vec4_t zzzz0_1 = VecShuffle(zzzz0, zzzz0, VecShuffleMask(1, 3, 0, 2));
            vec4_t wwww0_1 = VecShuffle(wwww0, wwww0, VecShuffleMask(1, 3, 0, 2));

            vec4_t xxxx1_1 = VecShuffle(xxxx1, xxxx1, VecShuffleMask(1, 3, 0, 2));
            vec4_t yyyy1_1 = VecShuffle(yyyy1, yyyy1, VecShuffleMask(1, 3, 0, 2));
            vec4_t zzzz1_1 = VecShuffle(zzzz1, zzzz1, VecShuffleMask(1, 3, 0, 2));
            vec4_t wwww1_1 = VecShuffle(wwww1, wwww1, VecShuffleMask(1, 3, 0, 2));

            vec4_t t = VecMul(zzzz1, VecRcp(_mm_sub_ps( zzzz1, zzzz0)));
            vec4_t new_xxxx0 = INTERSECT_EDGE(xxxx0, xxxx1, t);
            vec4_t new_yyyy0 = INTERSECT_EDGE(yyyy0, yyyy1, t);
            vec4_t new_wwww0 = INTERSECT_EDGE(wwww0, wwww1, t);

			t = VecMul(zzzz0_1, VecRcp(_mm_sub_ps( zzzz0_1, zzzz0)));
            vec4_t new_xxxx1 = INTERSECT_EDGE(xxxx0, xxxx0_1, t );
            vec4_t new_yyyy1 = INTERSECT_EDGE(yyyy0, yyyy0_1, t );
            vec4_t new_wwww1 = INTERSECT_EDGE(wwww0, wwww0_1, t );

			t = VecMul(zzzz1_1, VecRcp(_mm_sub_ps(zzzz1_1, zzzz0)));
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

			x_min = _mm_min_ps(_mm_min_ps(x0, x1), _mm_min_ps(x2, _mm_min_ps(x3, x4)));
			x_max = _mm_max_ps(_mm_max_ps(x0, x1), _mm_max_ps(x2, _mm_max_ps(x3, x4)));
			y_min = _mm_min_ps(_mm_min_ps(y0, y1), _mm_min_ps(y2, _mm_min_ps(y3, y4)));
			y_max = _mm_max_ps(_mm_max_ps(y0, y1), _mm_max_ps(y2, _mm_max_ps(y3, y4)));
#undef INTERSECT_EDGE
		}

        vec4_t min_0 = _mm_min_ps(VecShuffle(x_min, y_min, VecShuffleMask(0, 1, 0, 1)), VecShuffle(x_min, y_min, VecShuffleMask(2, 3, 2, 3)));
        vec4_t max_0 = _mm_max_ps(VecShuffle(x_max, y_max, VecShuffleMask(0, 1, 0, 1)), VecShuffle(x_max, y_max, VecShuffleMask(2, 3, 2, 3)));

        vec4_t min_1 = _mm_min_ps(VecShuffle(min_0, min_0, VecShuffleMask(0, 2, 0, 0)), VecShuffle(min_0, min_0, VecShuffleMask(1, 3, 0, 0)));
        vec4_t max_1 = _mm_max_ps(VecShuffle(max_0, max_0, VecShuffleMask(0, 2, 0, 0)), VecShuffle(max_0, max_0, VecShuffleMask(1, 3, 0, 0)));

		VecIntStore(bounds_array, VecFloat2Int(get_tile_bounds(min_1, max_1)));
		return mask == 15 && intersect_near == false;
	}

	template < bool select_tiles > void push_object_clipped( const vec4_t* matrix, const uint16_t* indices, int index_count, vec4_t* vertices, int vertex_count, int* bounds_array, int* flag )
	{
        constexpr uint32_t max_triangles_in_object = 1024;
		_MM_ALIGN16 vec4_t transformed_vertices[ max_triangles_in_object ];
		_MM_ALIGN16 vec4_t clipped_triangles[ max_triangles_in_object*3 ];

		int aligned_count = ( vertex_count + 3 ) & ~3;
		assert( aligned_count < max_triangles_in_object );
		for ( int i = 0; i < aligned_count; i += 4 )
			Vector3TransformCoord4Homogeneous( matrix, vertices + i, transformed_vertices + i );

        assert(index_count >= 12);

		int clipped_triangle_count = 0;
		for ( int i = 0; i < index_count; i += 12 )
		{
			assert( indices[ i + 0 ] < vertex_count );
			assert( indices[ i + 1 ] < vertex_count );
			assert( indices[ i + 2 ] < vertex_count );

            vec4_t v[12];

			v[0] = transformed_vertices[ indices[ i + 0 ] ];
			v[1] = transformed_vertices[ indices[ i + 1 ] ];
			v[2] = transformed_vertices[ indices[ i + 2 ] ];

			v[3] = transformed_vertices[ indices[ i + 3 ] ];
			v[4] = transformed_vertices[ indices[ i + 4 ] ];
			v[5] = transformed_vertices[ indices[ i + 5 ] ];

			v[6] = transformed_vertices[ indices[ i + 6 ] ];
			v[7] = transformed_vertices[ indices[ i + 7 ] ];
			v[8] = transformed_vertices[ indices[ i + 8 ] ];

			v[9] = transformed_vertices[ indices[ i + 9 ] ];
			v[10] = transformed_vertices[ indices[ i + 10 ] ];
			v[11] = transformed_vertices[ indices[ i + 11 ] ];

			clipped_triangle_count += clip_triangle( v, clipped_triangles + clipped_triangle_count * 3 );
		}
		assert( clipped_triangle_count < max_triangles_in_object );

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
		push_triangle_batched< select_tiles, false >( flag, clipped_triangles, clipped_triangle_count*3, 0, bounds_array );
	}

	template< int v, bool use_plane > __forceinline vec4_t intersectLineZ( vec4_t a, vec4_t b, vec4_t plane )
	{		
		//  t = (a.x - a.w) / (b.w - a.w - b.x + a.x);
        vec4_t bz = VecShuffle(b, b, VecShuffleMask(v, v, v, v));
        vec4_t az = VecShuffle(a, a, VecShuffleMask(v, v, v, v));
        vec4_t bw = _mm_mul_ps(plane, VecShuffle(b, b, VecShuffleMask(3, 3, 3, 3 )));
        vec4_t aw = _mm_mul_ps(plane, VecShuffle(a, a, VecShuffleMask(3, 3, 3, 3 )));
		return _mm_add_ps(b, _mm_mul_ps(_mm_sub_ps(a, b), use_plane ? VecMul(_mm_sub_ps( bw, bz ), VecRcp(_mm_add_ps(_mm_sub_ps(az, bz), _mm_sub_ps(bw, aw)))) : VecMul(bz, VecRcp(_mm_sub_ps(bz, az)))));
	}

	template< int vertex_component, bool cmp_func > __forceinline int clip_triangle( const vec4_t* input, int count, vec4_t* output, vec4_t plane )
	{
		int vertices = 0;
		for ( int i = 0; i < count; ++i )
		{
            vec4_t a0 = input[ i*3 + 0 ];
            vec4_t b0 = input[ i*3 + 1 ];
            vec4_t c0 = input[ i*3 + 2 ];
            vec4_t tmp = VecShuffle( a0, b0, VecShuffleMask( vertex_component, 0, vertex_component, 0 ) );
            vec4_t w = VecShuffle( tmp, c0, VecShuffleMask( 0, 2, vertex_component, vertex_component ) );

			int mask = cmp_func ? _mm_movemask_ps( _mm_cmpgt_ps( w, _mm_mul_ps( VecShuffle( VecShuffle( a0, b0, VecShuffleMask( 3, 3, 3, 3 ) ), c0, VecShuffleMask( 0, 2, 3, 3 ) ), plane ) ) ) & 7 : _mm_movemask_ps( _mm_cmplt_ps( w, _mm_setzero_ps() ) ) & 7;
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

                    vec4_t ba = intersectLineZ< vertex_component, cmp_func >( c0, a0, plane );
                    vec4_t bb = intersectLineZ< vertex_component, cmp_func >( c0, b0, plane );
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

	__forceinline int clip_triangle( __m128* v, vec4_t* dst )
	{
        vec4_t g_total_width_v = Vector4(g_total_width);
        vec4_t g_total_height_v = Vector4(g_total_height);

		int count = 4;
        vec4_t input_array[ 196 ], output_array[ 196 ];
		count = clip_triangle< 1, false >( v, count, input_array, _mm_setzero_ps() ); // y < 0
		count = clip_triangle< 0, false >( input_array, count, output_array, _mm_setzero_ps() ); // x < 0
		count = clip_triangle< 0, true >( output_array, count, input_array, g_total_width_v ); // x > 1280
		count = clip_triangle< 2, false >( input_array, count, output_array, _mm_setzero_ps() ); // z < 0
		count = clip_triangle< 1, true >( output_array, count, input_array, g_total_height_v ); // y > 720

		int aligned_count = 3 *( ( count + 1 ) & ~1 );
		for ( int i = 0; i < aligned_count; ++i)
            dst[i] = input_array[i];

		return count;
	}

	__forceinline vec4_t get_tile_bounds( vec4_t v_min, vec4_t v_max )
	{
        vec4_t minmax = _mm_movelh_ps( v_min, v_max ); // xyXY
        vec4_t tile_bounds = _mm_mul_ps( minmax, g_tile_size ); // x/w y/h X/w Y/h
		tile_bounds = _mm_add_ps( tile_bounds, g_almost_one );
		return _mm_max_ps( _mm_min_ps( tile_bounds, g_tile_bounds ), _mm_setzero_ps() );
	}

    void sort_triangles(uint64_t* triangles, uint32_t size, stl::vector<uint64_t>& temp);

    template <bool is_occluder> __forceinline bool draw_scanlines(Tile& tile, int& xs1, int& xs2, int y1, int y2, int xa1, int xa2, const __m128i* masks, int* flag)
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

            __m128i span = VecIntXor(masks[xb], masks[xe]);
            if (is_occluder)
            {
                tile.m_frame_buffer[scanline] = VecIntOr(tile.m_frame_buffer[scanline], span);
                uint64_t bit = _mm_movemask_epi8(VecIntCmpEqual(VecIntAnd(tile.m_frame_buffer[scanline], m_full_span), m_full_span)) == 65535;
                tile.m_mask |= bit << scanline;
            }
            else
            {
                if (_mm_movemask_epi8(VecIntCmpEqual(VecIntAnd(tile.m_frame_buffer[scanline], span), span)) != 65535)
                {
                    *flag = 1;
                    return true;
                }
            }
        }
        return false;
    }

    template < bool is_occluder > __forceinline void draw_4triangles(Tile& tile, const Triangle& tri)
    {
        __m128 vx0 = tri.x0, vx1 = tri.x1, vx2 = tri.x2;
        __m128 vy0 = tri.y0, vy1 = tri.y1, vy2 = tri.y2;

        _MM_ALIGN16 int iy0[4], iy1[4], iy2[4], ix0[4], ix1[4], ix2[4], dx1[4], dx2[4], dx3[4];

        vec4_t vdx1 = _mm_mul_ps(_mm_sub_ps(vx2, vx0), _mm_rcp_ps(_mm_sub_ps(vy2, vy0)));
        vec4_t vdx2 = _mm_mul_ps(_mm_sub_ps(vx1, vx0), _mm_rcp_ps(_mm_sub_ps(vy1, vy0)));
        vec4_t vdx3 = _mm_mul_ps(_mm_sub_ps(vx2, vx1), _mm_rcp_ps(_mm_sub_ps(vy2, vy1)));
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
public:

    bool        m_skip_full = false;
    uint32_t    m_triangles_total = 0;
    uint32_t    m_triangles_occluder_total = 0;
    uint32_t    m_triangles_occludee_total = 0;
    uint32_t    m_triangles_drawn_total = 0;
    uint32_t    m_triangles_drawn_occluder_total = 0;
    uint32_t    m_triangles_drawn_occludee_total = 0;
    uint32_t    m_triangles_skipped = 0;
    uint32_t    m_triangles_offscreen = 0;

    Rasterizer()
    {
        m_tiles.reserve(g_width*g_height);
        for (int j = 0; j < g_height; ++j)
            for (int i = 0; i < g_width; ++i)
                m_tiles.push_back(Tile(i, j));

        m_full_span = Vector4Int(0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff);

        m_triangles.resize(g_max_4triangles);
        m_sort.reserve(128*1024);
    }

    void begin(const Matrix& m)
    {
        m_transform = m;

        m_triangle_count = 0;
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

    void push_box(const Matrix& mat, int* flag)
    {
        _MM_ALIGN16 int bounds_array[4] = { 0 };
        vec4_t m[ 16 ];
        ExtractMatrix(mat * m_transform, m);

        m_triangles_total += 12;
        if (flag)
        {
            *flag = 0;
            m_triangles_occluder_total += 12;
        }
        else
            m_triangles_occludee_total += 12;

        vec4_t g_total_width_v = Vector4(g_total_width);
        vec4_t g_total_height_v = Vector4(g_total_height);

        vec4_t pt[ 4 ];
        pt[0] = Vector4(-1.f, -1.f, 1.f, 1.f); // xxXX
        pt[1] = Vector4(-1.f, 1.f, -1.f, 1.f); // yYyY
        pt[2] = Vector4(-1.f); // zzzz
        pt[3] = Vector4(1.f); // ZZZZ

        vec4_t xxxx0 = VecMad(m[8], pt[2], VecMad(m[4], pt[1], VecMad(m[0], pt[0], m[12])));
        vec4_t yyyy0 = VecMad(m[9], pt[2], VecMad(m[5], pt[1], VecMad(m[1], pt[0], m[13])));
        vec4_t zzzz0 = VecMad(m[10], pt[2], VecMad(m[6], pt[1], VecMad(m[2], pt[0], m[14])));
        vec4_t wwww0 = VecMad(m[11], pt[2], VecMad(m[7], pt[1], VecMad(m[3], pt[0], m[15])));

        vec4_t xxxx1 = VecMad(m[8], pt[3], VecMad(m[4], pt[1], VecMad(m[0], pt[0], m[12])));
        vec4_t yyyy1 = VecMad(m[9], pt[3], VecMad(m[5], pt[1], VecMad(m[1], pt[0], m[13])));
        vec4_t zzzz1 = VecMad(m[10], pt[3], VecMad(m[6], pt[1], VecMad(m[2], pt[0], m[14])));
        vec4_t wwww1 = VecMad(m[11], pt[3], VecMad(m[7], pt[1], VecMad(m[3], pt[0], m[15])));

        vec4_t v_mask00 = _mm_and_ps(_mm_cmpgt_ps(xxxx0, _mm_setzero_ps()), _mm_cmpgt_ps(xxxx1, _mm_setzero_ps()));
        vec4_t v_mask01 = _mm_and_ps(_mm_cmpgt_ps(yyyy0, _mm_setzero_ps()), _mm_cmpgt_ps(yyyy1, _mm_setzero_ps()));
        vec4_t v_mask10 = _mm_and_ps(_mm_cmplt_ps(xxxx0, _mm_mul_ps(wwww0, g_total_width_v)), _mm_cmplt_ps(xxxx1, _mm_mul_ps(wwww1, g_total_width_v )));
        vec4_t v_mask11 = _mm_and_ps(_mm_cmplt_ps(yyyy0, _mm_mul_ps(wwww0, g_total_height_v)), _mm_cmplt_ps(yyyy1, _mm_mul_ps(wwww1, g_total_height_v)));

        vec4_t v_mask0 = _mm_and_ps(v_mask00, v_mask10);
        vec4_t v_mask1 = _mm_and_ps(v_mask01, v_mask11);
        int mask = _mm_movemask_ps(_mm_and_ps(v_mask0, v_mask1));

        bool intersect_near = _mm_movemask_ps(_mm_and_ps(_mm_cmpgt_ps(zzzz0, _mm_setzero_ps()), _mm_cmpgt_ps(zzzz1, _mm_setzero_ps()))) != 15;

        // fast path, no clipping
        if (mask == 15 && intersect_near == false)
        {
            vec4_t xx0 = VecMul(xxxx0, VecRcp(wwww0));
            vec4_t yy0 = VecMul(yyyy0, VecRcp(wwww0));
            vec4_t xx1 = VecMul(xxxx1, VecRcp(wwww1));
            vec4_t yy1 = VecMul(yyyy1, VecRcp(wwww1));

            vec4_t x_min = _mm_min_ps(xx0, xx1);
            vec4_t x_max = _mm_max_ps(xx0, xx1);
            vec4_t y_min = _mm_min_ps(yy0, yy1);
            vec4_t y_max = _mm_max_ps(yy0, yy1);

            vec4_t min_0 = _mm_min_ps(VecShuffle(x_min, y_min, VecShuffleMask(0, 1, 0, 1)), VecShuffle(x_min, y_min, VecShuffleMask(2, 3, 2, 3)));
            vec4_t max_0 = _mm_max_ps(VecShuffle(x_max, y_max, VecShuffleMask(0, 1, 0, 1)), VecShuffle(x_max, y_max, VecShuffleMask(2, 3, 2, 3)));

            vec4_t min_1 = _mm_min_ps(VecShuffle(min_0, min_0, VecShuffleMask(0, 2, 0, 0)), VecShuffle(min_0, min_0, VecShuffleMask(1, 3, 0, 0)));
            vec4_t max_1 = _mm_max_ps(VecShuffle(max_0, max_0, VecShuffleMask(0, 2, 0, 0)), VecShuffle(max_0, max_0, VecShuffleMask(1, 3, 0, 0)));
            VecIntStore(bounds_array, VecFloat2Int(get_tile_bounds(min_1, max_1)));

            bool fast_path = bounds_array[0] + 1 == bounds_array[2] && bounds_array[1] + 1 == bounds_array[3];
            {
                // 0 2 4 6   1 3 6 4   2 0 5 7
                vec4_t x0 = VecShuffle(xx0, xx1, VecShuffleMask(0, 2, 0, 2));
                vec4_t y0 = VecShuffle(yy0, yy1, VecShuffleMask(0, 2, 0, 2));
                vec4_t w0 = VecShuffle(wwww0, wwww1, VecShuffleMask(0, 2, 0, 2));

                vec4_t x1 = VecShuffle(xx0, xx1, VecShuffleMask(1, 3, 2, 0));
                vec4_t y1 = VecShuffle(yy0, yy1, VecShuffleMask(1, 3, 2, 0));
                vec4_t w1 = VecShuffle(wwww0, wwww1, VecShuffleMask(1, 3, 2, 0));

                vec4_t x2 = VecShuffle(xx0, xx1, VecShuffleMask(2, 0, 1, 3));
                vec4_t y2 = VecShuffle(yy0, yy1, VecShuffleMask(2, 0, 1, 3));
                vec4_t w2 = VecShuffle(wwww0, wwww1, VecShuffleMask(2, 0, 1, 3));

                if (fast_path)
                    push_4triangles<false>(flag, bounds_array, x0, y0, w0, x1, y1, w1, x2, y2, w2);
                else
                    push_4triangles<true>(flag, bounds_array, x0, y0, w0, x1, y1, w1, x2, y2, w2);
            }
            {
                // 2 7 0 5   7 2 5 0   3 6 1 4
                vec4_t xtmp0 = _mm_movehl_ps(xx1, xx0); // 2 3 6 7
                vec4_t xtmp1 = _mm_movelh_ps(xx0, xx1); // 0 1 4 5
                vec4_t ytmp0 = _mm_movehl_ps(yy1, yy0); // 2 3 6 7
                vec4_t ytmp1 = _mm_movelh_ps(yy0, yy1); // 0 1 4 5
                vec4_t wtmp0 = _mm_movehl_ps(wwww1, wwww0); // 2 3 6 7
                vec4_t wtmp1 = _mm_movelh_ps(wwww0, wwww1); // 0 1 4 5

                vec4_t x0 = VecShuffle(xtmp0, xtmp1, VecShuffleMask(0, 3, 0, 3));
                vec4_t y0 = VecShuffle(ytmp0, ytmp1, VecShuffleMask(0, 3, 0, 3));
                vec4_t w0 = VecShuffle(wtmp0, wtmp1, VecShuffleMask(0, 3, 0, 3));

                vec4_t x1 = VecShuffle(xtmp0, xtmp1, VecShuffleMask(3, 0, 3, 0));
                vec4_t y1 = VecShuffle(ytmp0, ytmp1, VecShuffleMask(3, 0, 3, 0));
                vec4_t w1 = VecShuffle(wtmp0, wtmp1, VecShuffleMask(3, 0, 3, 0));

                vec4_t x2 = VecShuffle(xtmp0, xtmp1, VecShuffleMask(1, 2, 1, 2));
                vec4_t y2 = VecShuffle(ytmp0, ytmp1, VecShuffleMask(1, 2, 1, 2));
                vec4_t w2 = VecShuffle(wtmp0, wtmp1, VecShuffleMask(1, 2, 1, 2));

                if (fast_path)
                    push_4triangles<false>(flag, bounds_array, x0, y0, w0, x1, y1, w1, x2, y2, w2);
                else
                    push_4triangles<true>(flag, bounds_array, x0, y0, w0, x1, y1, w1, x2, y2, w2);
            }
            {
                // 1 6 0 7   6 1 3 4   2 5 7 0
                vec4_t xtmp0 = _mm_movelh_ps(VecShuffle(xx0, xx0, VecShuffleMask(1, 2, 1, 2)), VecShuffle(xx1, xx1, VecShuffleMask(1, 2, 1, 2))); // 1 2 5 6
                vec4_t xtmp1 = _mm_movelh_ps(VecShuffle(xx0, xx0, VecShuffleMask(0, 3, 0, 3)), VecShuffle(xx1, xx1, VecShuffleMask(3, 0, 3, 0))); // 0 3 7 4
                vec4_t ytmp0 = _mm_movelh_ps(VecShuffle(yy0, yy0, VecShuffleMask(1, 2, 1, 2)), VecShuffle(yy1, yy1, VecShuffleMask(1, 2, 1, 2))); // 1 2 5 6
                vec4_t ytmp1 = _mm_movelh_ps(VecShuffle(yy0, yy0, VecShuffleMask(0, 3, 0, 3)), VecShuffle(yy1, yy1, VecShuffleMask(3, 0, 3, 0))); // 0 3 7 4
                vec4_t wtmp0 = _mm_movelh_ps(VecShuffle(wwww0, wwww0, VecShuffleMask(1, 2, 1, 2)), VecShuffle(wwww1, wwww1, VecShuffleMask(1, 2, 1, 2))); // 1 2 5 6
                vec4_t wtmp1 = _mm_movelh_ps(VecShuffle(wwww0, wwww0, VecShuffleMask(0, 3, 0, 3)), VecShuffle(wwww1, wwww1, VecShuffleMask(3, 0, 3, 0))); // 0 3 7 4

                vec4_t x0 = VecShuffle(xtmp0, xtmp1, VecShuffleMask(0, 3, 0, 2));
                vec4_t y0 = VecShuffle(ytmp0, ytmp1, VecShuffleMask(0, 3, 0, 2));
                vec4_t w0 = VecShuffle(wtmp0, wtmp1, VecShuffleMask(0, 3, 0, 2));

                vec4_t x1 = VecShuffle(xtmp0, xtmp1, VecShuffleMask(3, 0, 1, 3));
                vec4_t y1 = VecShuffle(ytmp0, ytmp1, VecShuffleMask(3, 0, 1, 3));
                vec4_t w1 = VecShuffle(wtmp0, wtmp1, VecShuffleMask(3, 0, 1, 3));

                vec4_t x2 = VecShuffle(xtmp0, xtmp1, VecShuffleMask(1, 2, 2, 0));
                vec4_t y2 = VecShuffle(ytmp0, ytmp1, VecShuffleMask(1, 2, 2, 0));
                vec4_t w2 = VecShuffle(wtmp0, wtmp1, VecShuffleMask(1, 2, 2, 0));

                if (fast_path)
                    push_4triangles<false>(flag, bounds_array, x0, y0, w0, x1, y1, w1, x2, y2, w2);
                else
                    push_4triangles<true>(flag, bounds_array, x0, y0, w0, x1, y1, w1, x2, y2, w2);
            }
        }
        else
        {
        }
    }

    void push_object(const Matrix& m, const vec4_t& v_min, const vec4_t& v_max, const uint16_t* indices, int index_count, vec4_t* vertices, int vertex_count, int* flag)
    {
        _MM_ALIGN16 int bounds_array[4] = { 0 };
        vec4_t matrix[ 16 ];
        ExtractMatrix(m * m_transform, matrix);

        if (m_mt)
            __atomic_add_fetch(&m_triangles_total, index_count / 3, __ATOMIC_SEQ_CST);
        else
            m_triangles_total += index_count / 3;
        if (flag)
        {
            *flag = 0;
            if (m_mt)
                __atomic_add_fetch(&m_triangles_occluder_total, index_count / 3, __ATOMIC_SEQ_CST);
            else
                m_triangles_occluder_total += index_count / 3;
        }
        else
        {
            if (m_mt)
                __atomic_add_fetch(&m_triangles_occludee_total, index_count / 3, __ATOMIC_SEQ_CST);
            else
                m_triangles_occludee_total += index_count / 3;
        }

        bool inside = occlude_object(matrix, v_min, v_max, bounds_array);
        if ( bounds_array[0] == bounds_array[2] || bounds_array[1] == bounds_array[3] )
        {
            m_triangles_offscreen += index_count / 3;
            return;
        }

        if ( inside )
        {
            constexpr uint32_t max_triangles_in_object = 1024;
            _MM_ALIGN16 vec4_t transformed_vertices[ max_triangles_in_object ];

            int aligned_count = ( vertex_count + 3 ) & ~3;
            assert( aligned_count < max_triangles_in_object );
            for ( int i = 0; i < aligned_count; i += 4 )
                Vector3TransformCoord4Homogeneous( matrix, vertices + i, transformed_vertices + i );

            if ( bounds_array[0] + 2 > bounds_array[2] && bounds_array[1] + 2 > bounds_array[3] )
                push_triangle_batched< false, true >( flag, transformed_vertices, index_count, indices, bounds_array );
            else
                push_triangle_batched< true, true >( flag, transformed_vertices, index_count, indices, bounds_array );
        }
        else
        {
            if ( bounds_array[0] + 2 > bounds_array[2] && bounds_array[1] + 2 > bounds_array[3] )
                push_object_clipped<false>( matrix, indices, index_count, vertices, vertex_count, bounds_array, flag );
            else
                push_object_clipped<true>( matrix, indices, index_count, vertices, vertex_count, bounds_array, flag );
        }
    }

    const __m128i* get_framebuffer(int tile) const
    {
        return m_tiles[tile].m_frame_buffer;
    }

    const stl::vector<Tile>& get_tiles() const
    {
        return m_tiles;
    }

    void sort_triangles(uint32_t t, stl::vector<uint64_t>& temp)
    {
        auto & tile = m_tiles[t];
        if (tile.m_triangle_count)
            sort_triangles((uint64_t*)tile.m_triangles.data(), tile.m_triangle_count, temp);
    }

    void sort_triangles()
    {
        for (auto & tile : m_tiles)
            sort_triangles((uint64_t*)tile.m_triangles.data(), tile.m_triangle_count, m_sort);
    }

    void draw_triangles(uint32_t t)
    {
        auto & tile = m_tiles[t];

        assert(tile.m_triangle_count <= m_triangle_count);
        const uint64_t* tri = tile.m_triangles.data(), *tri_end = tri + tile.m_triangle_count;
        const Triangle* tri_data = m_triangles.data();
        while (tri != tri_end)
        {
            assert((*tri>>32) < m_triangle_count);
            auto & triangle = tri_data[*tri++>>32];
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
                        tile.m_triangles_skipped += bx::uint32_cntbits(tri_data[*tri++>>32].mask);
                    break;
                }
            }
        }
        if (m_mt)
        {
            __atomic_add_fetch(&m_triangles_drawn_total, tile.m_triangles_drawn_total, __ATOMIC_SEQ_CST);
            __atomic_add_fetch(&m_triangles_drawn_occluder_total, tile.m_triangles_drawn_occluder_total, __ATOMIC_SEQ_CST);
            __atomic_add_fetch(&m_triangles_drawn_occludee_total, tile.m_triangles_drawn_occludee_total, __ATOMIC_SEQ_CST);
            __atomic_add_fetch(&m_triangles_skipped, tile.m_triangles_skipped, __ATOMIC_SEQ_CST);
        }
        else
        {
            m_triangles_drawn_total += tile.m_triangles_drawn_total;
            m_triangles_drawn_occluder_total += tile.m_triangles_drawn_occluder_total;
            m_triangles_drawn_occludee_total += tile.m_triangles_drawn_occludee_total;
            m_triangles_skipped += tile.m_triangles_skipped;
        }
    }

    void draw_triangles()
    {
        for (auto & tile : m_tiles)
        {
            assert(tile.m_triangle_count <= m_triangle_count);
            const uint64_t* tri = tile.m_triangles.data(), *tri_end = tri + tile.m_triangle_count;
            const Triangle* tri_data = m_triangles.data();
            while (tri != tri_end)
            {
                assert((*tri>>32) < m_triangle_count);
                auto & triangle = tri_data[*tri++>>32];
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
                            tile.m_triangles_skipped += bx::uint32_cntbits(tri_data[*tri++>>32].mask);
                        break;
                    }
                }
            }
            m_triangles_drawn_total += tile.m_triangles_drawn_total;
            m_triangles_drawn_occluder_total += tile.m_triangles_drawn_occluder_total;
            m_triangles_drawn_occludee_total += tile.m_triangles_drawn_occludee_total;
            m_triangles_skipped += tile.m_triangles_skipped;
        }
    }

    void setMT(bool mt)
    {
        m_mt = mt;
    }
    
	static void Init();
};
