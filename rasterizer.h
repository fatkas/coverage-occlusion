#pragma once

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "rasterizer_math.h"

#include "rasterizer_tile.h"

struct _MM_ALIGN16 Rasterizer
{
	static const int g_width = 10;
	static const int g_total_width = g_width * Tile::g_tile_width;
	static const int g_total_height = Tile::g_tile_height;

    stl::vector<Tile>       m_tiles;
	Matrix			        m_transform;
    stl::vector<Triangle>   m_triangles;
    stl::vector<uint64_t>   m_sort;

    uint32_t                m_triangles_total = 0;
    uint32_t                m_triangles_occluder_total = 0;
    uint32_t                m_triangles_occludee_total = 0;
    uint32_t                m_triangles_drawn_total = 0;
    uint32_t                m_triangles_drawn_occluder_total = 0;
    uint32_t                m_triangles_drawn_occludee_total = 0;
    uint32_t                m_triangles_skipped = 0;

	Rasterizer( const Matrix& m ) : m_transform( m )
	{
        for (int i = 0; i < g_width; ++i)
            m_tiles.push_back(Tile(this, i));

        m_triangles.reserve(128*1024);
        m_sort.reserve(128*1024);
	}

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

	template < bool select_tiles, bool use_indices > inline void push_triangle_batched( int* flag, float* z, const vec2_t* src, int count, const unsigned short* indices, int* bounds_array )
	{
		assert( ( (count / 3) & 3 ) == 0 );

		_MM_ALIGN16 int transformed_bounds[ 4 ];

		const vec4_t local_fixed_point = Vector4(65536.f);
		for ( int i = 0; i < count; i += 12 )
		{
            Triangle t;

            t.flag = flag;

            float zz = flag ? z[0] : z[1];
            assert(zz >= 0.f && zz <= 65535.f);
            t.z = int(zz*256.f);

			#define IDX(num)( use_indices ? indices[ i + num ] : i + num )
            vec4_t v0_0 = VecLoadU( src + IDX(0) );
            vec4_t v0_1 = VecLoadU( src + IDX(3) );
            vec4_t v0_2 = VecLoadU( src + IDX(6) );
            vec4_t v0_3 = VecLoadU( src + IDX(9) );
            vec4_t tmp0 = _mm_unpacklo_ps( v0_0, v0_1 ); // x0_0 x1_0 y0_0 y1_0
            vec4_t tmp1 = _mm_unpacklo_ps( v0_2, v0_3 ); // x2_0 x3_0 y2_0 y3_0
			t.x0 = _mm_mul_ps( _mm_movelh_ps( tmp0, tmp1 ), local_fixed_point );
			t.y0 = _mm_movehl_ps( tmp1, tmp0 );
			
            vec4_t v1_0 = VecLoadU( src + IDX(1) );
            vec4_t v1_1 = VecLoadU( src + IDX(4) );
            vec4_t v1_2 = VecLoadU( src + IDX(7) );
            vec4_t v1_3 = VecLoadU( src + IDX(10) );
            vec4_t tmp2 = _mm_unpacklo_ps( v1_0, v1_1 ); // x0_1 x1_1 y0_1 y1_1
            vec4_t tmp3 = _mm_unpacklo_ps( v1_2, v1_3 ); // x2_1 x3_1 y2_1 y3_1
			t.x1 = _mm_mul_ps( _mm_movelh_ps( tmp2, tmp3 ), local_fixed_point );
			t.y1 = _mm_movehl_ps( tmp3, tmp2 );

            vec4_t v2_0 = VecLoadU( src + IDX(2) );
            vec4_t v2_1 = VecLoadU( src + IDX(5) );
            vec4_t v2_2 = VecLoadU( src + IDX(8) );
            vec4_t v2_3 = VecLoadU( src + IDX(11) );
            vec4_t tmp4 = _mm_unpacklo_ps( v2_0, v2_1 ); // x0_2 x1_2 y0_2 y1_2
            vec4_t tmp5 = _mm_unpacklo_ps( v2_2, v2_3 ); // x2_2 x3_2 y2_2 y3_2
			t.x2 = _mm_mul_ps( _mm_movelh_ps( tmp4, tmp5 ), local_fixed_point );
			t.y2 = _mm_movehl_ps( tmp5, tmp4 );
			#undef IDX			
	
			t.mask = _mm_movemask_ps( _mm_cmplt_ps( _mm_sub_ps( _mm_mul_ps( _mm_sub_ps( t.x1, t.x0 ), _mm_sub_ps( t.y2, t.y0 ) ), _mm_mul_ps( _mm_sub_ps( t.x2, t.x0 ), _mm_sub_ps( t.y1, t.y0 ) ) ), _mm_setzero_ps() ) );

			if ( t.mask == 0 )
				continue;

            uint64_t index = m_triangles.size();
            m_triangles.push_back(t);

			if ( select_tiles )
			{			
                vec4_t min_bound = _mm_min_ps( _mm_min_ps( _mm_min_ps( _mm_min_ps( v0_0, v1_0 ), v2_0 ), _mm_min_ps( _mm_min_ps( v0_1, v1_1 ), v2_1 ) ), _mm_min_ps( _mm_min_ps( _mm_min_ps( v0_2, v1_2 ), v2_2 ), _mm_min_ps( _mm_min_ps( v0_3, v1_3 ), v2_3 ) ) );
                vec4_t max_bound = _mm_max_ps( _mm_max_ps( _mm_max_ps( _mm_max_ps( v0_0, v1_0 ), v2_0 ), _mm_max_ps( _mm_max_ps( v0_1, v1_1 ), v2_1 ) ), _mm_max_ps( _mm_max_ps( _mm_max_ps( v0_2, v1_2 ), v2_2 ), _mm_max_ps( _mm_max_ps( v0_3, v1_3 ), v2_3 ) ) );

                vec4_t tile_bound = get_tile_bounds(min_bound, max_bound);
				VecIntStore( transformed_bounds, VecFloat2Int( tile_bound ) );
				assert( transformed_bounds[0] >= bounds_array[0] );
				assert( transformed_bounds[1] >= bounds_array[1] );
				assert( transformed_bounds[2] <= bounds_array[2] );
				assert( transformed_bounds[3] <= bounds_array[3] );

				for ( int x = transformed_bounds[0]; x < transformed_bounds[2]; ++x )
                    m_tiles[x].m_triangles.push_back((index<<32)|t.z);
			}
			else
			{
                assert(bounds_array[0] < g_width);
                m_tiles[bounds_array[0]].m_triangles.push_back((index<<32)|t.z);
			}
		}
	}

	bool occlude_object( const __m128* m, const __m128& v_min, const __m128& v_max, int* bounds_array, float* z )
	{
        vec4_t g_total_width_v = Vector4(g_total_width);
        vec4_t g_total_height_v = Vector4(g_total_height);

		#define INTERSECT_EDGE(a,b,c) _mm_add_ps( b, _mm_mul_ps( _mm_sub_ps( a, b ), t ) )
        vec4_t pt[ 4 ];
        vec4_t vTmp = _mm_unpacklo_ps( v_min, v_max );				// x, X, y, Y
		pt[0] = VecShuffle( v_min, v_max, _MM_SHUFFLE( 0, 0, 0, 0 ) ); // xxXX	
		pt[1] = VecShuffle( vTmp, vTmp, _MM_SHUFFLE( 2, 3, 2, 3 ) ); // yYyY
		pt[2] = VecShuffle( v_min, v_min, _MM_SHUFFLE( 2, 2, 2, 2 ) ); // zzzz
		pt[3] = VecShuffle( v_max, v_max, _MM_SHUFFLE( 2, 2, 2, 2 ) ); // ZZZZ

        vec4_t xxxx0 = VecMad( m[8], pt[2], VecMad( m[4], pt[1], VecMad( m[0], pt[0], m[12] ) ) );
        vec4_t yyyy0 = VecMad( m[9], pt[2], VecMad( m[5], pt[1], VecMad( m[1], pt[0], m[13] ) ) );
        vec4_t zzzz0 = VecMad( m[10], pt[2], VecMad( m[6], pt[1], VecMad( m[2], pt[0], m[14] ) ) );
        vec4_t wwww0 = VecMad( m[11], pt[2], VecMad( m[7], pt[1], VecMad( m[3], pt[0], m[15] ) ) );

        vec4_t xxxx1 = VecMad( m[8], pt[3], VecMad( m[4], pt[1], VecMad( m[0], pt[0], m[12] ) ) );
        vec4_t yyyy1 = VecMad( m[9], pt[3], VecMad( m[5], pt[1], VecMad( m[1], pt[0], m[13] ) ) );
        vec4_t zzzz1 = VecMad( m[10], pt[3], VecMad( m[6], pt[1], VecMad( m[2], pt[0], m[14] ) ) );
        vec4_t wwww1 = VecMad( m[11], pt[3], VecMad( m[7], pt[1], VecMad( m[3], pt[0], m[15] ) ) );

        vec4_t v_mask00 = _mm_and_ps( _mm_cmpgt_ps( xxxx0, _mm_setzero_ps() ), _mm_cmpgt_ps( xxxx1, _mm_setzero_ps() ) );
        vec4_t v_mask01 = _mm_and_ps( _mm_cmpgt_ps( yyyy0, _mm_setzero_ps() ), _mm_cmpgt_ps( yyyy1, _mm_setzero_ps() ) );
        vec4_t v_mask10 = _mm_and_ps( _mm_cmplt_ps( xxxx0, _mm_mul_ps( wwww0, g_total_width_v ) ), _mm_cmplt_ps( xxxx1, _mm_mul_ps( wwww1, g_total_width_v ) ) );
        vec4_t v_mask11 = _mm_and_ps( _mm_cmplt_ps( yyyy0, _mm_mul_ps( wwww0, g_total_height_v ) ), _mm_cmplt_ps( yyyy1, _mm_mul_ps( wwww1, g_total_height_v ) ) );

        vec4_t v_mask0 = _mm_and_ps( v_mask00, v_mask10 );
        vec4_t v_mask1 = _mm_and_ps( v_mask01, v_mask11 );
		int mask = _mm_movemask_ps( _mm_and_ps( v_mask0, v_mask1 ) );

        vec4_t min_w = _mm_min_ps(wwww0, wwww1);
        vec4_t max_w = _mm_max_ps(wwww0, wwww1);
		min_w = _mm_min_ps( min_w, VecShuffle( min_w, min_w, VecShuffleMask( 2, 3, 2, 3 ) ) );
		min_w = _mm_min_ps( min_w, VecShuffle( min_w, min_w, VecShuffleMask( 1, 1, 1, 1 ) ) );
		max_w = _mm_max_ps( max_w, VecShuffle( max_w, max_w, VecShuffleMask( 2, 3, 2, 3 ) ) );		
		max_w = _mm_max_ps( max_w, VecShuffle( max_w, max_w, VecShuffleMask( 1, 1, 1, 1 ) ) );
		_mm_store_ps( z, _mm_unpacklo_ps( min_w, max_w ) );

		bool intersect_near = _mm_movemask_ps( _mm_and_ps( _mm_cmpgt_ps( zzzz0, _mm_setzero_ps() ), _mm_cmpgt_ps( zzzz1, _mm_setzero_ps() ) ) ) != 15;

        vec4_t x_min, x_max, y_min, y_max;
		if ( intersect_near == false )
		{
            vec4_t x0 = _mm_div_ps( xxxx0, wwww0 );
            vec4_t y0 = _mm_div_ps( yyyy0, wwww0 );
            vec4_t x1 = _mm_div_ps( xxxx1, wwww1 );
            vec4_t y1 = _mm_div_ps( yyyy1, wwww1 );

			x_min = _mm_min_ps( x0, x1 );
			x_max = _mm_max_ps( x0, x1 );
			y_min = _mm_min_ps( y0, y1 );
			y_max = _mm_max_ps( y0, y1 );
		}
		else
		{
            vec4_t xxxx0_1 = VecShuffle( xxxx0, xxxx0, VecShuffleMask( 1, 3, 0, 2 ) );
            vec4_t yyyy0_1 = VecShuffle( yyyy0, yyyy0, VecShuffleMask( 1, 3, 0, 2 ) );
            vec4_t zzzz0_1 = VecShuffle( zzzz0, zzzz0, VecShuffleMask( 1, 3, 0, 2 ) );
            vec4_t wwww0_1 = VecShuffle( wwww0, wwww0, VecShuffleMask( 1, 3, 0, 2 ) );

            vec4_t xxxx1_1 = VecShuffle( xxxx1, xxxx1, VecShuffleMask( 1, 3, 0, 2 ) );
            vec4_t yyyy1_1 = VecShuffle( yyyy1, yyyy1, VecShuffleMask( 1, 3, 0, 2 ) );
            vec4_t zzzz1_1 = VecShuffle( zzzz1, zzzz1, VecShuffleMask( 1, 3, 0, 2 ) );
            vec4_t wwww1_1 = VecShuffle( wwww1, wwww1, VecShuffleMask( 1, 3, 0, 2 ) );

            vec4_t t = _mm_div_ps( zzzz1, _mm_sub_ps( zzzz1, zzzz0 ) );
            vec4_t new_xxxx0 = INTERSECT_EDGE( xxxx0, xxxx1, t );
            vec4_t new_yyyy0 = INTERSECT_EDGE( yyyy0, yyyy1, t );
            vec4_t new_wwww0 = INTERSECT_EDGE( wwww0, wwww1, t );

			t = _mm_div_ps( zzzz0_1, _mm_sub_ps( zzzz0_1, zzzz0 ) );
            vec4_t new_xxxx1 = INTERSECT_EDGE( xxxx0, xxxx0_1, t );
            vec4_t new_yyyy1 = INTERSECT_EDGE( yyyy0, yyyy0_1, t );
            vec4_t new_wwww1 = INTERSECT_EDGE( wwww0, wwww0_1, t );

			t = _mm_div_ps( zzzz1_1, _mm_sub_ps( zzzz1_1, zzzz0 ) );
            vec4_t new_xxxx2 = INTERSECT_EDGE( xxxx1, xxxx1_1, t );
            vec4_t new_yyyy2 = INTERSECT_EDGE( yyyy1, yyyy1_1, t );
            vec4_t new_wwww2 = INTERSECT_EDGE( wwww1, wwww1_1, t );

            vec4_t x0 = _mm_div_ps( xxxx0, wwww0 );
            vec4_t y0 = _mm_div_ps( yyyy0, wwww0 );
            vec4_t x1 = _mm_div_ps( xxxx1, wwww1 );
            vec4_t y1 = _mm_div_ps( yyyy1, wwww1 );
            vec4_t x2 = _mm_div_ps( new_xxxx0, new_wwww0 );
            vec4_t y2 = _mm_div_ps( new_yyyy0, new_wwww0 );
            vec4_t x3 = _mm_div_ps( new_xxxx1, new_wwww1 );
            vec4_t y3 = _mm_div_ps( new_yyyy1, new_wwww1 );
            vec4_t x4 = _mm_div_ps( new_xxxx2, new_wwww2 );
            vec4_t y4 = _mm_div_ps( new_yyyy2, new_wwww2 );

			x_min = _mm_min_ps( _mm_min_ps( x0, x1 ), _mm_min_ps( x2, _mm_min_ps( x3, x4 ) ) );
			x_max = _mm_max_ps( _mm_max_ps( x0, x1 ), _mm_max_ps( x2, _mm_max_ps( x3, x4 ) ) );
			y_min = _mm_min_ps( _mm_min_ps( y0, y1 ), _mm_min_ps( y2, _mm_min_ps( y3, y4 ) ) );
			y_max = _mm_max_ps( _mm_max_ps( y0, y1 ), _mm_max_ps( y2, _mm_max_ps( y3, y4 ) ) );
		}

        vec4_t min_0 = _mm_min_ps( VecShuffle( x_min, y_min, VecShuffleMask( 0, 1, 0, 1 ) ), VecShuffle( x_min, y_min, VecShuffleMask( 2, 3, 2, 3 ) ) );
        vec4_t max_0 = _mm_max_ps( VecShuffle( x_max, y_max, VecShuffleMask( 0, 1, 0, 1 ) ), VecShuffle( x_max, y_max, VecShuffleMask( 2, 3, 2, 3 ) ) );

        vec4_t min_1 = _mm_min_ps( VecShuffle( min_0, min_0, VecShuffleMask( 0, 2, 0, 0 ) ), VecShuffle( min_0, min_0, VecShuffleMask( 1, 3, 0, 0 ) ) );
        vec4_t max_1 = _mm_max_ps( VecShuffle( max_0, max_0, VecShuffleMask( 0, 2, 0, 0 ) ), VecShuffle( max_0, max_0, VecShuffleMask( 1, 3, 0, 0 ) ) );

		VecIntStore( bounds_array, VecFloat2Int( get_tile_bounds( min_1, max_1 ) ) );
		return mask == 15 && intersect_near == false;

		#undef INTERSECT_EDGE
	}

	void push_object( const Matrix& m, const vec4_t& v_min, const vec4_t& v_max, const uint16_t* indices, int index_count, vec4_t* vertices, int vertex_count, int* flag )
	{
		_MM_ALIGN16 int bounds_array[4];
        vec4_t matrix[ 16 ];
		ExtractMatrix(m * m_transform, matrix);

        m_triangles_total += index_count / 3;
        if (flag)
        {
            *flag = 0;
            m_triangles_occluder_total += index_count / 3;
        }
        else
            m_triangles_occludee_total += index_count / 3;

		_MM_ALIGN16 float z[4];
		bool inside = occlude_object(matrix, v_min, v_max, bounds_array, z);
		if ( bounds_array[0] == bounds_array[2] || bounds_array[1] == bounds_array[3] )
            return;

		if ( inside )
		{
            constexpr uint32_t max_triangles_in_object = 1024;
			_MM_ALIGN16 vec2_t transformed_vertices[ max_triangles_in_object ];

			int aligned_count = ( vertex_count + 3 ) & ~3;
			assert( aligned_count < max_triangles_in_object );
			for ( int i = 0; i < aligned_count; i += 4 )
				Vector3TransformCoord4( matrix, vertices + i, transformed_vertices + i );

			if ( bounds_array[0] + 2 > bounds_array[2] && bounds_array[1] + 2 > bounds_array[3] )
				push_triangle_batched< false, true >( flag, z, transformed_vertices, index_count, indices, bounds_array );
			else
				push_triangle_batched< true, true >( flag, z, transformed_vertices, index_count, indices, bounds_array );
		}
		else
		{
			if ( bounds_array[0] + 2 > bounds_array[2] && bounds_array[1] + 2 > bounds_array[3] )
				push_object_clipped<false>( matrix, z, indices, index_count, vertices, vertex_count, bounds_array, flag );
			else
				push_object_clipped<true>( matrix, z, indices, index_count, vertices, vertex_count, bounds_array, flag );
		}
	}

	template < bool select_tiles > void push_object_clipped( const vec4_t* matrix, float* z, const uint16_t* indices, int index_count, vec4_t* vertices, int vertex_count, int* bounds_array, int* flag )
	{
        constexpr uint32_t max_triangles_in_object = 1024;
		_MM_ALIGN16 vec4_t transformed_vertices[ max_triangles_in_object ];
		_MM_ALIGN16 vec2_t clipped_triangles[ max_triangles_in_object*3 ];

		int aligned_count = ( vertex_count + 3 ) & ~3;
		assert( aligned_count < max_triangles_in_object );
		for ( int i = 0; i < aligned_count; i += 4 )
			Vector3TransformCoord4Homogeneous( matrix, vertices + i, transformed_vertices + i );

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
		vec2_t& v0 = clipped_triangles[ clipped_triangle_count * 3 - 3 ];
		vec2_t& v1 = clipped_triangles[ clipped_triangle_count * 3 - 2 ];
		vec2_t& v2 = clipped_triangles[ clipped_triangle_count * 3 - 1 ];
		for ( int i = 0; i < tris_to_pad; ++i, ++clipped_triangle_count )
		{
			clipped_triangles[ clipped_triangle_count*3 + 0 ] = v0;
			clipped_triangles[ clipped_triangle_count*3 + 1 ] = v1;
			clipped_triangles[ clipped_triangle_count*3 + 2 ] = v2;
		}
		push_triangle_batched< select_tiles, false >( flag, z, clipped_triangles, clipped_triangle_count*3, 0, bounds_array );
	}

	template< int v, bool use_plane > __forceinline vec4_t intersectLineZ( vec4_t a, vec4_t b, vec4_t plane )
	{		
		//  t = (a.x - a.w) / (b.w - a.w - b.x + a.x);
        vec4_t bz = VecShuffle( b, b, VecShuffleMask( v, v, v, v ) );
        vec4_t az = VecShuffle( a, a, VecShuffleMask( v, v, v, v ) );
        vec4_t bw = _mm_mul_ps( plane, VecShuffle( b, b, VecShuffleMask( 3, 3, 3, 3 ) ) );
        vec4_t aw = _mm_mul_ps( plane, VecShuffle( a, a, VecShuffleMask( 3, 3, 3, 3 ) ) );
		return _mm_add_ps( b, _mm_mul_ps( _mm_sub_ps( a, b ), use_plane ? _mm_div_ps( _mm_sub_ps( bw, bz ), _mm_add_ps( _mm_sub_ps( az, bz ), _mm_sub_ps( bw, aw ) ) ) : _mm_div_ps( bz, _mm_sub_ps( bz, az ) ) ) );
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

	__forceinline int clip_triangle( __m128* v, vec2_t* dst )
	{
        vec4_t g_total_width_v = Vector4(g_total_width);
        vec4_t g_total_height_v = Vector4(g_total_height);

		#define PROJECT(v0, v1) _mm_div_ps( _mm_movelh_ps( v0, v1 ), ( VecShuffle( v0, v1, VecShuffleMask( 3, 3, 3, 3 ) ) ) )
		int count = 4;
        vec4_t input_array[ 196 ], output_array[ 196 ];
		count = clip_triangle< 1, false >( v, count, input_array, _mm_setzero_ps() ); // y < 0
		count = clip_triangle< 0, false >( input_array, count, output_array, _mm_setzero_ps() ); // x < 0
		count = clip_triangle< 0, true >( output_array, count, input_array, g_total_width_v ); // x > 1280
		count = clip_triangle< 2, false >( input_array, count, output_array, _mm_setzero_ps() ); // z < 0
		count = clip_triangle< 1, true >( output_array, count, input_array, g_total_height_v ); // y > 720

		int aligned_count = 3 *( ( count + 1 ) & ~1 );
		for ( int i = 0; i < aligned_count; i += 6 )
		{
			VecStoreU( dst + i + 0, PROJECT( input_array[ i + 0 ], input_array[ i + 1 ] ) );
			VecStoreU( dst + i + 2, PROJECT( input_array[ i + 2 ], input_array[ i + 3 ] ) );
			VecStoreU( dst + i + 4, PROJECT( input_array[ i + 4 ], input_array[ i + 5 ] ) );
		}
		return count;

		#undef PROJECT
	}

	__forceinline vec4_t get_tile_bounds( vec4_t v_min, vec4_t v_max )
	{
        vec4_t tile_size = Vector4(1.f / (float)Tile::g_tile_width, 1.f / (float)Tile::g_tile_height, 1.f / (float)Tile::g_tile_width, 1.f / (float)Tile::g_tile_height);
        vec4_t almost_one = Vector4(0.f, 0.f, 0.9999f, 0.9999f);
        vec4_t g_tile_bounds = Vector4((float)Rasterizer::g_width, 1.f, (float)Rasterizer::g_width, 1.f);

        vec4_t minmax = _mm_movelh_ps( v_min, v_max ); // xyXY
        vec4_t tile_bounds = _mm_mul_ps( minmax, tile_size ); // x/w y/h X/w Y/h
		tile_bounds = _mm_add_ps( tile_bounds, almost_one );
		return _mm_max_ps( _mm_min_ps( tile_bounds, g_tile_bounds ), _mm_setzero_ps() );
	}

    void sort_triangles(int tile)
    {
        m_tiles[tile].sort_triangles();
    }

    void draw_triangles(int tile)
    {
        m_tiles[tile].draw_triangles();
        m_triangles_drawn_total += m_tiles[tile].m_triangles_drawn_total;
        m_triangles_drawn_occluder_total += m_tiles[tile].m_triangles_drawn_occluder_total;
        m_triangles_drawn_occludee_total += m_tiles[tile].m_triangles_drawn_occludee_total;
        m_triangles_skipped += m_tiles[tile].m_triangles_skipped;
    }
    
	static void Init();
};
