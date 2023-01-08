#pragma once

#include <assert.h>
#include "rasterizer_math.h"

#include <tinystl/allocator.h>
#include <tinystl/vector.h>

namespace stl = tinystl;

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

struct Rasterizer;

struct _MM_ALIGN16 Tile
{
    static const int        g_tile_height = 720;
    static const int        g_tile_width = 128;

    __m128i                 m_frame_buffer[ g_tile_height ];
    stl::vector<uint64_t>   m_triangles;
    stl::vector<__m128i>    m_shifts;
    Rasterizer*             m_rasterizer = nullptr;
    uint32_t                m_triangles_drawn_total = 0;
    uint32_t                m_triangles_drawn_occluder_total = 0;
    uint32_t                m_triangles_drawn_occludee_total = 0;
    uint32_t                m_triangles_skipped = 0;

    Tile(Rasterizer* rasterizer, int x);

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

    template <bool is_occluder> __forceinline bool draw_scanlines(int& xs1, int xs2, int y1, int y2, int xa1, int xa2, const __m128i* masks, int* flag)
    {
        assert((is_occluder && !flag) || (flag && !is_occluder));
        for (int scanline = y1; scanline < y2; ++scanline)
        {
            int xb = xs1 >> 16;
            int xe = xs2 >> 16;
            xs1 += xa1;
            xs2 += xa2;

            __m128i span = VecIntXor(masks[xb], masks[xe]);
            if (is_occluder)
                m_frame_buffer[scanline] = VecIntOr(m_frame_buffer[scanline], span);
            else
            {
                if (_mm_movemask_epi8(VecIntCmpEqual(VecIntAnd(m_frame_buffer[scanline], span), span)) != 65535)
                {
                    *flag = 1;
                    return true;
                }
            }
        }
        return false;
    }
		
	template < bool is_occluder > __forceinline void draw_4triangles(const Triangle& tri)
	{
		__m128 vx0 = tri.x0, vx1 = tri.x1, vx2 = tri.x2;
		__m128 vy0 = tri.y0, vy1 = tri.y1, vy2 = tri.y2;

		sort(vx0, vx1, vy0, vy1);
		sort(vx1, vx2, vy1, vy2);
		sort(vx0, vx1, vy0, vy1);

		_MM_ALIGN16 int iy0[4], iy1[4], iy2[4], ix0[4], ix1[4], dx1[4], dx2[4], dx3[4];
		VecIntStore(iy0, VecFloat2Int(vy0));
		VecIntStore(iy1, VecFloat2Int(vy1));
		VecIntStore(iy2, VecFloat2Int(vy2));
		VecIntStore(ix0, VecFloat2Int(vx0));
		VecIntStore(ix1, VecFloat2Int(vx1));

		VecIntStore(dx1, VecFloat2Int(_mm_mul_ps(_mm_sub_ps(vx2, vx0), _mm_rcp_ps(_mm_sub_ps(vy2, vy0)))));
		VecIntStore(dx2, VecFloat2Int(_mm_mul_ps(_mm_sub_ps(vx1, vx0), _mm_rcp_ps(_mm_sub_ps(vy1, vy0)))));
		VecIntStore(dx3, VecFloat2Int(_mm_mul_ps(_mm_sub_ps(vx2, vx1), _mm_rcp_ps(_mm_sub_ps(vy2, vy1)))));

        const __m128i span = Vector4Int(0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff);
		for (size_t i = 0, mask = 1; i < 4; ++i, mask <<= 1)
		{
			if ((tri.mask & mask ) == 0)
				continue;

            bool skip = true;
            for (int y = iy0[i]; y < iy2[i]; ++y)
            {
                skip &= _mm_movemask_epi8(VecIntCmpEqual(VecIntAnd(m_frame_buffer[y], span), span)) == 65535;
                if (!skip)
                    break;
            }
            if (skip)
            {
                ++m_triangles_skipped;
                return;
            }

            m_triangles_drawn_total++;
            if (is_occluder)
                m_triangles_drawn_occluder_total++;
            else
                m_triangles_drawn_occludee_total++;

			int xs1 = ix0[i], xs2 = ix0[i], xs3 = ix1[i];
            if (draw_scanlines<is_occluder>(xs1, xs2, iy0[i], iy1[i], dx1[i], dx2[i], m_shifts.data(), tri.flag))
                return;
            if (draw_scanlines<is_occluder>(xs1, xs3, iy1[i], iy2[i], dx1[i], dx3[i], m_shifts.data(), tri.flag))
                return;
		}
	}

    void sort_triangles();
    void draw_triangles();
};
