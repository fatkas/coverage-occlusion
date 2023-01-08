
#include <assert.h>
#include "rasterizer.h"
#include "rasterizer_tile.h"

#include <algorithm>

inline static int Shift( int val )
{
    if( val > 31 )
        return 0;
    if( val < 0 )
        return 0xffffffff;
    return 0xffffffff >> val;
}

Tile::Tile(Rasterizer* rasterizer, int x)
    : m_rasterizer(rasterizer)
{
    for ( int i = 0; i < g_tile_height; ++i )
        m_frame_buffer[i] = VecIntZero();

    m_shifts.resize(4096);
    for (int i = 0; i < x*128; ++i )
        m_shifts[i] = Vector4Int( 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff );
    for (int i = (x+1)*128; i < 4096; ++i )
        m_shifts[i] = VecIntZero();
    for (int i = 0; i < 128; ++i )
        m_shifts[x*128 + i] = _mm_set_epi32( Shift(i), Shift(i - 32), Shift(i - 64), Shift(i - 96) );

    m_triangles.reserve(128*1024);
}

void Tile::draw_triangles()
{
    for (auto & tri : m_triangles)
    {
        auto & triangle = m_rasterizer->m_triangles[tri>>32];
        if (triangle.flag)
            draw_4triangles<false>(triangle);
        else
            draw_4triangles<true>(triangle);
    }
}

void Tile::sort_triangles()
{
#if 0
    // const Triangle* tris = m_rasterizer->m_triangles.data();
    std::sort(m_triangles.begin(), m_triangles.end(), [](uint64_t a, uint64_t b)
    {
        return (a&0xffffffff) < (b&0xffffffff);
        // return tris[a].z < tris[b].z;
    });
#else
    constexpr unsigned max_bytes = 3;
    
	unsigned int histograms[256*max_bytes];

    memset(histograms, 0, sizeof(histograms));

	unsigned int* h0 = histograms;
	unsigned int* h1 = histograms + 256;
	unsigned int* h2 = histograms + 256*2;
	unsigned int* h3 = histograms + 256*3;

    m_rasterizer->m_sort.resize(m_triangles.size());

    // const Triangle* tris = m_rasterizer->m_triangles.data();
    uint64_t* src = m_triangles.data(), *src_end = m_triangles.data() + m_triangles.size(), *dst = m_rasterizer->m_sort.data();
    while (src < src_end)
	{
#define _0(h) ((h) & 255)
#define _1(h) (((h) >> 8) & 255)
#define _2(h) (((h) >> 16) & 255)
#define _3(h) ((h) >> 24)

		unsigned int h = (*src) & 0xffffffff;
		h0[_0(h)]++; 
		h1[_1(h)]++; 
		h2[_2(h)]++;
        if (max_bytes>3)
        {
            h3[_3(h)]++;
            ++src;
        }
        else
            *dst++ = *src++;
	}

	unsigned int sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;
	unsigned int tsum;

	for (unsigned int i = 0; i < 256; i += 4)
	{
		tsum = h0[i+0] + sum0; h0[i+0] = sum0 - 1; sum0 = tsum;
		tsum = h0[i+1] + sum0; h0[i+1] = sum0 - 1; sum0 = tsum;
		tsum = h0[i+2] + sum0; h0[i+2] = sum0 - 1; sum0 = tsum;
		tsum = h0[i+3] + sum0; h0[i+3] = sum0 - 1; sum0 = tsum;

		tsum = h1[i+0] + sum1; h1[i+0] = sum1 - 1; sum1 = tsum;
		tsum = h1[i+1] + sum1; h1[i+1] = sum1 - 1; sum1 = tsum;
		tsum = h1[i+2] + sum1; h1[i+2] = sum1 - 1; sum1 = tsum;
		tsum = h1[i+3] + sum1; h1[i+3] = sum1 - 1; sum1 = tsum;

		tsum = h2[i+0] + sum2; h2[i+0] = sum2 - 1; sum2 = tsum;
		tsum = h2[i+1] + sum2; h2[i+1] = sum2 - 1; sum2 = tsum;
		tsum = h2[i+2] + sum2; h2[i+2] = sum2 - 1; sum2 = tsum;
		tsum = h2[i+3] + sum2; h2[i+3] = sum2 - 1; sum2 = tsum;

        if (max_bytes>3)
        {
            tsum = h3[i+0] + sum3; h3[i+0] = sum3 - 1; sum3 = tsum;
            tsum = h3[i+1] + sum3; h3[i+1] = sum3 - 1; sum3 = tsum;
            tsum = h3[i+2] + sum3; h3[i+2] = sum3 - 1; sum3 = tsum;
            tsum = h3[i+3] + sum3; h3[i+3] = sum3 - 1; sum3 = tsum;
        }
	}
#define RADIX_PASS(src, src_end, dst, hist, func) \
	for (auto i = src; i != src_end; ++i) \
	{ \
		unsigned int h = (*i) & 0xffffffff; \
		dst[++hist[func(h)]] = *i; \
	}

    if (max_bytes>3)
    {
        RADIX_PASS(m_triangles.data(), m_triangles.data() + m_triangles.size(), m_rasterizer->m_sort.data(), h0, _0);
        RADIX_PASS(m_rasterizer->m_sort.data(), m_rasterizer->m_sort.data() + m_rasterizer->m_sort.size(), m_triangles.data(), h1, _1);
        RADIX_PASS(m_triangles.data(), m_triangles.data() + m_triangles.size(), m_rasterizer->m_sort.data(), h2, _2);
        RADIX_PASS(m_rasterizer->m_sort.data(), m_rasterizer->m_sort.data() + m_rasterizer->m_sort.size(), m_triangles.data(), h3, _3);
    }
    else
    {
        RADIX_PASS(m_rasterizer->m_sort.data(), m_rasterizer->m_sort.data() + m_rasterizer->m_sort.size(), m_triangles.data(), h0, _0);
        RADIX_PASS(m_triangles.data(), m_triangles.data() + m_triangles.size(), m_rasterizer->m_sort.data(), h1, _1);
        RADIX_PASS(m_rasterizer->m_sort.data(), m_rasterizer->m_sort.data() + m_rasterizer->m_sort.size(), m_triangles.data(), h2, _2);
    }
#endif
}
