
#include "rasterizer.h"
#include <algorithm>

// 11bit radix http://stereopsis.com/radix.html
#define USE_STD_SORT 0
#define USE_11BIT_RADIX 1

void Rasterizer::sort_triangles(SortKey* triangles, uint32_t size, stl::vector<SortKey>& temp)
{
#if USE_STD_SORT
    std::sort(triangles, triangles + size, [](const SortKey& a, const SortKey& b)
    {
        return a.z < b.z;
    });
#elif USE_11BIT_RADIX
    constexpr uint32_t histogramSize = 2048;

    unsigned int histograms[histogramSize * 2];
    memset(histograms, 0, sizeof(histograms));

    unsigned int* h0 = histograms;
    unsigned int* h1 = histograms + histogramSize;

    assert(temp.size() >= size);

#define _0(h) ((h) & 0x7FF)
#define _1(h) (((h) >> 11) & 0x7FF)

    SortKey* src = triangles, *src_end = triangles + size;
    while (src < src_end)
    {
        unsigned int h = src->z;
        h0[_0(h)]++;
        h1[_1(h)]++;
        ++src;
    }

    unsigned int sum0 = 0, sum1 = 0;
    unsigned int tsum;

    for (unsigned int i = 0; i < histogramSize; i += 4)
    {
        tsum = h0[i+0] + sum0; h0[i+0] = sum0 - 1; sum0 = tsum;
        tsum = h0[i+1] + sum0; h0[i+1] = sum0 - 1; sum0 = tsum;
        tsum = h0[i+2] + sum0; h0[i+2] = sum0 - 1; sum0 = tsum;
        tsum = h0[i+3] + sum0; h0[i+3] = sum0 - 1; sum0 = tsum;

        tsum = h1[i+0] + sum1; h1[i+0] = sum1 - 1; sum1 = tsum;
        tsum = h1[i+1] + sum1; h1[i+1] = sum1 - 1; sum1 = tsum;
        tsum = h1[i+2] + sum1; h1[i+2] = sum1 - 1; sum1 = tsum;
        tsum = h1[i+3] + sum1; h1[i+3] = sum1 - 1; sum1 = tsum;
    }

#define RADIX_PASS(src, src_end, dst, hist, func) \
    for (auto i = src; i != src_end; ++i) \
    { \
        unsigned int h = i->z; \
        dst[++hist[func(h)]] = *i; \
    }

    RADIX_PASS(triangles, triangles + size, temp.data(), h0, _0);
    RADIX_PASS(temp.data(), temp.data() + size, triangles, h1, _1);

#undef RADIX_PASS
#undef _0
#undef _1

#else
    constexpr unsigned max_bytes = 3;

    unsigned int histograms[256*max_bytes];

    memset(histograms, 0, sizeof(histograms));

    unsigned int* h0 = histograms;
    unsigned int* h1 = histograms + 256;
    unsigned int* h2 = histograms + 256*2;
    unsigned int* h3 = histograms + 256*3;

    assert(temp.size() >= size);

    SortKey* src = triangles, *src_end = triangles + size, *dst = temp.data();
    while (src < src_end)
    {
#define _0(h) ((h) & 255)
#define _1(h) (((h) >> 8) & 255)
#define _2(h) (((h) >> 16) & 255)
#define _3(h) ((h) >> 24)

        unsigned int h = src->z;
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
        unsigned int h = i->z; \
        dst[++hist[func(h)]] = *i; \
    }

    if (max_bytes>3)
    {
        RADIX_PASS(triangles, triangles + size, temp.data(), h0, _0);
        RADIX_PASS(temp.data(), temp.data() + size, triangles, h1, _1);
        RADIX_PASS(triangles, triangles + size, temp.data(), h2, _2);
        RADIX_PASS(temp.data(), temp.data() + size, triangles, h3, _3);
    }
    else
    {
        RADIX_PASS(temp.data(), temp.data() + size, triangles, h0, _0);
        RADIX_PASS(triangles, triangles + size, temp.data(), h1, _1);
        RADIX_PASS(temp.data(), temp.data() + size, triangles, h2, _2);
    }
#endif
}

void Rasterizer::sort_triangles(uint32_t tile_index, uint32_t thread_index)
{
    auto & tile = m_data.data.tiles[tile_index];
    if (tile.triangle_index_count)
        sort_triangles(tile.triangle_index_data, tile.triangle_index_count, m_thread_data[thread_index].data.tiles[tile_index].triangle_indices);
}

void Rasterizer::sort_triangles()
{
    for (uint32_t i = 0; i < g_width*g_height; ++i)
        sort_triangles(i, 0);
}
