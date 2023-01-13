
#include "rasterizer.h"

void Rasterizer::push_box(const Matrix& mat, int* flag)
{
#if 0
    ALIGN16 int bounds_array[4] = { 0 };
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

    vec4_t v_mask00 = VecAnd(VecCmpGt(xxxx0, VecZero()), VecCmpGt(xxxx1, VecZero()));
    vec4_t v_mask01 = VecAnd(VecCmpGt(yyyy0, VecZero()), VecCmpGt(yyyy1, VecZero()));
    vec4_t v_mask10 = VecAnd(VecCmpLt(xxxx0, VecMul(wwww0, g_total_width_v)), VecCmpLt(xxxx1, VecMul(wwww1, g_total_width_v )));
    vec4_t v_mask11 = VecAnd(VecCmpLt(yyyy0, VecMul(wwww0, g_total_height_v)), VecCmpLt(yyyy1, VecMul(wwww1, g_total_height_v)));

    vec4_t v_mask0 = VecAnd(v_mask00, v_mask10);
    vec4_t v_mask1 = VecAnd(v_mask01, v_mask11);
    int mask = VecMask(VecAnd(v_mask0, v_mask1));

    bool intersect_near = VecMask(VecAnd(VecCmpGt(zzzz0, VecZero()), VecCmpGt(zzzz1, VecZero()))) != 15;

    // fast path, no clipping
    if (mask == 15 && intersect_near == false)
    {
        vec4_t xx0 = VecMul(xxxx0, VecRcp(wwww0));
        vec4_t yy0 = VecMul(yyyy0, VecRcp(wwww0));
        vec4_t xx1 = VecMul(xxxx1, VecRcp(wwww1));
        vec4_t yy1 = VecMul(yyyy1, VecRcp(wwww1));

        vec4_t x_min = VecMin(xx0, xx1);
        vec4_t x_max = VecMax(xx0, xx1);
        vec4_t y_min = VecMin(yy0, yy1);
        vec4_t y_max = VecMax(yy0, yy1);

        vec4_t min_0 = VecMin(VecShuffle(x_min, y_min, VecShuffleMask(0, 1, 0, 1)), VecShuffle(x_min, y_min, VecShuffleMask(2, 3, 2, 3)));
        vec4_t max_0 = VecMax(VecShuffle(x_max, y_max, VecShuffleMask(0, 1, 0, 1)), VecShuffle(x_max, y_max, VecShuffleMask(2, 3, 2, 3)));

        vec4_t min_1 = VecMin(VecShuffle(min_0, min_0, VecShuffleMask(0, 2, 0, 0)), VecShuffle(min_0, min_0, VecShuffleMask(1, 3, 0, 0)));
        vec4_t max_1 = VecMax(VecShuffle(max_0, max_0, VecShuffleMask(0, 2, 0, 0)), VecShuffle(max_0, max_0, VecShuffleMask(1, 3, 0, 0)));
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
            vec4_t xtmp0 = VecMoveHL(xx0, xx1); // 2 3 6 7
            vec4_t xtmp1 = VecMoveLH(xx0, xx1); // 0 1 4 5
            vec4_t ytmp0 = VecMoveHL(yy0, yy1); // 2 3 6 7
            vec4_t ytmp1 = VecMoveLH(yy0, yy1); // 0 1 4 5
            vec4_t wtmp0 = VecMoveHL(wwww0, wwww1); // 2 3 6 7
            vec4_t wtmp1 = VecMoveLH(wwww0, wwww1); // 0 1 4 5

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
            vec4_t xtmp0 = VecMoveLH(VecShuffle(xx0, xx0, VecShuffleMask(1, 2, 1, 2)), VecShuffle(xx1, xx1, VecShuffleMask(1, 2, 1, 2))); // 1 2 5 6
            vec4_t xtmp1 = VecMoveLH(VecShuffle(xx0, xx0, VecShuffleMask(0, 3, 0, 3)), VecShuffle(xx1, xx1, VecShuffleMask(3, 0, 3, 0))); // 0 3 7 4
            vec4_t ytmp0 = VecMoveLH(VecShuffle(yy0, yy0, VecShuffleMask(1, 2, 1, 2)), VecShuffle(yy1, yy1, VecShuffleMask(1, 2, 1, 2))); // 1 2 5 6
            vec4_t ytmp1 = VecMoveLH(VecShuffle(yy0, yy0, VecShuffleMask(0, 3, 0, 3)), VecShuffle(yy1, yy1, VecShuffleMask(3, 0, 3, 0))); // 0 3 7 4
            vec4_t wtmp0 = VecMoveLH(VecShuffle(wwww0, wwww0, VecShuffleMask(1, 2, 1, 2)), VecShuffle(wwww1, wwww1, VecShuffleMask(1, 2, 1, 2))); // 1 2 5 6
            vec4_t wtmp1 = VecMoveLH(VecShuffle(wwww0, wwww0, VecShuffleMask(0, 3, 0, 3)), VecShuffle(wwww1, wwww1, VecShuffleMask(3, 0, 3, 0))); // 0 3 7 4

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
#endif
}
