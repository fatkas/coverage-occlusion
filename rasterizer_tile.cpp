
#include <assert.h>
#include "rasterizer_tile.h"

inline static int Shift( int val )
{
    if( val > 31 )
        return 0;
    if( val < 0 )
        return 0xffffffff;
    return 0xffffffff >> val;
}

Tile::Tile(int x, int y)
    : m_x(x)
    , m_y(y)
{
    for ( int i = 0; i < g_tile_height; ++i )
        m_frame_buffer[i] = VecIntZero();

    m_shifts.resize(4096);
    for (int i = 0; i < x*128; ++i )
        m_shifts[i] = Vector4Int( 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff );
    for (int i = (x+1)*128; i < 4096; ++i )
        m_shifts[i] = VecIntZero();
    for (int i = 0; i < 128; ++i )
        m_shifts[x*128 + i] = Vector4Int(Shift(i - 96), Shift(i - 64), Shift(i - 32), Shift(i));
}
