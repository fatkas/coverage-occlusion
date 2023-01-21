
#include <assert.h>
#include "rasterizer_tile.h"

Tile::Tile(int x, int y)
    : m_x(x)
{
    m_y = Vector4(y*g_tile_height);
    clear();
}
