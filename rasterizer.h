#pragma once

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "rasterizer_math.h"

#include "rasterizer_tile.h"

#include <bx/uint32_t.h>

#define USE_PACKED_TRIANGLES 1
#define USE_STATS 0

struct ALIGN16 Rasterizer
{
    static constexpr int g_width = 8;
    static constexpr int g_height = 20;
    static constexpr int g_total_width = g_width * Tile::g_tile_width;
    static constexpr int g_total_height = g_height * Tile::g_tile_height;
    static constexpr int g_max_masks_per_tile = g_width * Tile::g_tile_width; // should match max x

#if USE_PACKED_TRIANGLES
    typedef TrianglePacked TriangleType;
    static constexpr int g_fixed_point_bits = 6;
#else
    typedef Triangle TriangleType;
    static constexpr int g_fixed_point_bits = 16;
#endif

    struct SortKey
    {
        uint32_t z;
        uint32_t index;
    };

    struct TrianagleData
    {
        struct TileData
        {
            stl::vector<SortKey>  triangle_indices;
            SortKey*              triangle_index_data = nullptr;
            uint32_t              triangle_index_count = 0;
        };
        stl::vector<TriangleType> triangles;
        TriangleType*             triangle_data = nullptr;
        uint32_t                  triangle_count = 0;
        TileData                  tiles[g_width*g_height];
    };

    struct ThreadData
    {
        TrianagleData        data;
        stl::vector<vec4_t>  vertices;
        stl::vector<SortKey> sort;

        void clear();
    };
private:

    Matrix			        m_transform;
    vec4i_t                 m_full_span;
    vec4_t                  m_tile_size;
    vec4_t                  m_almost_one;
    vec4_t                  m_tile_bounds;

    ThreadData              m_data;
    stl::vector<Tile>       m_tiles;
    stl::vector<ThreadData> m_thread_data;
    stl::vector<uint32_t*>  m_flags;
    stl::vector<vec4i_t>    m_masks;
    uint32_t                m_flag_count = 1;

    bool                    m_mt = false;

    inline void push_4triangles(TrianagleData& data, uint32_t flag, int* bounds_array,
                                const vec4_t* x, const vec4_t* y, const vec4_t* w, bool select_tiles);

    inline void push_triangle_batched(TrianagleData& data, uint32_t flag, const vec4_t* src, int count, const uint16_t* indices,
                                      int* bounds_array, bool select_tiles, bool use_indices);

    bool occlude_object(const vec4_t* m, vec4_t v_min, vec4_t v_max, int* bounds_array);

    void push_object_clipped(ThreadData& data, const uint16_t* indices, int index_count,
                             const vec4_t* vertices, int vertex_count, int* bounds_array, uint32_t flag, bool select_tiles);

    void sort_triangles(SortKey* triangles, uint32_t size, stl::vector<SortKey>& temp);

    __forceinline vec4_t get_tile_bounds( vec4_t v_min, vec4_t v_max );

    __forceinline bool draw_scanlines(Tile& tile, int& xs1, int& xs2, int y1, int y2, int xa1, int xa2, const vec4i_t* masks, uint32_t* flag);

    __forceinline void draw_4triangles(Tile& tile, const TriangleType& tri, uint32_t* flags, const vec4i_t* masks);

    void flush_thread_data(ThreadData& thread_data);
public:

    struct Object
    {
        Matrix transform;
        vec4_t bound_min;
        vec4_t bound_max;
        const uint16_t* indices = nullptr;
        uint32_t index_count = 0;
        const vec4_t* vertices = nullptr;
        uint32_t vertex_count = 0;
        uint32_t* visibility = nullptr;
    };

    bool        m_skip_full = true;
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
    }

    void begin(const Matrix& m);

    void push_box(const Matrix& mat, int* flag);

    void push_objects(const Object* objects, uint32_t object_count, uint32_t thread_index = 0);

    const vec4i_t* get_framebuffer(int tile) const
    {
        return m_tiles[tile].m_frame_buffer;
    }

    const stl::vector<Tile>& get_tiles() const
    {
        return m_tiles;
    }

    void sort_triangles(uint32_t tile, uint32_t thread_index = 0);
    void sort_triangles();

    void draw_triangles(uint32_t tile);
    void draw_triangles();

    void setMT(bool mt)
    {
        m_mt = mt;
    }
    
	void Init(uint32_t num_threads);
};
