#pragma once

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "rasterizer_math.h"

#include "rasterizer_tile.h"

#include <bx/uint32_t.h>

#define USE_PACKED_TRIANGLES 1
#define USE_STATS 0
#define USE_NORMAL_MASKS 1

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
        uint64_t z : 22;
        uint64_t flag : 21; // up to 512k objects pushed
        uint64_t index : 21; // up to 512k triangles
    };
    static_assert(sizeof(SortKey)==8, "the sort key should be 8 bytes");

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
        // this is where we put triangles/tile indices in worker/main threads
        TrianagleData           data;
        // temp thread local data that holds tranformed/clipped vertices
        stl::vector<vec4_t>     vertices;
        stl::vector<vec2_t>     positions;
        stl::vector<uint32_t>   depths;
        // temp thead local data that holds clipped indices
        stl::vector<uint16_t>   indices;
        // temp thread local data that's used for radix sort
        stl::vector<SortKey>    sort;

        void clear();
    };
private:

    Matrix			        m_transform;
    vec4i_t                 m_full_span;
    vec4_t                  m_tile_size;
    vec4_t                  m_almost_one;
    vec4_t                  m_tile_bounds;
    vec4_t                  m_tile_height_v;
    vec4_t                  m_camera_position;
    vec4_t                  m_fixed_point;
    vec4_t                  m_inv_fixed_point;
    vec4_t                  m_total_size;

    ThreadData              m_data;
    stl::vector<Tile>       m_tiles;
    stl::vector<ThreadData> m_thread_data;
    stl::vector<uint32_t*>  m_flags;
    stl::vector<vec4i_t>    m_masks;
    uint32_t                m_flag_count = 1;

    bool                    m_mt = false;

    inline void push_4triangles(TrianagleData& data, uint32_t flag, int* bounds_array,
                                vec2_t v0[4], vec2_t v1[4], vec2_t v2[4], uint32_t group_w);

    inline void push_triangle_batched(TrianagleData& data, uint32_t flag, const vec2_t* src, const uint32_t* w, const uint16_t* indices, int count,
                                      const uint8_t* normal_masks, uint32_t normal_mask, int* bounds_array, bool select_tiles);

    bool occlude_object(const vec4_t* m, vec4_t v_min, vec4_t v_max, int* bounds_array);

    void push_object_clipped(ThreadData& data, const uint16_t* indices, int index_count,
                             vec4_t* transformed_vertices, uint32_t vertex_count,
                             int* bounds_array, uint32_t flag, bool select_tiles);

    void sort_triangles(SortKey* triangles, uint32_t size, stl::vector<SortKey>& temp);

    __forceinline vec4_t get_tile_bounds(vec4_t minmax);

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

        // https://cgvr.cs.uni-bremen.de/teaching/cg_literatur/backface_normal_masks.pdf
        // instead of calculating triangle orientation lets use triangle masks
        // since the normal mask is uint8, store up to 8 normals
        const uint8_t* normal_masks = nullptr;
        vec4_t* normals = nullptr;
        uint8_t normal_count = 0;
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
    uint32_t    m_triangles_backface = 0;

    Rasterizer()
    {
    }

    void begin(const Matrix& m, vec4_t cam);

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

    const TrianagleData::TileData* get_tile_data() const
    {
        return m_data.data.tiles;
    }

    void sort_triangles(uint32_t tile, uint32_t thread_index = 0);
    void sort_triangles();

    void draw_triangles(uint32_t tile);
    void draw_triangles();

    uint32_t get_total_groups() const
    {
        return m_data.data.triangle_count;
    }

    void setMT(bool mt)
    {
        m_mt = mt;
    }
    
	void Init(uint32_t num_threads);
};
