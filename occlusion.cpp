/*
 * Copyright 2011-2022 Branimir Karadzic. All rights reserved.
 * License: https://github.com/bkaradzic/bgfx/blob/master/LICENSE
 */

#include "common.h"
#include "bgfx_utils.h"

#include <bx/uint32_t.h>
#include <bx/thread.h>
#include <bx/os.h>
#include "imgui/imgui.h"
#include "camera.h"

#include "rasterizer.h"

#include "TaskScheduler.h"

#include <bgfx/embedded_shader.h>

// embedded shaders
#include "vs_occlusion.bin.h"
#include "fs_occlusion.bin.h"

#include "vs_fullscreen.bin.h"
#include "fs_fullscreen.bin.h"

#include <string>
#include <thread>

namespace
{

static const bgfx::EmbeddedShader s_embeddedShaders[] =
{
	BGFX_EMBEDDED_SHADER(vs_occlusion),
	BGFX_EMBEDDED_SHADER(fs_occlusion),

    BGFX_EMBEDDED_SHADER(vs_fullscreen),
    BGFX_EMBEDDED_SHADER(fs_fullscreen),

	BGFX_EMBEDDED_SHADER_END()
};

struct PosColorVertex
{
	float m_x;
	float m_y;
	float m_z;
	uint32_t m_abgr;

	static void init()
	{
		ms_layout
			.begin()
			.add(bgfx::Attrib::Position, 3, bgfx::AttribType::Float)
			.add(bgfx::Attrib::Color0,   4, bgfx::AttribType::Uint8, true)
			.end();
	}

	static bgfx::VertexLayout ms_layout;
};

struct PosVertex
{
    float m_x;
    float m_y;
    float m_z;

    static void init()
    {
        ms_layout
            .begin()
            .add(bgfx::Attrib::Position,  3, bgfx::AttribType::Float)
            .end();
    }

    static bgfx::VertexLayout ms_layout;
};

bgfx::VertexLayout PosVertex::ms_layout;
bgfx::VertexLayout PosColorVertex::ms_layout;

static PosColorVertex s_cubeVertices[8] =
{
	{-1.0f,  1.0f,  1.0f, 0xff000000 },
	{ 1.0f,  1.0f,  1.0f, 0xff0000ff },
	{-1.0f, -1.0f,  1.0f, 0xff00ff00 },
	{ 1.0f, -1.0f,  1.0f, 0xff00ffff },
	{-1.0f,  1.0f, -1.0f, 0xffff0000 },
	{ 1.0f,  1.0f, -1.0f, 0xffff00ff },
	{-1.0f, -1.0f, -1.0f, 0xffffff00 },
	{ 1.0f, -1.0f, -1.0f, 0xffffffff },
};

static vec4_t s_cubeVerticesSIMD[8] =
{
    {-1.0f,  1.0f,  1.0f, 1.f },
    { 1.0f,  1.0f,  1.0f, 1.f },
    {-1.0f, -1.0f,  1.0f, 1.f },
    { 1.0f, -1.0f,  1.0f, 1.f },
    {-1.0f,  1.0f, -1.0f, 1.f },
    { 1.0f,  1.0f, -1.0f, 1.f },
    {-1.0f, -1.0f, -1.0f, 1.f },
    { 1.0f, -1.0f, -1.0f, 1.f },
};

static const uint16_t s_cubeIndices[36] =
{
	0, 1, 2, // 0
	1, 3, 2,
	4, 6, 5, // 2
	5, 6, 7,
	0, 2, 4, // 4
	4, 2, 6,
	1, 5, 3, // 6
	5, 7, 3,
	0, 4, 1, // 8
	4, 5, 1,
	2, 3, 6, // 10
	6, 3, 7,
};

static const float s_mod[6][3] =
{
	{ 1.0f, 1.0f, 1.0f },
	{ 1.0f, 0.0f, 0.0f },
	{ 0.0f, 1.0f, 0.0f },
	{ 0.0f, 0.0f, 1.0f },
	{ 1.0f, 1.0f, 0.0f },
	{ 0.0f, 1.0f, 1.0f },
};

static void screenSpaceQuad(bgfx::Encoder* encoder)
{
    if (3 == bgfx::getAvailTransientVertexBuffer(3, PosVertex::ms_layout) )
    {
        bgfx::TransientVertexBuffer vb;
        bgfx::allocTransientVertexBuffer(&vb, 3, PosVertex::ms_layout);
        PosVertex* vertex = (PosVertex*)vb.data;

        const float zz = 0.0f;

        vertex[0].m_x = -1.f;
        vertex[0].m_y = -1.f;
        vertex[0].m_z = zz;

        vertex[1].m_x = 3.f;
        vertex[1].m_y = -1.f;
        vertex[1].m_z = zz;

        vertex[2].m_x = -1.f;
        vertex[2].m_y = 3.f;
        vertex[2].m_z = zz;

        encoder->setVertexBuffer(0, &vb);
    }
}

class ExampleOcclusionCulling : public entry::AppI
{
public:
    ExampleOcclusionCulling(const char* _name, const char* _description, const char* _url)
		: entry::AppI(_name, _description, _url)
        , m_DrawTasks(this)
        , m_SortTasks(this)
        , m_PushTasks(this)
	{
	}

	void init(int32_t _argc, const char* const* _argv, uint32_t _width, uint32_t _height) override
	{
		Args args(_argc, _argv);

		m_width  = _width;
		m_height = _height;
		m_debug  = BGFX_DEBUG_NONE;
		m_reset  = BGFX_RESET_NONE;

		m_scrollArea = 0;
		m_dim        = 1;
		m_maxDim     = 40;
		m_transform  = 1;

		m_timeOffset = bx::getHPCounter();
        m_timeStop = m_timeOffset;

		m_deltaTimeNs    = 0;
		m_deltaTimeAvgNs = 0;
		m_numFrames      = 0;

		bgfx::Init init;
		init.type     = args.m_type;
		init.vendorId = args.m_pciId;
		init.platformData.nwh  = entry::getNativeWindowHandle(entry::kDefaultWindowHandle);
		init.platformData.ndt  = entry::getNativeDisplayHandle();
		init.resolution.width  = m_width;
		init.resolution.height = m_height;
		init.resolution.reset  = m_reset;
		bgfx::init(init);

		const bgfx::Caps* caps = bgfx::getCaps();
		m_maxDim = (int32_t)bx::pow(float(caps->limits.maxDrawCalls), 1.0f/3.0f);

		// Enable debug text.
		bgfx::setDebug(m_debug);

		// Set view 0 clear state.
		bgfx::setViewClear(0
			, BGFX_CLEAR_COLOR|BGFX_CLEAR_DEPTH
			, 0x303030ff
			, 1.0f
			, 0
			);

		// Create vertex stream declaration.
		PosColorVertex::init();
        PosVertex::init();

		bgfx::RendererType::Enum type = bgfx::getRendererType();

		// Create program from shaders.
		m_program = bgfx::createProgram(
			  bgfx::createEmbeddedShader(s_embeddedShaders, type, "vs_occlusion")
			, bgfx::createEmbeddedShader(s_embeddedShaders, type, "fs_occlusion")
			, true /* destroy shaders when program is destroyed */
			);

        m_fullscreen = bgfx::createProgram(
              bgfx::createEmbeddedShader(s_embeddedShaders, type, "vs_fullscreen")
            , bgfx::createEmbeddedShader(s_embeddedShaders, type, "fs_fullscreen")
            , true /* destroy shaders when program is destroyed */
            );
        s_texColor = bgfx::createUniform("s_texColor", bgfx::UniformType::Sampler);

		// Create static vertex buffer.
		m_vbh = bgfx::createVertexBuffer(
			  bgfx::makeRef(s_cubeVertices, sizeof(s_cubeVertices) )
			, PosColorVertex::ms_layout
			);

		// Create static index buffer.
		m_ibh = bgfx::createIndexBuffer(bgfx::makeRef(s_cubeIndices, sizeof(s_cubeIndices) ) );

		// Imgui.
		imguiCreate();

        cameraCreate();
        cameraSetPosition({ 0.0f, 0.0f, -35.0f });
        cameraSetVerticalAngle(0.f);

        debug_coverage = bgfx::createTexture2D(Rasterizer::g_total_width, Rasterizer::g_total_height, false, 1, bgfx::TextureFormat::R8);

        m_Scheduler.Initialize(std::thread::hardware_concurrency()-1);
        int workersCount = (int)m_Scheduler.GetNumTaskThreads();
        printf("Scheduler started, %d workers\n", workersCount);

        m_Rasterizer.Init(workersCount);
	}

	int shutdown() override
	{
		// Cleanup.
        cameraDestroy();
		imguiDestroy();
		bgfx::destroy(m_ibh);
		bgfx::destroy(m_vbh);
		bgfx::destroy(m_program);
        bgfx::destroy(m_fullscreen);
        bgfx::destroy(debug_coverage);
        bgfx::destroy(s_texColor);

		// Shutdown bgfx.
		bgfx::shutdown();

		return 0;
	}

	void submit()
	{
		bgfx::Encoder* encoder = bgfx::begin();

		if (NULL != encoder)
		{
			for (uint32_t i = 0; i < uint32_t(m_dim)*uint32_t(m_dim)*uint32_t(m_dim); ++i)
			{
                if (m_Visibility[i] == 0)
                    continue;

				encoder->setTransform((float*)&m_Objects[i].transform);
                encoder->setVertexBuffer(0, m_vbh);
                encoder->setIndexBuffer(m_ibh);

                encoder->setState(
                                  BGFX_STATE_DEFAULT
                                  | (m_Wireframe ? BGFX_STATE_PT_LINES | BGFX_STATE_LINEAA | BGFX_STATE_BLEND_ALPHA : 0)
                                  );
                encoder->submit(0, m_program);
			}

            if (m_ShowCoverage && m_Occlusion)
            {
                encoder->setTexture(0, s_texColor, debug_coverage);
                encoder->setState(BGFX_STATE_WRITE_RGB | BGFX_STATE_WRITE_A);
                screenSpaceQuad(encoder);
                encoder->submit(0, m_fullscreen);
            }

			bgfx::end(encoder);
		}
	}

    struct PushTask : enki::ITaskSet
    {
        PushTask(ExampleOcclusionCulling* parent)
            : m_parent(parent)
        {
        }

        void setCount(uint32_t count)
        {
            this->m_SetSize = count;
            this->m_MinRange = 1024;
        }

        void ExecuteRange(enki::TaskSetPartition range, uint32_t thread_index) override
        {
            m_parent->m_Rasterizer.push_objects(m_parent->m_Objects.data() + range.start, range.end - range.start, thread_index);
        }

        ExampleOcclusionCulling* m_parent = nullptr;
    };

    struct SortTask : enki::ITaskSet
    {
        SortTask(ExampleOcclusionCulling* parent)
            : m_parent(parent)
        {
            this->m_SetSize = Rasterizer::g_width*Rasterizer::g_height;
        }

        void ExecuteRange(enki::TaskSetPartition range, uint32_t thread_index) override
        {
            for (uint32_t index = range.start; index < range.end; ++index)
            {
                m_parent->m_Rasterizer.sort_triangles(index, thread_index);
            }
        }

        ExampleOcclusionCulling* m_parent = nullptr;
    };

    struct DrawTask : enki::ITaskSet
    {
        DrawTask(ExampleOcclusionCulling* parent)
            : m_parent(parent)
        {
            this->m_SetSize = Rasterizer::g_width*Rasterizer::g_height;
        }

        void ExecuteRange(enki::TaskSetPartition range, uint32_t) override
        {
            for (uint32_t index = range.start; index < range.end; ++index)
            {
                m_parent->m_Rasterizer.draw_triangles(index);
            }
        }

        ExampleOcclusionCulling* m_parent = nullptr;
    };

	bool update() override
	{
		if (!entry::processEvents(m_width, m_height, m_debug, m_reset, &m_mouseState) )
		{
			int64_t now = bx::getHPCounter();
			static int64_t last = now;
			const int64_t hpFreq = bx::getHPFrequency();
			const int64_t frameTime = now - last;
			last = now;
			const double freq = double(hpFreq);
			const double toMs = 1000.0/freq;

			m_deltaTimeNs += frameTime*1000000/hpFreq;

            int old_transform = m_transform;

			imguiBeginFrame(m_mouseState.m_mx
				,  m_mouseState.m_my
				, (m_mouseState.m_buttons[entry::MouseButton::Left  ] ? IMGUI_MBUT_LEFT   : 0)
				| (m_mouseState.m_buttons[entry::MouseButton::Right ] ? IMGUI_MBUT_RIGHT  : 0)
				| (m_mouseState.m_buttons[entry::MouseButton::Middle] ? IMGUI_MBUT_MIDDLE : 0)
				,  m_mouseState.m_mz
				, uint16_t(m_width)
				, uint16_t(m_height)
				);

			showExampleDialog(this);

			ImGui::SetNextWindowPos(
				  ImVec2((float)m_width - (float)m_width / 4.0f - 10.0f, 10.0f)
				, ImGuiCond_FirstUseEver
				);
			ImGui::SetNextWindowSize(
				  ImVec2((float)m_width / 4.0f, m_height - 20.f)
				, ImGuiCond_FirstUseEver
				);
			ImGui::Begin("Settings"
				, NULL
				, 0
				);

            uint32_t visible = 0;
            for (auto & v : m_Visibility)
                visible += v != 0;

			ImGui::RadioButton("Rotate",&m_transform,0);
			ImGui::RadioButton("No rotate",&m_transform,1);
			ImGui::Separator();

            ImGui::Checkbox("Wireframe", &m_Wireframe);
            ImGui::Checkbox("Enable occlusion", &m_Occlusion);
            ImGui::Checkbox("Show coverage buffer", &m_ShowCoverage);
            ImGui::Checkbox("Skip full tiles", &m_Rasterizer.m_skip_full);
            ImGui::Checkbox("Use fast path for box", &m_UseBox);
            ImGui::Checkbox("Enable MT submit/draw/push", &m_MT);
			ImGui::SliderInt("Dim", &m_dim, 1, m_maxDim);
            ImGui::SliderFloat("Spacing", &m_Spacing, 50.f, 100.f);
			ImGui::Text("Draw calls: %d", m_dim*m_dim*m_dim);
            ImGui::Text("Draw calls visible: %d", visible);

			ImGui::Separator();
            if (ImGui::TreeNode("Occlusion time", "Occlusion time %f", double(occlusion_push_time + occlusion_sort_time + occlusion_draw_time)*toMs))
            {
                ImGui::Text("occlusion push %0.6f [ms]", double(occlusion_push_time)*toMs);
                ImGui::Text("occlusion sort %0.6f [ms]", double(occlusion_sort_time)*toMs);
                ImGui::Text("occlusion draw %0.6f [ms]", double(occlusion_draw_time)*toMs);
                ImGui::TreePop();
            }

            ImGui::Text("total triangles %d", m_Rasterizer.m_triangles_total);
            ImGui::Text("total occluder triangles %d", m_Rasterizer.m_triangles_occluder_total);
            ImGui::Text("total occludee triangles %d", m_Rasterizer.m_triangles_occludee_total);
            ImGui::Text("total drawn triangles %d", m_Rasterizer.m_triangles_drawn_total);
            ImGui::Text("total drawn occluder triangles %d", m_Rasterizer.m_triangles_drawn_occluder_total);
            ImGui::Text("total drawn occludee triangles %d", m_Rasterizer.m_triangles_drawn_occludee_total);
            ImGui::Text("total skipped triangles %d", m_Rasterizer.m_triangles_skipped);
            ImGui::Text("total offscreen triangles %d", m_Rasterizer.m_triangles_offscreen);

            const auto& tiles = m_Rasterizer.get_tiles();
            if (ImGui::TreeNode("Tiles", "Tiles (%d)", (uint32_t)tiles.size()))
            {
                for (auto & t : tiles)
                {
                    if (ImGui::TreeNode(std::to_string(t.m_x + t.m_y*Rasterizer::g_width).c_str(), "Tile %d (%d/%d) %s", t.m_x, /*(uint32_t)t.m_triangle_count*/0, t.m_triangles_drawn_total, t.m_mask == ~0u ? "full" : ""))
                    {
                        ImGui::Text("total sorted triangles %d", /*(uint32_t)t.m_triangles.size()*/0);
                        ImGui::Text("total drawn triangles %d", t.m_triangles_drawn_total);
                        ImGui::Text("total drawn occluder triangles %d", t.m_triangles_drawn_occluder_total);
                        ImGui::Text("total drawn occludee triangles %d", t.m_triangles_drawn_occludee_total);
                        ImGui::Text("total skipped triangles %d", t.m_triangles_skipped);
                        ImGui::Text("mask %llu", t.m_mask);
                        ImGui::TreePop();
                    }
                }
                ImGui::TreePop();
            }

			ImGui::End();

			imguiEndFrame();

            if (old_transform == 0 && m_transform == 1)
                m_timeStop = now;
            if (m_transform == 1)
                now = m_timeStop;

            cameraUpdate(m_deltaTimeNs/1000000.f, m_mouseState, ImGui::MouseOverArea());

			// const bx::Vec3 at  = { 0.0f, 0.0f,   0.0f };
			// const bx::Vec3 eye = { 0.0f, 0.0f, -35.0f };

			float view[16];
            cameraGetViewMtx(view);
			// bx::mtxLookAt(view, eye, at);

			const bgfx::Caps* caps = bgfx::getCaps();
			float proj[16];
			bx::mtxProj(proj, 60.0f, float(m_width)/float(m_height), 0.1f, 100.0f, caps->homogeneousDepth);

			// Set view and projection matrix for view 0.
			bgfx::setViewTransform(0, view, proj);

			// Set view 0 default viewport.
			bgfx::setViewRect(0, 0, 0, uint16_t(m_width), uint16_t(m_height) );

			// This dummy draw call is here to make sure that view 0 is cleared
			// if no other draw calls are submitted to view 0.
			bgfx::touch(0);

            {
                float time = (float)( (now-m_timeOffset)/freq);

                const float* mod = s_mod[0];

                float mtxS[16];
                const float scale = 0.25f;
                bx::mtxScale(mtxS, scale, scale, scale);

                const float step = m_Spacing * 0.01f;
                float pos[3];
                pos[0] = -step*m_dim / 2.0f;
                pos[1] = -step*m_dim / 2.0f;
                pos[2] = -15.0;

                uint32_t max_drawcalls = m_dim*m_dim*m_dim;
                m_Objects.resize(max_drawcalls*2);
                m_Visibility.resize(max_drawcalls);
                for (uint32_t zz = 0; zz < uint32_t(m_dim); ++zz)
                {
                    for (uint32_t yy = 0; yy < uint32_t(m_dim); ++yy)
                    {
                        for (uint32_t xx = 0; xx < uint32_t(m_dim); ++xx)
                        {
                            float mtxR[16];
                            bx::mtxRotateXYZ(mtxR
                                , (time + xx*0.21f)*mod[0]
                                , (time + yy*0.37f)*mod[1]
                                , (time + zz*0.13f)*mod[2]
                                );

                            float mtx[16];
                            bx::mtxMul(mtx, mtxS, mtxR);

                            mtx[12] = pos[0] + float(xx)*step;
                            mtx[13] = pos[1] + float(yy)*step;
                            mtx[14] = pos[2] + float(zz)*step;

                            uint32_t idx = xx + yy*m_dim + zz*m_dim*m_dim;
                            m_Objects[idx].transform = MatrixSet(mtx);
                            m_Objects[idx].indices = s_cubeIndices;
                            m_Objects[idx].index_count = sizeof(s_cubeIndices) / sizeof(s_cubeIndices[0]);
                            m_Objects[idx].vertices = s_cubeVerticesSIMD;
                            m_Objects[idx].vertex_count = sizeof(s_cubeVerticesSIMD) / sizeof(s_cubeVerticesSIMD[0]);
                            m_Objects[idx].visibility = &m_Visibility[xx + yy*m_dim + zz*m_dim*m_dim];
                            m_Objects[idx].bound_min = {-1.f, -1.f, -1.f, 1.f};
                            m_Objects[idx].bound_max = {1.f, 1.f, 1.f, 1.f};

                            m_Objects[idx+max_drawcalls].transform = MatrixSet(mtx);
                            m_Objects[idx+max_drawcalls].indices = s_cubeIndices;
                            m_Objects[idx+max_drawcalls].index_count = sizeof(s_cubeIndices) / sizeof(s_cubeIndices[0]);
                            m_Objects[idx+max_drawcalls].vertices = s_cubeVerticesSIMD;
                            m_Objects[idx+max_drawcalls].vertex_count = sizeof(s_cubeVerticesSIMD) / sizeof(s_cubeVerticesSIMD[0]);
                            m_Objects[idx+max_drawcalls].visibility = nullptr;
                            m_Objects[idx+max_drawcalls].bound_min = {-1.f, -1.f, -1.f, 1.f};
                            m_Objects[idx+max_drawcalls].bound_max = {1.f, 1.f, 1.f, 1.f};
                        }
                    }
                }
            }

            Matrix view_mat = MatrixSet(view), proj_mat = MatrixSet(proj);
            m_Rasterizer.begin(view_mat * proj_mat * MatrixScaling(0.5f, -0.5f, 1.0f) * MatrixTranslation(Vector4( .5f, 0.5f, 0.0f, 1.0f )) * MatrixScaling( (float)Rasterizer::g_total_width, (float)Rasterizer::g_total_height, 1.0f));
            if (m_Occlusion)
            {
                m_Rasterizer.setMT(m_MT);

                int64_t occlusion_start = bx::getHPCounter();
                if (m_MT)
                {
                    m_PushTasks.setCount(m_Objects.size());
                    m_Scheduler.AddTaskSetToPipe(&m_PushTasks);
                    m_Scheduler.WaitforAll();
                }
                else
                {
                    m_Rasterizer.push_objects(m_Objects.data(), m_Objects.size());
                }
                int64_t occlusion_mid = bx::getHPCounter();
                occlusion_push_time = occlusion_mid - occlusion_start;

                {
                    occlusion_start = bx::getHPCounter();
                    if (m_MT)
                    {
                        m_Scheduler.AddTaskSetToPipe(&m_SortTasks);
                        m_Scheduler.WaitforAll();
                    }
                    else
                        m_Rasterizer.sort_triangles();
                    occlusion_mid = bx::getHPCounter();
                    if (m_MT)
                    {
                        m_Scheduler.AddTaskSetToPipe(&m_DrawTasks);
                        m_Scheduler.WaitforAll();
                    }
                    else
                        m_Rasterizer.draw_triangles();
                    int64_t occlusion_end = bx::getHPCounter();

                    occlusion_sort_time = occlusion_mid - occlusion_start;
                    occlusion_draw_time = occlusion_end - occlusion_mid;
                }

                if (m_ShowCoverage)
                {
                    stl::vector<uint8_t> data(Rasterizer::g_total_width*Rasterizer::g_total_height);
                    for ( int i = 0; i < Rasterizer::g_height; ++i )
                    {
                        for ( int j = 0; j < Rasterizer::g_width; ++j )
                        {
                            for ( int y = 0; y < Tile::g_tile_height; ++y )
                                for ( int x = 0; x < Tile::g_tile_width; ++x )
                                {
                                    vec4i_t buf = m_Rasterizer.get_framebuffer(j + i*Rasterizer::g_width)[y];
                                    unsigned int mask = ( (unsigned int*)( &buf ) )[ x >> 5 ];
                                    int bit = mask & ( 1 << ( x & 31 ) );
                                    data[j*Tile::g_tile_width + 127 - x + (y + i*Tile::g_tile_height)*Rasterizer::g_total_width] = bit ? 255 : 0;
                                }
                        }
                    }
                    const bgfx::Memory* mem = bgfx::makeRef(data.data(), (uint32_t)data.size());
                    bgfx::updateTexture2D(debug_coverage, 0, 0, 0, 0, Rasterizer::g_total_width, Rasterizer::g_total_height, mem, Rasterizer::g_total_width);
                }
            }
            else
            {
                occlusion_push_time = occlusion_sort_time = occlusion_draw_time = 0;
                for (uint32_t i = 0; i < uint32_t(m_dim)*uint32_t(m_dim)*uint32_t(m_dim); ++i )
                    m_Visibility[i] = 1;
            }

			submit();

			// Advance to next frame. Rendering thread will be kicked to
			// process submitted rendering primitives.
			bgfx::frame();

			return true;
		}

		return false;
	}

	entry::MouseState m_mouseState;

	uint32_t m_width;
	uint32_t m_height;
	uint32_t m_debug;
	uint32_t m_reset;

    float    m_Spacing = 60.f;

    stl::vector<Rasterizer::Object> m_Objects;
    stl::vector<uint32_t> m_Visibility;

    bool     m_Wireframe = false;
    bool     m_Occlusion = true;
    bool     m_ShowCoverage = false;
    bool     m_UseBox = false;
    bool     m_MT = false;

    int64_t  occlusion_push_time = 0;
    int64_t  occlusion_sort_time = 0;
    int64_t  occlusion_draw_time = 0;

	int32_t  m_scrollArea;
	int32_t  m_dim;
	int32_t  m_maxDim;
	int32_t  m_transform;

	int64_t  m_timeOffset;
    int64_t  m_timeStop = 0;

	int64_t  m_deltaTimeNs;
	int64_t  m_deltaTimeAvgNs;
	int64_t  m_numFrames;

    Rasterizer m_Rasterizer;

	bgfx::ProgramHandle m_program;
    bgfx::ProgramHandle m_fullscreen;
	bgfx::VertexBufferHandle m_vbh;
	bgfx::IndexBufferHandle  m_ibh;
    bgfx::UniformHandle s_texColor;

    enki::TaskScheduler m_Scheduler;
    DrawTask m_DrawTasks;
    SortTask m_SortTasks;
    PushTask m_PushTasks;

    bgfx::TextureHandle debug_coverage;
};

} // namespace

ENTRY_IMPLEMENT_MAIN(
	  ExampleOcclusionCulling
	, "occlusion-cullsing"
	, "Occlusion culling demo"
	, "no link for now"
	);
