set(BGFX_ROOT_DIR ${PROJECT_SOURCE_DIR})

SET(BIMG_SOURCES
  ${BGFX_ROOT_DIR}/bimg/src/image.cpp
  ${BGFX_ROOT_DIR}/bimg/src/image_gnf.cpp
  ${BGFX_ROOT_DIR}/bimg/src/image_decode.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_weight_quant_xfer_tables.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_quantization.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_averages_and_directions.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_ideal_endpoints_and_weights.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_pick_best_endpoint_format.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_find_best_partitioning.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_partition_tables.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_compress_symbolic.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_mathlib.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_weight_align.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_block_sizes.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_compute_variance.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_percentile_tables.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_integer_sequence.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_symbolic_physical.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_color_quantize.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_color_unquantize.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_platform_isa_detection.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_image.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_compress_symbolic.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_decompress_symbolic.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_mathlib_softfloat.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/source/astcenc_entry.cpp
  ${BGFX_ROOT_DIR}/bimg/3rdparty/tinyexr/deps/miniz/miniz.c
)

set(BX_SOURCES
  ${BGFX_ROOT_DIR}/bx/src/amalgamated.cpp
)

set(BGFX_SOURCES
  ${BGFX_ROOT_DIR}/bgfx/src/amalgamated.mm

  ${BGFX_ROOT_DIR}/bgfx/3rdparty/meshoptimizer/src/indexcodec.cpp
  ${BGFX_ROOT_DIR}/bgfx/3rdparty/meshoptimizer/src/vertexcodec.cpp
  ${BGFX_ROOT_DIR}/bgfx/examples/common/camera.cpp
  ${BGFX_ROOT_DIR}/bgfx/examples/common/bgfx_utils.cpp
  ${BGFX_ROOT_DIR}/bgfx/examples/common/example-glue.cpp
  ${BGFX_ROOT_DIR}/bgfx/examples/common/entry/cmd.cpp
  ${BGFX_ROOT_DIR}/bgfx/examples/common/entry/dialog.cpp
  ${BGFX_ROOT_DIR}/bgfx/examples/common/entry/entry.cpp
  ${BGFX_ROOT_DIR}/bgfx/examples/common/entry/input.cpp
  ${BGFX_ROOT_DIR}/bgfx/examples/common/imgui/imgui.cpp

  ${BGFX_ROOT_DIR}/bgfx/3rdparty/dear-imgui/imgui.cpp
  ${BGFX_ROOT_DIR}/bgfx/3rdparty/dear-imgui/imgui_draw.cpp
  ${BGFX_ROOT_DIR}/bgfx/3rdparty/dear-imgui/imgui_demo.cpp
  ${BGFX_ROOT_DIR}/bgfx/3rdparty/dear-imgui/imgui_widgets.cpp
  ${BGFX_ROOT_DIR}/bgfx/3rdparty/dear-imgui/imgui_tables.cpp
)

if (MSVC)
  set(BGFX_SOURCE_FILES
    ${BX_SOURCES}
    ${BIMG_SOURCES}
    ${BGFX_SOURCES}
    ${BGFX_ROOT_DIR}/bgfx/examples/common/entry/entry_windows.cpp
  )
else()
  set(BGFX_SOURCE_FILES
    ${BX_SOURCES}
    ${BIMG_SOURCES}
    ${BGFX_SOURCES}
    ${BGFX_ROOT_DIR}/bgfx/examples/common/entry/entry_osx.mm
  )
endif()

add_library(bgfx-static STATIC ${BGFX_SOURCE_FILES})

add_custom_command(TARGET bgfx-static
    PRE_BUILD
    COMMAND cd ${BGFX_ROOT_DIR}/bgfx && make -j8 shaderc
)

target_compile_features(bgfx-static PRIVATE cxx_std_17)

# to correctly set __cplusplus macro
if (MSVC)
  target_compile_options(bgfx-static PUBLIC "/Zc:__cplusplus")
  target_include_directories(bgfx-static PUBLIC
    ${BGFX_ROOT_DIR}/bx/include/compat/msvc/
  )
else()
  target_include_directories(bgfx-static PUBLIC
    ${BGFX_ROOT_DIR}/bx/include/compat/osx/
  )
endif()

target_include_directories(bgfx-static PRIVATE
  ${BGFX_ROOT_DIR}/bimg/3rdparty/astc-encoder/include/
  ${BGFX_ROOT_DIR}/bimg/3rdparty
  ${BGFX_ROOT_DIR}/bimg/3rdparty/tinyexr/deps/miniz
  ${BGFX_ROOT_DIR}/bimg/include/
)

add_definitions(
  -D__STDC_LIMIT_MACROS=1
  -D__STDC_FORMAT_MACROS=1
  -D__STDC_CONSTANT_MACROS=1
)

target_compile_definitions(bgfx-static PRIVATE
  "$<$<CONFIG:Debug>:BX_CONFIG_DEBUG=1>"
  "$<$<CONFIG:Release>:BX_CONFIG_DEBUG=0>"
)

target_include_directories(bgfx-static PRIVATE
  ${BGFX_ROOT_DIR}/bx/3rdparty/
)

target_include_directories(bgfx-static PUBLIC
  ${BGFX_ROOT_DIR}/bx/include/
)

target_include_directories(bgfx-static PUBLIC
  ${BGFX_ROOT_DIR}/bgfx/include/
)

target_include_directories(bgfx-static PRIVATE
  ${BGFX_ROOT_DIR}/bgfx/3rdparty/
  ${BGFX_ROOT_DIR}/bgfx/3rdparty/dxsdk/include
  ${BGFX_ROOT_DIR}/bgfx/3rdparty/khronos
  ${BGFX_ROOT_DIR}/bgfx/src/
)

set_source_files_properties(${BGFX_ROOT_DIR}/bimg/3rdparty/tinyexr/deps/miniz/miniz.c
    PROPERTIES LANGUAGE C
)
