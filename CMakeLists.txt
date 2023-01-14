cmake_minimum_required(VERSION 3.0)

set(CMAKE_XCODE_ATTRIBUTE_ARCHS "$(ARCHS_STANDARD)" CACHE STRING "" FORCE)
set(CMAKE_OSX_ARCHITECTURES "$(ARCHS_STANDARD)" CACHE STRING "" FORCE)
set(CMAKE_XCODE_ATTRIBUTE_DEBUG_INFORMATION_FORMAT "dwarf")
set(CMAKE_XCODE_ATTRIBUTE_DEPLOYMENT_POSTPROCESSING "NO")

project(coverage-occlusion VERSION 0.18 LANGUAGES CXX)

SET(CMAKE_CONFIGURATION_TYPES Debug Release)
enable_language(C)
enable_language(CXX)

set(SOURCES
bgfx.cmake
occlusion.cpp
rasterizer.h
rasterizer.cpp
rasterizer_box.cpp
rasterizer_sort.cpp
rasterizer_draw.cpp
rasterizer_math.h
rasterizer_math_sse.inl
rasterizer_math_neon.inl
rasterizer_tile.h
rasterizer_tile.cpp
)

set(SOURCES_ENKI
enkiTS/src/TaskScheduler.h
enkiTS/src/TaskScheduler.cpp
)

if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-g -ffast-math -Wall -Wextra -Wshadow -Wno-missing-field-initializers)
endif()

include(bgfx.cmake)

add_executable(demo ${SOURCES} ${SOURCES_ENKI})
target_link_libraries(demo bgfx-static)

target_compile_features(demo PRIVATE cxx_std_17)

target_compile_definitions(demo PRIVATE
  "ENTRY_CONFIG_IMPLEMENT_MAIN=1"
  "$<$<CONFIG:Debug>:BX_CONFIG_DEBUG=1>"
  "$<$<CONFIG:Release>:BX_CONFIG_DEBUG=0>"
)

if (MSVC)
else()
  target_link_libraries(demo
    "-framework Metal"
    "-framework Carbon"
    "-framework QuartzCore"
    "-framework AppKit"
    "-framework Cocoa"
    "-framework IOKit"
  )
  set_target_properties(demo
                      PROPERTIES
                      XCODE_ATTRIBUTE_ONLY_ACTIVE_ARCH "YES"
  )
endif()

add_custom_command(TARGET demo
    PRE_BUILD
    COMMAND cd ${BGFX_ROOT_DIR} && make rebuild
)

target_include_directories(demo PUBLIC
  ${BGFX_ROOT_DIR}/bx/include/
  ${BGFX_ROOT_DIR}/bgfx/include/
  ${BGFX_ROOT_DIR}/bgfx/3rdparty/
  ${BGFX_ROOT_DIR}/bimg/include/
  ${BGFX_ROOT_DIR}/bgfx/examples/common/
  ${CMAKE_CURRENT_SOURCE_DIR}/enkiTS/src
)
