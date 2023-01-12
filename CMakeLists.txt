cmake_minimum_required(VERSION 3.0)

project(coverage-occlusion VERSION 0.18 LANGUAGES CXX)

SET(CMAKE_CONFIGURATION_TYPES Debug Release)

set(SOURCES
bgfx.cmake
occlusion.cpp
rasterizer.h
rasterizer.cpp
rasterizer_math.h
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
    add_compile_options(-Wall -Wextra -Wshadow -Wno-missing-field-initializers)
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
