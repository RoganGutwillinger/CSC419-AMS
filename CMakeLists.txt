cmake_minimum_required(VERSION 3.1)
project(ams)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

### libIGL options: choose between header only and compiled static library
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" OFF)
option(LIBIGL_GLFW       "Use GLFW"           ON)

include(libigl)

# Add your project files
include_directories("include/")
file(GLOB SRCFILES src/*.cpp)

add_library(core ${SRCFILES})
target_link_libraries(core igl::core igl::glfw )

add_executable(ams "main.cpp")
target_link_libraries(ams core igl::core igl::opengl igl::glfw)