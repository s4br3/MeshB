cmake_minimum_required(VERSION 3.14)
project(MeshB LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ------------------------------------------------------------------------------
# Dependencies (Header/Lightweight Libraries via FetchContent)
# ------------------------------------------------------------------------------
include(FetchContent)

set(ENABLE_BUILD_SHARED ON CACHE BOOL "" FORCE)
set(ENABLE_SOLVER OFF CACHE BOOL "" FORCE)
set(ENABLE_GAMEPAD OFF CACHE BOOL "" FORCE)
set(ENABLE_DOMHEX OFF CACHE BOOL "" FORCE)
set(ENABLE_QUADMESHINGTOOLS OFF CACHE BOOL "" FORCE)
set(ENABLE_WINSLOWUNTANGLER OFF CACHE BOOL "" FORCE)
set(CDT_USE_AS_COMPILED_LIBRARY OFF CACHE BOOL "" FORCE)

FetchContent_Declare(
    CDT
    GIT_REPOSITORY https://github.com/artem-ogre/CDT.git
    GIT_TAG        763e826
)
FetchContent_Declare(
    OpenBEM
    GIT_REPOSITORY https://github.com/shash-sharma/OpenBEM.git
    GIT_TAG        0b06400
)

FetchContent_MakeAvailable(CDT OpenBEM)

# ------------------------------------------------------------------------------
# System Dependencies (Installed Gmsh & OpenMP)
# ------------------------------------------------------------------------------
find_package(OpenMP REQUIRED)

# Locate system-installed Gmsh dynamic/static library and C++ header
find_library(GMSH_LIBRARY NAMES gmsh REQUIRED)
find_path(GMSH_INCLUDE_DIR NAMES gmsh.h PATH_SUFFIXES gmsh REQUIRED)

# ------------------------------------------------------------------------------
# Core Library Target
# ------------------------------------------------------------------------------
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

add_library(meshb_core STATIC
    ${SOURCE_DIR}/bem_interop.cpp
    ${SOURCE_DIR}/boolean_ops.cpp
    ${SOURCE_DIR}/bvh_collisions.cpp
    ${SOURCE_DIR}/gmsh_interop.cpp
    ${SOURCE_DIR}/math_utils.cpp
    ${SOURCE_DIR}/mesh_clean.cpp
    ${SOURCE_DIR}/mesh_types.cpp
    ${SOURCE_DIR}/raycast.cpp
    ${SOURCE_DIR}/triangulation.cpp
    ${SOURCE_DIR}/vector.cpp
)

target_include_directories(meshb_core PUBLIC
    ${SOURCE_DIR}
    ${OpenBEM_SOURCE_DIR}
    ${OpenBEM_SOURCE_DIR}/source
    ${cdt_SOURCE_DIR}/CDT/include
    ${GMSH_INCLUDE_DIR}
)

target_link_libraries(meshb_core PUBLIC
    ${GMSH_LIBRARY}
    OpenMP::OpenMP_CXX
)

if (TARGET shared)
    target_link_libraries(meshb_core PUBLIC shared)
endif()

# ------------------------------------------------------------------------------
# Main Executable
# ------------------------------------------------------------------------------
add_executable(main ${SOURCE_DIR}/main.cpp)
target_link_libraries(main PRIVATE meshb_core)

# ------------------------------------------------------------------------------
# Testing
# ------------------------------------------------------------------------------
enable_testing()

macro(add_mesh_test TEST_NAME TEST_FILE)
    add_executable(${TEST_NAME} test/${TEST_FILE})
    target_link_libraries(${TEST_NAME} PRIVATE meshb_core)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endmacro()

add_mesh_test(test_vector test_vector.cpp)
add_mesh_test(test_math_utils test_math_utils.cpp)
add_mesh_test(test_bvh test_bvh.cpp)
add_mesh_test(test_mesh_clean test_mesh_clean.cpp)
add_mesh_test(test_triangulation test_triangulation.cpp)
add_mesh_test(test_raycast test_raycast.cpp)