cmake_minimum_required(VERSION 3.14)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

FetchContent_Declare(
    OpenBEM
    GIT_REPOSITORY https://github.com/shash-sharma/OpenBEM.git
    GIT_TAG        0b06400
)
FetchContent_MakeAvailable(OpenBEM)

find_package(OpenMP REQUIRED)

find_path(GMSH_INCLUDE_DIR NAMES gmsh.h)
find_library(GMSH_LIBRARY NAMES gmsh)

if(NOT GMSH_INCLUDE_DIR OR NOT GMSH_LIBRARY)
    message(FATAL_ERROR "Gmsh headers or library not found! Ensure libgmsh-dev is installed.")
endif()

set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

add_library(meshb_core STATIC
    ${SOURCE_DIR}/bem_interop.cpp
    ${SOURCE_DIR}/boolean_ops.cpp
    ${SOURCE_DIR}/bvh_collisions.cpp
    ${SOURCE_DIR}/mesh_io.cpp
    ${SOURCE_DIR}/math_utils.cpp
    ${SOURCE_DIR}/mesh_clean.cpp
    ${SOURCE_DIR}/mesh_types.cpp
    ${SOURCE_DIR}/classify.cpp
    ${SOURCE_DIR}/triangulation.cpp
    ${SOURCE_DIR}/vector.cpp
    ${SOURCE_DIR}/geom_2d.cpp
    ${SOURCE_DIR}/geom_3d.cpp
)

target_include_directories(meshb_core PUBLIC
    ${SOURCE_DIR}
    ${OpenBEM_SOURCE_DIR}
    ${OpenBEM_SOURCE_DIR}/source
    ${GMSH_INCLUDE_DIR}
)

target_link_libraries(meshb_core PUBLIC
    OpenMP::OpenMP_CXX
    ${GMSH_LIBRARY}
)

add_executable(main ${SOURCE_DIR}/main.cpp)
target_link_libraries(main PRIVATE meshb_core)

enable_testing()
macro(add_mesh_test TEST_NAME TEST_FILE)
    add_executable(${TEST_NAME} test/${TEST_FILE})
    target_link_libraries(${TEST_NAME} PRIVATE meshb_core)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endmacro()

add_mesh_test(test_vector test_vector.cpp)
add_mesh_test(test_math_utils test_math_utils.cpp)
add_mesh_test(test_bvh test_bvh.cpp)
add_mesh_test(test_geom_2d test_geom_2d.cpp)
add_mesh_test(test_geom_3d test_geom_3d.cpp)
add_mesh_test(test_mesh_clean test_mesh_clean.cpp)
add_mesh_test(test_triangulation test_triangulation.cpp)
add_mesh_test(test_classify test_classify.cpp)