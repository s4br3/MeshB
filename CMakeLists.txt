cmake_minimum_required(VERSION 3.14)
project(MeshB LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(FetchContent)
set(ENABLE_BUILD_SHARED ON CACHE BOOL "" FORCE)
set(CDT_USE_AS_COMPILED_LIBRARY OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
    bvh
    GIT_REPOSITORY https://github.com/madmann91/bvh.git
    GIT_TAG        master
)
FetchContent_Declare(
    CDT
    GIT_REPOSITORY https://github.com/artem-ogre/CDT.git
    GIT_TAG        master
)
FetchContent_Declare(
    OpenBEM
    GIT_REPOSITORY https://github.com/shash-sharma/OpenBEM.git
    GIT_TAG        main
)
FetchContent_Declare(
    gmsh
    GIT_REPOSITORY https://github.com/live-clones/gmsh.git
    GIT_TAG        master
)

FetchContent_MakeAvailable(bvh CDT OpenBEM gmsh)

find_package(OpenMP REQUIRED)

add_executable(main
    main.cpp
    boolean_ops.cpp
    bvh_collisions.cpp
    math_utils.cpp
    mesh_types.cpp
    triangulation.cpp
    raycast.cpp
    bem_interop.cpp
    gmsh_interop.cpp
    mesh_clean.cpp
)
target_include_directories(main PRIVATE
    ${OpenBEM_SOURCE_DIR}
    ${OpenBEM_SOURCE_DIR}/source
    ${cdt_SOURCE_DIR}/CDT/include
    ${gmsh_SOURCE_DIR}
)
target_link_libraries(main PRIVATE
    OpenMP::OpenMP_CXX
)

if (TARGET shared)
    target_link_libraries(main PRIVATE shared)
endif()

if (TARGET bvh)
    target_link_libraries(main PRIVATE bvh)
endif()