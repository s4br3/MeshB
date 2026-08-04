cmake_minimum_required(VERSION 3.14)
project(MeshB LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(FetchContent)
FetchContent_Declare(
    bvh
    GIT_REPOSITORY https://github.com/madmann91/bvh.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(bvh)
FetchContent_Declare(
    CDT
    GIT_REPOSITORY https://github.com/artem-ogre/CDT.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(CDT)
FetchContent_Declare(
    OpenBEM
    GIT_REPOSITORY https://github.com/shash-sharma/OpenBEM.git
    GIT_TAG        main
)
FetchContent_MakeAvailable(OpenBEM)
find_package(OpenMP REQUIRED)
add_executable(main main.cpp)
if (TARGET CDT)
    target_link_libraries(main PRIVATE CDT)
endif()
if (TARGET CDT::CDT)
    target_link_libraries(main PRIVATE CDT::CDT)
endif()
if (TARGET bvh)
    target_link_libraries(main PRIVATE bvh)
endif()
target_include_directories(main PRIVATE
    ${OpenBEM_SOURCE_DIR}
    ${OpenBEM_SOURCE_DIR}/source
    ${cdt_SOURCE_DIR}/CDT/include
)
target_link_libraries(main PRIVATE OpenMP::OpenMP_CXX)