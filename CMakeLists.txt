cmake_minimum_required(VERSION 3.21)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set the project name
project(mesh_booleans)

set(TBB_TEST OFF CACHE BOOL " " FORCE)
set(TBB_EXAMPLES OFF CACHE BOOL " " FORCE)
add_subdirectory(arrangements/external/oneTBB)

# add the executable
add_executable(${PROJECT_NAME}
    main.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC
    ./
    code/
    arrangements/code/
)

target_link_libraries(${PROJECT_NAME} tbb)
target_compile_definitions(${PROJECT_NAME} PUBLIC TBB_PARALLEL=1)

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/abseil-cpp/)
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/oneTBB/)

set(cinolib_DIR ${PROJECT_SOURCE_DIR}/arrangements/external/Cinolib)
set(CINOLIB_USES_OPENGL_GLFW_IMGUI ON)
set(CINOLIB_USES_SHEWCHUK_PREDICATES ON)
set(CINOLIB_USES_INDIRECT_PREDICATES ON)

find_package(cinolib REQUIRED)
target_link_libraries(${PROJECT_NAME} cinolib)