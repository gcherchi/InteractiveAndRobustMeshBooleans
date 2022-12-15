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
add_executable(${PROJECT_NAME} main.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC
    ./
    code/
    arrangements/code/
)

set(cinolib_DIR ${PROJECT_SOURCE_DIR}/arrangements/external/Cinolib)
set(CINOLIB_USES_OPENGL_GLFW_IMGUI ON)
set(CINOLIB_USES_SHEWCHUK_PREDICATES ON)
set(CINOLIB_USES_INDIRECT_PREDICATES ON)

find_package(cinolib REQUIRED)

target_link_libraries(${PROJECT_NAME} cinolib)
target_link_libraries(${PROJECT_NAME} tbb)
target_compile_definitions(${PROJECT_NAME} PUBLIC TBB_PARALLEL=1)
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/abseil-cpp/)
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/oneTBB/)


add_executable(${PROJECT_NAME}_rotation main-rotation.cpp)

target_include_directories(${PROJECT_NAME}_rotation PUBLIC
    ./
    code/
    arrangements/code/
)

target_link_libraries(${PROJECT_NAME}_rotation cinolib)
target_link_libraries(${PROJECT_NAME}_rotation tbb)
target_compile_definitions(${PROJECT_NAME}_rotation PUBLIC TBB_PARALLEL=1)
target_include_directories(${PROJECT_NAME}_rotation PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/abseil-cpp/)
target_include_directories(${PROJECT_NAME}_rotation PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/oneTBB/)

add_executable(${PROJECT_NAME}_arap main-arap.cpp)

target_include_directories(${PROJECT_NAME}_arap PUBLIC
    ./
    code/
    arrangements/code/
)

target_link_libraries(${PROJECT_NAME}_arap cinolib)
target_link_libraries(${PROJECT_NAME}_arap tbb)
target_compile_definitions(${PROJECT_NAME}_arap PUBLIC TBB_PARALLEL=1)
target_include_directories(${PROJECT_NAME}_arap PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/abseil-cpp/)
target_include_directories(${PROJECT_NAME}_arap PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/oneTBB/)

add_executable(${PROJECT_NAME}_stencil main-stencil.cpp)

target_include_directories(${PROJECT_NAME}_stencil PUBLIC
    ./
    code/
    arrangements/code/
)

target_link_libraries(${PROJECT_NAME}_stencil cinolib)
target_link_libraries(${PROJECT_NAME}_stencil tbb)
target_compile_definitions(${PROJECT_NAME}_stencil PUBLIC TBB_PARALLEL=1)
target_include_directories(${PROJECT_NAME}_stencil PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/abseil-cpp/)
target_include_directories(${PROJECT_NAME}_stencil PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/oneTBB/)


add_executable(${PROJECT_NAME}_inputcheck main-inputcheck.cpp)

target_include_directories(${PROJECT_NAME}_inputcheck PUBLIC
        ./
        code/
        arrangements/code/
        )

target_link_libraries(${PROJECT_NAME}_inputcheck cinolib)
target_link_libraries(${PROJECT_NAME}_inputcheck tbb)
target_compile_definitions(${PROJECT_NAME}_inputcheck PUBLIC TBB_PARALLEL=1)
target_include_directories(${PROJECT_NAME}_inputcheck PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/abseil-cpp/)
target_include_directories(${PROJECT_NAME}_inputcheck PUBLIC ${PROJECT_SOURCE_DIR}/arrangements/external/oneTBB/)
