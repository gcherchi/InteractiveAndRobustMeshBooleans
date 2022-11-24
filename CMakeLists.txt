cmake_minimum_required(VERSION 3.21)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# https://cmake.org/cmake/help/git-stage/module/FetchContent.html
include(FetchContent)

# set the project name
project(mesh_booleans)

FetchContent_Declare(arrangements GIT_REPOSITORY "https://github.com/gcherchi/FastAndRobustMeshArrangements.git")
FetchContent_MakeAvailable(arrangements)

# add the executable
add_executable(${PROJECT_NAME} 
    main.cpp
	code/booleans.cpp
	code/foctree.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC ${arrangements_SOURCE_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC ${arrangements_SOURCE_DIR}/external/Cinolib/include)
target_include_directories(${PROJECT_NAME} PUBLIC ${arrangements_SOURCE_DIR}/external/Cinolib/external/eigen)