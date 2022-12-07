cmake_minimum_required(VERSION 3.21)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set the project name
project(mesh_booleans)

set(TBB_TEST OFF CACHE BOOL " " FORCE)
set(TBB_EXAMPLES OFF CACHE BOOL " " FORCE)
add_subdirectory(arrangements/external/oneTBB)

set(cinolib_DIR ${PROJECT_SOURCE_DIR}/arrangements/external/Cinolib)
set(CINOLIB_USES_OPENGL_GLFW_IMGUI ON)
set(CINOLIB_USES_SHEWCHUK_PREDICATES ON)
set(CINOLIB_USES_INDIRECT_PREDICATES ON)

find_package(cinolib REQUIRED)

set(IRMB_INCLUDE_DIRS
	${PROJECT_SOURCE_DIR}
	${PROJECT_SOURCE_DIR}/code
	${PROJECT_SOURCE_DIR}/arrangements/code
	${PROJECT_SOURCE_DIR}/arrangements/external/abseil-cpp
	${PROJECT_SOURCE_DIR}/arrangements/external/oneTBB
	)

set(IRMB_LIBS cinolib tbb)

# add the executable
add_executable(${PROJECT_NAME} main.cpp)
target_include_directories(${PROJECT_NAME} PUBLIC ${IRMB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${IRMB_LIBS})
target_compile_definitions(${PROJECT_NAME} PUBLIC TBB_PARALLEL=1)

add_executable(${PROJECT_NAME}_rotation main-rotation.cpp)
target_include_directories(${PROJECT_NAME}_rotation PUBLIC ${IRMB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME}_rotation ${IRMB_LIBS})
target_compile_definitions(${PROJECT_NAME}_rotation PUBLIC TBB_PARALLEL=1)

add_executable(${PROJECT_NAME}_arap main-arap.cpp)
target_include_directories(${PROJECT_NAME}_arap PUBLIC ${IRMB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME}_arap ${IRMB_LIBS})
target_compile_definitions(${PROJECT_NAME}_arap PUBLIC TBB_PARALLEL=1)

add_executable(${PROJECT_NAME}_stencil main-stencil.cpp)
target_include_directories(${PROJECT_NAME}_stencil PUBLIC ${IRMB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME}_stencil ${IRMB_LIBS})
target_compile_definitions(${PROJECT_NAME}_stencil PUBLIC TBB_PARALLEL=1)

