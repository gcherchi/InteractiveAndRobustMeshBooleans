cmake_minimum_required(VERSION 3.21)

# Set the target architecture.
# All modern x86/x64 processors support AVX2.
# Older x86/x64 processors may support SSE2 but not AVX2.
# Very old x86/x64 processors, or non x86/x64
# processors, do not support any of the two.
set(ENABLE_SSE2 True)
set(ENABLE_AVX2 True)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set the project name
project(mesh_booleans)

# External dependencies
set(TBB_TEST OFF CACHE BOOL " " FORCE)
set(TBB_EXAMPLES OFF CACHE BOOL " " FORCE)
add_subdirectory(arrangements/external/oneTBB)

# cinolib configuration
set(cinolib_DIR ${PROJECT_SOURCE_DIR}/arrangements/external/Cinolib)
set(CINOLIB_USES_OPENGL_GLFW_IMGUI   ON)
set(CINOLIB_USES_CGAL_GMP_MPFR       ON)
set(CINOLIB_USES_SHEWCHUK_PREDICATES ON)
set(CINOLIB_USES_INDIRECT_PREDICATES ON)

find_package(CGAL REQUIRED COMPONENTS Core)
find_package(cinolib REQUIRED)

# Common include directories
set(COMMON_INCLUDE_DIRS
        ./
        code/
        arrangements/code/
        rationals_code/
        ${PROJECT_SOURCE_DIR}/arrangements/external/abseil-cpp/
        ${PROJECT_SOURCE_DIR}/arrangements/external/oneTBB/
)

# Common source files
set(RATIONALS_SOURCES
        rationals_code/intersect_custom.h
        rationals_code/intersect_custom.cpp
        rationals_code/intersect_point_rationals.h
        rationals_code/intersect_point_rationals.cpp
)

# Common libraries and definitions
set(COMMON_LIBS cinolib tbb CGAL::CGAL)
set(COMMON_DEFINES TBB_PARALLEL=1)

# Function to reduce duplication
function(add_mesh_executable name main_cpp)
        add_executable(${name} ${main_cpp} ${RATIONALS_SOURCES})
        target_include_directories(${name} PUBLIC ${COMMON_INCLUDE_DIRS})
        target_link_libraries(${name} ${COMMON_LIBS})
        target_compile_definitions(${name} PUBLIC ${COMMON_DEFINES})
endfunction()

# Add all executables
add_mesh_executable(${PROJECT_NAME}                     main.cpp)
add_mesh_executable(${PROJECT_NAME}_rotation            main-rotation.cpp)
add_mesh_executable(${PROJECT_NAME}_arap                main-arap.cpp)
add_mesh_executable(${PROJECT_NAME}_stencil             main-stencil.cpp)
add_mesh_executable(${PROJECT_NAME}_inputcheck          main-inputcheck.cpp)
add_mesh_executable(${PROJECT_NAME}_test-rationals      main-test-rationals.cpp)
add_mesh_executable(${PROJECT_NAME}_test                main_test.cpp)
add_mesh_executable(${PROJECT_NAME}_visualizer          main_visualizer.cpp)
add_mesh_executable(${PROJECT_NAME}_test_paper          main_test_paper.cpp)

# Compiler-specific options
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        # grant IEEE 754 compliance
        target_compile_options(${PROJECT_NAME} PUBLIC "/fp:strict")
        # use intrinsic functions
        target_compile_options(${PROJECT_NAME} PUBLIC "/Oi")
        # set target architecture
        if(ENABLE_AVX2)
                target_compile_options(${PROJECT_NAME} PUBLIC "/arch:AVX2")
        elseif(ENABLE_SSE2)
                target_compile_options(${PROJECT_NAME} PUBLIC "/arch:SSE2")
        endif()
        # reserve enough stack size
        target_link_options(${PROJECT_NAME} PUBLIC "/STACK:8421376")
        # turn off annoying warnings
        target_compile_options(${PROJECT_NAME} PUBLIC "/D _CRT_SECURE_NO_WARNINGS")
else()

        # set standard optimization level
        #target_compile_options(${PROJECT_NAME} PUBLIC -O3)
        # reserve enough stack size
        target_compile_options(${PROJECT_NAME} PUBLIC -Wl,-z,stacksize=8421376)
        # grant IEEE 754 compliance
        target_compile_options(${PROJECT_NAME} PUBLIC -frounding-math)
        # set target architecture
        if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|amd64")
                if(ENABLE_AVX2)
                        target_compile_options(${PROJECT_NAME} PUBLIC "-mavx2")
                elseif(ENABLE_SSE2)
                        target_compile_options(${PROJECT_NAME} PUBLIC "-msse2")
                endif()
        endif ()
endif()
