cmake_minimum_required(VERSION 3.27)
project(triangle_intersections)

set(CMAKE_CXX_STANDARD 20)

set(cinolib_DIR "/home/luca/lib/cinolib")


set(CINOLIB_USES_OPENGL_GLFW_IMGUI    ON  ) # used for rendering with OpenGL
set(CINOLIB_USES_TETGEN               OFF ) # used for polygon triangulation
set(CINOLIB_USES_TRIANGLE             OFF ) # used for tetrahedralization
set(CINOLIB_USES_SHEWCHUK_PREDICATES  ON  ) # used for exact geometric tests on input points
set(CINOLIB_USES_INDIRECT_PREDICATES  ON ) # used for exact geometric tests on implicit points
set(CINOLIB_USES_GRAPH_CUT            OFF ) # used for graph clustering
set(CINOLIB_USES_BOOST                OFF ) # used for 2D polygon operations (e.g. thickening, clipping, 2D booleans...)
set(CINOLIB_USES_VTK                  OFF ) # used just to support VTK file formats
set(CINOLIB_USES_SPECTRA              OFF ) # used for matrix eigendecomposition
set(CINOLIB_USES_CGAL_GMP_MPFR        ON ) # used for rational numbers with lazy kernel

set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)

find_package(cinolib REQUIRED)
find_package(CGAL REQUIRED)
add_executable(
        ${PROJECT_NAME} main.cpp
        orientTemplated.h
        tri_tri_intersect.h
        printIntersections.cpp
        printIntersections.h

)

target_link_libraries(${CMAKE_PROJECT_NAME} CGAL::CGAL cinolib)
