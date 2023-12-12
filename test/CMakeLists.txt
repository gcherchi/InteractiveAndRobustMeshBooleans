cmake_minimum_required(VERSION 3.21)

project(test_code)

set(cinolib_DIR "/Users/gianmarco/Code/Cinolib")
set(CINOLIB_USES_OPENGL_GLFW_IMGUI ON)
set(CINOLIB_USES_SHEWCHUK_PREDICATES ON)

find_package(cinolib REQUIRED)

add_executable(${PROJECT_NAME}
        main.cpp
)

target_link_libraries(${PROJECT_NAME} cinolib)