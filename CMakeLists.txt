cmake_minimum_required(VERSION 3.21)

# set the project name
project(mesh_booleans)

# --- Project-wide Settings ---
# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set the target architecture flags.
set(ENABLE_SSE2 True)
set(ENABLE_AVX2 True)


# --- Dependency: oneTBB ---
# Correctly build oneTBB as a subdirectory
set(TBB_TEST OFF CACHE BOOL "Disable TBB tests" FORCE)
set(TBB_EXAMPLES OFF CACHE BOOL "Disable TBB examples" FORCE)
add_subdirectory(arrangements/external/oneTBB)


# --- Dependency: Abseil ---
# CORRECT WAY TO ADD ABSEIL: Build it as a subdirectory.
# CMake will configure it, build its libraries, and create targets for us to link against.
add_subdirectory(arrangements/external/abseil-cpp)


# --- Dependency: Cinolib ---
# Your Cinolib configuration looks good.
set(cinolib_DIR ${PROJECT_SOURCE_DIR}/arrangements/external/Cinolib)
set(CINOLIB_USES_OPENGL_GLFW_IMGUI ON)
set(CINOLIB_USES_SHEWCHUK_PREDICATES ON)
find_package(cinolib REQUIRED)


# --- Define all your executables ---
# We list all targets here to avoid repetition later.
set(ALL_TARGETS
    ${PROJECT_NAME}
    ${PROJECT_NAME}_rotation
    ${PROJECT_NAME}_arap
    ${PROJECT_NAME}_stencil
    ${PROJECT_NAME}_inputcheck
)

add_executable(${PROJECT_NAME}           main.cpp)
add_executable(${PROJECT_NAME}_rotation  main-rotation.cpp)
add_executable(${PROJECT_NAME}_arap      main-arap.cpp)
add_executable(${PROJECT_NAME}_stencil   main-stencil.cpp)
add_executable(${PROJECT_NAME}_inputcheck main-inputcheck.cpp)


# --- Configure ALL targets in a loop to avoid repetition and bugs ---
foreach(TARGET ${ALL_TARGETS})
    # --- Include Directories ---
    # These are common to all targets.
    target_include_directories(${TARGET} PUBLIC
        ./
        code/
        arrangements/code/
        arrangements/external/Indirect_Predicates/include/
    )

    # --- Link Libraries ---
    # Link against the libraries created by our dependencies.
    # Note how we link to specific absl::* targets now.
    target_link_libraries(${TARGET} PUBLIC
        cinolib
        tbb
        absl::time          # For absl::CivilTime and related functions
        absl::base          # Often a required base library
        absl::strings       # Useful for string manipulations
        # Add other absl::* targets if needed, e.g., absl::container
    )

    # --- Compile Definitions ---
    # These definitions are applied to all targets.
    target_compile_definitions(${TARGET} PUBLIC
        TBB_PARALLEL=1
        # [FIX] Add NOMINMAX to prevent conflicts with windows.h macros (min/max)
        $<$<PLATFORM_ID:Windows>:NOMINMAX>
        $<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>
    )

    # --- Compiler-specific Options ---
    # This block now correctly applies to EACH target in the loop.
    if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        target_compile_options(${TARGET} PUBLIC "/fp:strict" "/Oi")
        if(ENABLE_AVX2)
            target_compile_options(${TARGET} PUBLIC "/arch:AVX2")
        elseif(ENABLE_SSE2)
            target_compile_options(${TARGET} PUBLIC "/arch:SSE2")
        endif()
        target_link_options(${TARGET} PUBLIC "/STACK:8421376")
    else() # For GCC/Clang
        target_compile_options(${TARGET} PUBLIC -O2 -frounding-math)
        target_link_options(${TARGET} PUBLIC -Wl,-z,stack-size=8421376)
        if(ENABLE_AVX2)
            target_compile_options(${TARGET} PUBLIC "-mavx2")
        elseif(ENABLE_SSE2)
            target_compile_options(${TARGET} PUBLIC "-msse2")
        endif()
    endif()
endforeach()
