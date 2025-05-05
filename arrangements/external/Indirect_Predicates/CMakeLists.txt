cmake_minimum_required(VERSION 3.10)

# Set the target architecture.
# All modern x86/x64 processors support AVX2.
# Older x86/x64 processors may support SSE2 but not AVX2.
# Very old x86/x64 processors, or non x86/x64 
# processors, do not support any of the two.
# Modern ARM processors support NEON.

# --- Uncomment max one of the following lines ---

#							#Generic build
#set(ENABLE_SSE2 True)		#Build for x86-64 with SSE2
set(ENABLE_AVX2 True)		#Build for x86-64 with AVX2
#set(ENABLE_ARMNEON True)	#Build for ARM with NEON

# --------- DON NOT EDIT BELOW THIS LINE! --------

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# disable SSE and AVX on MacOS
if(APPLE)
set(ENABLE_SSE2 False)
set(ENABLE_AVX2 False)
endif()

# AVX includes SSE -> make that explicit
if(ENABLE_AVX2)
set(ENABLE_SSE2 True)
endif()

# Check that no incompatible architectures are slelected
if(ENABLE_ARMNEON AND ENABLE_SSE2)
message(FATAL_ERROR "ARMNEON and SSE2/AVX2 cannot be both set!")
endif()

# set the project name
project(indirectPredicates)

# add the executable
add_executable(${PROJECT_NAME} 
    test.cpp
)

# Compiler-specific options
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	# grant IEEE 754 compliance; remove security controls; turn off annoying warnings
	target_compile_options(${PROJECT_NAME} PUBLIC "/fp:strict;/GS-;/D _CRT_SECURE_NO_WARNINGS")
	# set target architecture
	if(ENABLE_AVX2)
		target_compile_options(${PROJECT_NAME} PUBLIC "/arch:AVX2")
	elseif(ENABLE_SSE2)
		target_compile_options(${PROJECT_NAME} PUBLIC "/arch:SSE2")
	elseif(ENABLE_ARMNEON)
		target_compile_options(${PROJECT_NAME} PUBLIC "/arch:ARMV8.0")
	endif()
else()
	# set max optimization
	target_compile_options(${PROJECT_NAME} PUBLIC -O3)
	# grant IEEE 754 compliance
	target_compile_options(${PROJECT_NAME} PUBLIC -frounding-math)
	# set target architecture
	if(ENABLE_AVX2)
		target_compile_options(${PROJECT_NAME} PUBLIC "-mavx2")
		target_compile_options(${PROJECT_NAME} PUBLIC "-mfma")
	elseif(ENABLE_SSE2)
		target_compile_options(${PROJECT_NAME} PUBLIC "-msse2")
	elseif(ENABLE_ARMNEON)
		target_compile_options(${PROJECT_NAME} PUBLIC "-march=armv8-a+simd")
	endif()
endif()

# Public include directory
target_include_directories(${PROJECT_NAME} PUBLIC 
	.
	include
)

if(ENABLE_ARMNEON)
	include(FetchContent)

	FetchContent_Declare(
		simde
		URL https://github.com/simd-everywhere/simde/archive/refs/heads/master.zip)
	FetchContent_Populate(simde)
	message(STATUS "simde is available in " ${simde_SOURCE_DIR})

	target_include_directories(${PROJECT_NAME} PUBLIC ${simde_SOURCE_DIR}/simde)
	target_compile_options(${PROJECT_NAME} PUBLIC -DSIMDE_ENABLE_NATIVE_ALIASES)
endif()
