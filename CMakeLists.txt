cmake_minimum_required(VERSION 3.20)
project(excavator_core VERSION 3.0.0 LANGUAGES CXX)

# C++17 required for modern features
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find dependencies
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# Optional: OpenMP for parallelization
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
else()
    message(STATUS "OpenMP not found - building without parallelization support")
endif()

# Source files
set(CPP_SOURCES
    cpp/src/bindings.cpp
    # Will add more as we implement:
    # cpp/src/common.cpp
    # cpp/src/hslm.cpp
    # cpp/src/fastcall.cpp
)

# Create Python module
pybind11_add_module(_excavator_core ${CPP_SOURCES})

# Include directories
target_include_directories(_excavator_core PRIVATE
    cpp/include
    ${Python_INCLUDE_DIRS}
)

# Compiler warnings
if(MSVC)
    target_compile_options(_excavator_core PRIVATE /W4)
else()
    target_compile_options(_excavator_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Optimization flags for Release build
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    if(NOT MSVC)
        target_compile_options(_excavator_core PRIVATE
            -O3
            -march=native
            -ffast-math
        )
    endif()
    message(STATUS "Building in Release mode with optimizations")
endif()

# Link OpenMP if available
if(OpenMP_CXX_FOUND)
    target_link_libraries(_excavator_core PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(_excavator_core PRIVATE EXCAVATOR_HAS_OPENMP)
endif()

# Install target
install(TARGETS _excavator_core LIBRARY DESTINATION excavator2)

# Print configuration summary
message(STATUS "=== EXCAVATOR2 Build Configuration ===")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Python: ${Python_VERSION}")
message(STATUS "pybind11: ${pybind11_VERSION}")
message(STATUS "OpenMP: ${OpenMP_CXX_FOUND}")
message(STATUS "======================================")
