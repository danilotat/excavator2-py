cmake_minimum_required(VERSION 3.20)
project(excavator2 LANGUAGES CXX)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_core MODULE cpp/bindings/module.cpp cpp/fastcall.cpp cpp/hslm.cpp cpp/bindings/hslm.cpp)
target_include_directories(_core PRIVATE cpp/include)
target_compile_features(_core PRIVATE cxx_std_17)
set_target_properties(_core PROPERTIES CXX_EXTENSIONS OFF)

# Compatibility first: no fast-math, FP reassociation, or host-specific ISA.
if(MSVC)
    target_compile_options(_core PRIVATE /fp:strict)
else()
    target_compile_options(_core PRIVATE -fno-fast-math -ffp-contract=off)
endif()

install(TARGETS _core LIBRARY DESTINATION excavator2 RUNTIME DESTINATION excavator2)
