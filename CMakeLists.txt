cmake_minimum_required(VERSION 3.20)
project(excavator2 LANGUAGES CXX)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_core MODULE cpp/bindings/module.cpp cpp/fastcall.cpp cpp/hslm.cpp cpp/bindings/hslm.cpp cpp/bindings/reads.cpp)
# Use the same HTSlib headers/library as the pinned runtime pysam dependency.
execute_process(COMMAND "${Python_EXECUTABLE}" -c
    "import pysam; print(';'.join(pysam.get_include()))"
    OUTPUT_VARIABLE PYSAM_INCLUDES OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
execute_process(COMMAND "${Python_EXECUTABLE}" -c
    "import pysam.libchtslib; print(pysam.libchtslib.__file__)"
    OUTPUT_VARIABLE PYSAM_HTSLIB OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
target_include_directories(_core PRIVATE cpp/include ${PYSAM_INCLUDES})
# Link by filename on ELF to avoid recording the temporary build environment path.
get_filename_component(PYSAM_LIBRARY_DIR "${PYSAM_HTSLIB}" DIRECTORY)
get_filename_component(PYSAM_LIBRARY_NAME "${PYSAM_HTSLIB}" NAME)
if(APPLE)
    target_link_libraries(_core PRIVATE "${PYSAM_HTSLIB}")
    set_target_properties(_core PROPERTIES INSTALL_RPATH "@loader_path/../pysam")
else()
    target_link_directories(_core PRIVATE "${PYSAM_LIBRARY_DIR}")
    target_link_libraries(_core PRIVATE "-l:${PYSAM_LIBRARY_NAME}")
    set_target_properties(_core PROPERTIES INSTALL_RPATH "$ORIGIN/../pysam")
endif()
target_compile_features(_core PRIVATE cxx_std_17)
set_target_properties(_core PROPERTIES CXX_EXTENSIONS OFF)

# Compatibility first: no fast-math, FP reassociation, or host-specific ISA.
if(MSVC)
    target_compile_options(_core PRIVATE /fp:strict)
else()
    target_compile_options(_core PRIVATE -fno-fast-math -ffp-contract=off)
endif()

install(TARGETS _core LIBRARY DESTINATION excavator2 RUNTIME DESTINATION excavator2)
