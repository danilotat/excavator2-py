/**
 * @file bindings.cpp
 * @brief Python bindings for EXCAVATOR2 C++ modules
 *
 * This file defines the pybind11 interface that exposes C++ functionality
 * to Python. Currently contains a minimal "hello world" implementation.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>

namespace py = pybind11;

/**
 * Simple test function to verify the build system works
 */
std::string hello_excavator() {
    return "Hello from EXCAVATOR2 C++ module!";
}

/**
 * Get the C++ module version
 */
std::string get_version() {
    return "3.0.0";
}

/**
 * Check if OpenMP support is enabled
 */
bool has_openmp() {
#ifdef EXCAVATOR_HAS_OPENMP
    return true;
#else
    return false;
#endif
}

/**
 * Main pybind11 module definition
 *
 * This module will eventually contain:
 * - hslm: HSLM segmentation algorithm
 * - fastcall: FastCall CNV calling algorithm
 */
PYBIND11_MODULE(_excavator_core, m) {
    m.doc() = "EXCAVATOR2 C++ core algorithms";

    // Module version and info
    m.attr("__version__") = get_version();
    m.def("hello", &hello_excavator, "Test function to verify C++ module works");
    m.def("has_openmp", &has_openmp, "Check if OpenMP support is enabled");

    // TODO: Add HSLM module
    // py::module_ hslm = m.def_submodule("hslm", "HSLM segmentation algorithm");
    // ...

    // TODO: Add FastCall module
    // py::module_ fastcall = m.def_submodule("fastcall", "FastCall CNV calling algorithm");
    // ...
}
