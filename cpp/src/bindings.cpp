/**
 * @file bindings.cpp
 * @brief Python bindings for EXCAVATOR2 C++ modules
 *
 * This file defines the pybind11 interface that exposes C++ functionality
 * to Python.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <string>
#include <vector>

#include "excavator/common.hpp"
#include "excavator/hslm.hpp"
#include "excavator/fastcall.hpp"

namespace py = pybind11;
using namespace excavator;

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

#ifdef EXCAVATOR_HAS_OPENMP
#include <omp.h>
#endif

/**
 * Get number of OpenMP threads
 */
int get_num_threads() {
#ifdef EXCAVATOR_HAS_OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

/**
 * Main pybind11 module definition
 */
PYBIND11_MODULE(_excavator_core, m) {
    m.doc() = R"pbdoc(
        EXCAVATOR2 C++ core algorithms
        ------------------------------

        This module provides the core computational algorithms for EXCAVATOR2:

        - hslm: Heterogeneous Shifting Level Model segmentation
        - fastcall: FastCall CNV classification (coming soon)

        These algorithms are implemented in C++ for performance and exposed
        to Python via pybind11.
    )pbdoc";

    // Module version and info
    m.attr("__version__") = get_version();
    m.def("hello", &hello_excavator, "Test function to verify C++ module works");
    m.def("has_openmp", &has_openmp, "Check if OpenMP support is enabled");
    m.def("get_num_threads", &get_num_threads, "Get number of OpenMP threads");

    // ==================== HSLM Submodule ====================
    py::module_ hslm_mod = m.def_submodule("hslm", R"pbdoc(
        HSLM (Heterogeneous Shifting Level Model) segmentation algorithm.

        This module implements an inhomogeneous Hidden Markov Model for
        detecting copy number breakpoints in log2-ratio data. The transition
        probabilities depend on genomic distance between consecutive observations.

        Example:
            >>> from excavator2._excavator_core import hslm
            >>> params = hslm.HSLMParameters()
            >>> segmenter = hslm.HSLM(params)
            >>> result = segmenter.segment(log2_ratios, positions)
            >>> print(f"Found {result.n_segments} segments")
    )pbdoc");

    // HSLMParameters struct
    py::class_<hslm::HSLMParameters>(hslm_mod, "HSLMParameters", R"pbdoc(
        Parameters for the HSLM segmentation algorithm.

        Attributes:
            omega (float): Variance partitioning factor (0-1). Default: 0.1
                           Fraction of variance attributed to states vs noise.
            theta (float): Base transition probability. Default: 1e-5
            d_norm (float): Distance normalization factor. Default: 1e6
            step_eta (float): Distance step for eta computation (bp). Default: 200000
            n_states (int): Number of hidden states. Default: 21
            min_segment_size (int): Minimum segment size to keep. Default: 1
    )pbdoc")
        .def(py::init<>())
        .def_readwrite("omega", &hslm::HSLMParameters::omega)
        .def_readwrite("theta", &hslm::HSLMParameters::theta)
        .def_readwrite("d_norm", &hslm::HSLMParameters::d_norm)
        .def_readwrite("step_eta", &hslm::HSLMParameters::step_eta)
        .def_readwrite("n_states", &hslm::HSLMParameters::n_states)
        .def_readwrite("min_segment_size", &hslm::HSLMParameters::min_segment_size)
        .def("__repr__", [](const hslm::HSLMParameters& p) {
            return "<HSLMParameters omega=" + std::to_string(p.omega) +
                   " theta=" + std::to_string(p.theta) +
                   " n_states=" + std::to_string(p.n_states) + ">";
        });

    // HSLMResult struct
    py::class_<hslm::HSLMResult>(hslm_mod, "HSLMResult", R"pbdoc(
        Result of HSLM segmentation.

        Attributes:
            breakpoints (list[int]): Indices where segments change (0-indexed).
                                     Includes 0 at start and len(data) at end.
            segment_means (list[float]): Mean value for each segment.
            state_path (list[int]): State sequence from Viterbi algorithm.
            state_values (list[float]): The actual values corresponding to each state.
            n_segments (int): Number of segments detected.
            success (bool): Whether the algorithm succeeded.
            error_message (str): Error message if failed.
    )pbdoc")
        .def(py::init<>())
        .def_readonly("breakpoints", &hslm::HSLMResult::breakpoints)
        .def_readonly("segment_means", &hslm::HSLMResult::segment_means)
        .def_readonly("state_path", &hslm::HSLMResult::state_path)
        .def_readonly("state_values", &hslm::HSLMResult::state_values)
        .def_readonly("n_segments", &hslm::HSLMResult::n_segments)
        .def_readonly("success", &hslm::HSLMResult::success)
        .def_readonly("error_message", &hslm::HSLMResult::error_message)
        .def("__repr__", [](const hslm::HSLMResult& r) {
            return "<HSLMResult n_segments=" + std::to_string(r.n_segments) +
                   " success=" + (r.success ? "True" : "False") + ">";
        });

    // HSLM class
    py::class_<hslm::HSLM>(hslm_mod, "HSLM", R"pbdoc(
        HSLM (Heterogeneous Shifting Level Model) segmentation algorithm.

        This class implements an inhomogeneous HMM for detecting breakpoints
        in log2-ratio data from copy number analysis.

        Example:
            >>> params = HSLMParameters()
            >>> params.omega = 0.1
            >>> params.theta = 1e-5
            >>> segmenter = HSLM(params)
            >>> result = segmenter.segment(log2_ratios, positions)

        Args:
            params: HSLMParameters object with algorithm settings
    )pbdoc")
        .def(py::init<const hslm::HSLMParameters&>(),
             py::arg("params") = hslm::HSLMParameters())
        .def("segment", &hslm::HSLM::segment,
             py::arg("log2_ratios"),
             py::arg("positions"),
             R"pbdoc(
                Run segmentation on a single chromosome.

                Args:
                    log2_ratios: Vector of log2 ratio values (test/control)
                    positions: Genomic positions corresponding to each log2 ratio

                Returns:
                    HSLMResult containing breakpoints and segment information
             )pbdoc")
        .def("segment_multi", &hslm::HSLM::segment_multi,
             py::arg("data_matrix"),
             py::arg("positions"),
             R"pbdoc(
                Run segmentation on multiple sequences (multi-sample mode).

                Args:
                    data_matrix: List of log2 ratio vectors, one per sample
                    positions: Genomic positions (same for all samples)

                Returns:
                    HSLMResult containing breakpoints
             )pbdoc")
        .def_property("params",
            &hslm::HSLM::params,
            &hslm::HSLM::set_params,
            "HSLM parameters");

    // Convenience function for single-call segmentation
    hslm_mod.def("segment", [](
        const std::vector<double>& log2_ratios,
        const std::vector<int64_t>& positions,
        double omega,
        double theta,
        double step_eta,
        int n_states,
        int min_segment_size
    ) {
        hslm::HSLMParameters params;
        params.omega = omega;
        params.theta = theta;
        params.step_eta = step_eta;
        params.n_states = n_states;
        params.min_segment_size = min_segment_size;
        hslm::HSLM segmenter(params);
        return segmenter.segment(log2_ratios, positions);
    },
    py::arg("log2_ratios"),
    py::arg("positions"),
    py::arg("omega") = 0.1,
    py::arg("theta") = 1e-5,
    py::arg("step_eta") = 200000.0,
    py::arg("n_states") = 21,
    py::arg("min_segment_size") = 1,
    R"pbdoc(
        Convenience function for single-call HSLM segmentation.

        Args:
            log2_ratios: Vector of log2 ratio values
            positions: Genomic positions
            omega: Variance partitioning factor (default: 0.1)
            theta: Base transition probability (default: 1e-5)
            step_eta: Distance step for eta (default: 200000)
            n_states: Number of hidden states (default: 21)
            min_segment_size: Minimum segment size to keep (default: 1)

        Returns:
            HSLMResult with segmentation results
    )pbdoc");

    // ==================== Common utilities ====================
    m.def("elnsum", &elnsum,
          py::arg("x"), py::arg("y"),
          R"pbdoc(
              Log-sum-exp: compute log(exp(x) + exp(y)) in a numerically stable way.

              Args:
                  x: First log-space value
                  y: Second log-space value

              Returns:
                  log(exp(x) + exp(y))
          )pbdoc");

    m.def("logsumexp", &logsumexp,
          py::arg("v"),
          R"pbdoc(
              Log-sum-exp for a vector: compute log(sum(exp(v[i]))) in a numerically stable way.

              Args:
                  v: Vector of log-space values

              Returns:
                  log(sum(exp(v[i])))
          )pbdoc");

    // ==================== FastCall Submodule ====================
    py::module_ fastcall_mod = m.def_submodule("fastcall", R"pbdoc(
        FastCall CNV classification algorithm.

        This module implements an Expectation-Maximization algorithm using a
        5-state Gaussian mixture model to classify segments into copy number states:
          - CN=0 (homozygous deletion)
          - CN=1 (heterozygous deletion)
          - CN=2 (normal/diploid)
          - CN=3 (single copy gain)
          - CN=4+ (amplification)

        Example:
            >>> from excavator2._excavator_core import fastcall
            >>> params = fastcall.FastCallParameters()
            >>> caller = fastcall.FastCall(params)
            >>> result = caller.call(segment_means)
            >>> for call in result.calls:
            ...     print(f"CN call: {call.cn_call}, prob: {call.probability:.3f}")
    )pbdoc");

    // FastCallParameters struct
    py::class_<fastcall::FastCallParameters>(fastcall_mod, "FastCallParameters", R"pbdoc(
        Parameters for the FastCall CNV classification algorithm.

        Attributes:
            cellularity (float): Tumor purity (0.0-1.0). Default: 1.0
            thrd (float): Lower threshold for normal state. Default: 0.5
            thru (float): Upper threshold for normal state. Default: 0.35
            min_exons (int): Minimum exons per segment. Default: 4
            max_iterations (int): Maximum EM iterations. Default: 1000
            convergence (float): Convergence threshold. Default: 1e-5
    )pbdoc")
        .def(py::init<>())
        .def_readwrite("cellularity", &fastcall::FastCallParameters::cellularity)
        .def_readwrite("thrd", &fastcall::FastCallParameters::thrd)
        .def_readwrite("thru", &fastcall::FastCallParameters::thru)
        .def_readwrite("min_exons", &fastcall::FastCallParameters::min_exons)
        .def_readwrite("max_iterations", &fastcall::FastCallParameters::max_iterations)
        .def_readwrite("convergence", &fastcall::FastCallParameters::convergence)
        .def("__repr__", [](const fastcall::FastCallParameters& p) {
            return "<FastCallParameters cellularity=" + std::to_string(p.cellularity) +
                   " thrd=" + std::to_string(p.thrd) +
                   " thru=" + std::to_string(p.thru) + ">";
        });

    // SegmentCall struct
    py::class_<fastcall::SegmentCall>(fastcall_mod, "SegmentCall", R"pbdoc(
        Copy number call for a single segment.

        Attributes:
            cn_call (int): Relative CN call (-2, -1, 0, +1, +2)
            absolute_cn (int): Absolute copy number (0, 1, 2, 3, 4+)
            probability (float): Posterior probability of the call
            state_index (int): Index of the most likely state (0-4)
            segment_mean (float): Mean log2 ratio of the segment
    )pbdoc")
        .def(py::init<>())
        .def_readonly("cn_call", &fastcall::SegmentCall::cn_call)
        .def_readonly("absolute_cn", &fastcall::SegmentCall::absolute_cn)
        .def_readonly("probability", &fastcall::SegmentCall::probability)
        .def_readonly("state_index", &fastcall::SegmentCall::state_index)
        .def_readonly("segment_mean", &fastcall::SegmentCall::segment_mean)
        .def("__repr__", [](const fastcall::SegmentCall& c) {
            return "<SegmentCall CN=" + std::to_string(c.absolute_cn) +
                   " call=" + std::to_string(c.cn_call) +
                   " prob=" + std::to_string(c.probability) + ">";
        });

    // FastCallResult struct
    py::class_<fastcall::FastCallResult>(fastcall_mod, "FastCallResult", R"pbdoc(
        Result of FastCall CNV classification.

        Attributes:
            calls (list[SegmentCall]): Calls for each segment
            state_means (list[float]): Fitted means for each state
            state_sds (list[float]): Fitted standard deviations
            state_priors (list[float]): Fitted prior probabilities
            iterations (int): Number of EM iterations
            converged (bool): Whether EM converged
            success (bool): Whether the algorithm succeeded
            error_message (str): Error message if failed
    )pbdoc")
        .def(py::init<>())
        .def_readonly("calls", &fastcall::FastCallResult::calls)
        .def_readonly("state_means", &fastcall::FastCallResult::state_means)
        .def_readonly("state_sds", &fastcall::FastCallResult::state_sds)
        .def_readonly("state_priors", &fastcall::FastCallResult::state_priors)
        .def_readonly("iterations", &fastcall::FastCallResult::iterations)
        .def_readonly("converged", &fastcall::FastCallResult::converged)
        .def_readonly("success", &fastcall::FastCallResult::success)
        .def_readonly("error_message", &fastcall::FastCallResult::error_message)
        .def("__repr__", [](const fastcall::FastCallResult& r) {
            return "<FastCallResult n_segments=" + std::to_string(r.calls.size()) +
                   " iterations=" + std::to_string(r.iterations) +
                   " converged=" + (r.converged ? "True" : "False") + ">";
        });

    // FastCall class
    py::class_<fastcall::FastCall>(fastcall_mod, "FastCall", R"pbdoc(
        FastCall CNV classification algorithm.

        Uses a 5-state Gaussian mixture model with EM to classify segments
        into copy number states.

        Example:
            >>> params = FastCallParameters()
            >>> params.cellularity = 0.8
            >>> caller = FastCall(params)
            >>> result = caller.call(segment_means)

        Args:
            params: FastCallParameters object with algorithm settings
    )pbdoc")
        .def(py::init<const fastcall::FastCallParameters&>(),
             py::arg("params") = fastcall::FastCallParameters())
        .def("call", &fastcall::FastCall::call,
             py::arg("segment_means"),
             py::arg("segment_sds") = std::vector<double>(),
             R"pbdoc(
                Call copy number states for segments.

                Args:
                    segment_means: Mean log2 ratios for each segment
                    segment_sds: Standard deviations for each segment (optional)

                Returns:
                    FastCallResult containing CN calls and probabilities
             )pbdoc")
        .def_property("params",
            &fastcall::FastCall::params,
            &fastcall::FastCall::set_params,
            "FastCall parameters");

    // Convenience function for single-call classification
    fastcall_mod.def("call", [](
        const std::vector<double>& segment_means,
        double cellularity,
        double thrd,
        double thru,
        int max_iterations
    ) {
        fastcall::FastCallParameters params;
        params.cellularity = cellularity;
        params.thrd = thrd;
        params.thru = thru;
        params.max_iterations = max_iterations;
        fastcall::FastCall caller(params);
        return caller.call(segment_means);
    },
    py::arg("segment_means"),
    py::arg("cellularity") = 1.0,
    py::arg("thrd") = 0.5,
    py::arg("thru") = 0.35,
    py::arg("max_iterations") = 1000,
    R"pbdoc(
        Convenience function for single-call FastCall classification.

        Args:
            segment_means: Vector of segment mean log2 ratios
            cellularity: Tumor purity (default: 1.0)
            thrd: Lower threshold for normal state (default: 0.5)
            thru: Upper threshold for normal state (default: 0.35)
            max_iterations: Maximum EM iterations (default: 1000)

        Returns:
            FastCallResult with CN calls
    )pbdoc");
}
