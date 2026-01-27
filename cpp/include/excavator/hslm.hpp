#pragma once

#include <vector>
#include <cstdint>
#include <string>
#include "common.hpp"

namespace excavator {
namespace hslm {

/**
 * Parameters for the Heterogeneous Shifting Level Model (HSLM) algorithm.
 *
 * The HSLM is an inhomogeneous Hidden Markov Model where transition probabilities
 * depend on genomic distance between consecutive observations.
 */
struct HSLMParameters {
    double omega = 0.1;         // Variance partitioning: fraction of variance attributed to states
    double theta = 1e-5;        // Base transition probability (eta)
    double d_norm = 1e6;        // Distance normalization factor
    double step_eta = 200000.0; // Distance step for eta computation (in base pairs)
    int n_states = 21;          // Number of hidden states (default: 21 for range -1 to 1 by 0.1)
    int min_segment_size = 1;   // Minimum segment size to keep (filter window)
};

/**
 * Result of HSLM segmentation.
 */
struct HSLMResult {
    std::vector<int> breakpoints;           // Indices where segments change (0-indexed)
    std::vector<double> segment_means;      // Mean value for each segment
    std::vector<int> state_path;            // State sequence from Viterbi algorithm
    std::vector<double> state_values;       // The actual values corresponding to each state
    int n_segments;                         // Number of segments detected
    bool success;                           // Whether the algorithm succeeded
    std::string error_message;              // Error message if failed
};

/**
 * Estimated parameters from data.
 */
struct EstimatedParameters {
    std::vector<double> mi;        // Mean for each sequence (typically 0)
    std::vector<double> smu;       // State standard deviation
    std::vector<double> sepsilon;  // Noise standard deviation
};

/**
 * HSLM (Heterogeneous Shifting Level Model) segmentation algorithm.
 *
 * This is a port of the Fortran implementation from FastJointSLMLibraryI.f
 * and the R wrapper from LibraryJSLMIn.R.
 *
 * The algorithm detects change points (breakpoints) in log2-ratio data
 * using an inhomogeneous HMM where transition probabilities depend on
 * genomic distance.
 *
 * Reference: EXCAVATOR2 paper, D'Aurizio et al., NAR 2016
 */
class HSLM {
public:
    explicit HSLM(const HSLMParameters& params = HSLMParameters());

    /**
     * Run segmentation on a single chromosome.
     *
     * @param log2_ratios Vector of log2 ratio values (test/control)
     * @param positions   Genomic positions corresponding to each log2 ratio
     * @return HSLMResult containing breakpoints and segment information
     */
    HSLMResult segment(
        const std::vector<double>& log2_ratios,
        const std::vector<int64_t>& positions
    );

    /**
     * Run segmentation on multiple sequences (multi-sample mode).
     *
     * @param data_matrix Matrix of log2 ratios, shape (n_sequences, n_positions)
     * @param positions   Genomic positions
     * @return HSLMResult containing breakpoints
     */
    HSLMResult segment_multi(
        const std::vector<std::vector<double>>& data_matrix,
        const std::vector<int64_t>& positions
    );

    // Getters/setters for parameters
    const HSLMParameters& params() const { return params_; }
    void set_params(const HSLMParameters& params) { params_ = params; }

private:
    HSLMParameters params_;

    /**
     * Estimate initial parameters from data.
     * Ported from ParamEstSeq in LibraryJSLMIn.R
     */
    EstimatedParameters estimate_parameters(
        const std::vector<std::vector<double>>& data_matrix
    );

    /**
     * Estimate state means (muk matrix).
     * Ported from MukEst in LibraryJSLMIn.R
     */
    std::vector<std::vector<double>> estimate_state_means(
        const std::vector<std::vector<double>>& data_matrix
    );

    /**
     * Compute distance-dependent transition probabilities.
     * eta_vec[i] = theta + (1 - theta) * exp(log(theta) / cov_pos_norm[i])
     */
    std::vector<double> compute_eta_vector(
        const std::vector<int64_t>& positions
    );

    /**
     * Compute HMM transition and emission matrices.
     * Ported from TRANSEMISI subroutine in FastJointSLMLibraryI.f
     *
     * @param muk       State means matrix (n_sequences x n_states)
     * @param mi        Data means for each sequence
     * @param eta_vec   Distance-dependent transition probabilities
     * @param data      Data matrix (n_sequences x n_positions)
     * @param smu       State standard deviations
     * @param sepsilon  Noise standard deviations
     * @param P         Output: transition matrix (n_states x n_states * n_covariates)
     * @param emission  Output: emission matrix (n_states x n_positions)
     */
    void compute_transition_emission(
        const std::vector<std::vector<double>>& muk,
        const std::vector<double>& mi,
        const std::vector<double>& eta_vec,
        const std::vector<std::vector<double>>& data,
        const std::vector<double>& smu,
        const std::vector<double>& sepsilon,
        std::vector<std::vector<double>>& P,
        std::vector<std::vector<double>>& emission
    );

    /**
     * Viterbi algorithm to find most likely state sequence.
     * Ported from BIOVITERBII subroutine in FastJointSLMLibraryI.f
     *
     * @param etav      Initial state log-probabilities (uniform)
     * @param P         Transition matrix (time-varying)
     * @param emission  Emission matrix
     * @return          Most likely state path (0-indexed)
     */
    std::vector<int> viterbi(
        const std::vector<double>& etav,
        const std::vector<std::vector<double>>& P,
        const std::vector<std::vector<double>>& emission
    );

    /**
     * Extract breakpoints from state path.
     * Ported from SortState in LibraryJSLMIn.R
     */
    std::vector<int> extract_breakpoints(const std::vector<int>& state_path);

    /**
     * Compute segment means from breakpoints.
     * Ported from SegResults in LibraryJSLMIn.R
     */
    std::vector<double> compute_segment_means(
        const std::vector<double>& data,
        const std::vector<int>& breakpoints
    );

    /**
     * Filter short segments.
     * Ported from FilterSeg in LibraryJSLMIn.R
     */
    std::vector<int> filter_segments(
        const std::vector<int>& breakpoints,
        int min_size
    );
};

} // namespace hslm
} // namespace excavator
