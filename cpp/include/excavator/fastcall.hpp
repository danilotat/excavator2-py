#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <string>
#include "common.hpp"

namespace excavator {
namespace fastcall {

/**
 * Number of copy number states in FastCall.
 * States correspond to:
 *   0: CN=0 (homozygous deletion)
 *   1: CN=1 (heterozygous deletion)
 *   2: CN=2 (normal/diploid)
 *   3: CN=3 (single copy gain)
 *   4: CN=4+ (amplification)
 */
constexpr int N_STATES = 5;

/**
 * Default means for each CN state (log2 ratio values).
 * From LibraryFastCall.R: muvec <- c(-3, -1, 0, 0.58, 1)
 */
constexpr std::array<double, N_STATES> DEFAULT_MEANS = {-3.0, -1.0, 0.0, 0.58, 1.0};

/**
 * Default priors for each CN state.
 * From LibraryFastCall.R: prior <- c(0.05, 0.1, 0.7, 0.1, 0.05)
 */
constexpr std::array<double, N_STATES> DEFAULT_PRIORS = {0.05, 0.1, 0.7, 0.1, 0.05};

/**
 * Parameters for the FastCall algorithm.
 */
struct FastCallParameters {
    double cellularity = 1.0;     // Tumor purity (0.0-1.0), affects expected log2 ratios
    double thrd = 0.5;            // Lower threshold for normal state (deletion boundary)
    double thru = 0.35;           // Upper threshold for normal state (duplication boundary)
    int min_exons = 4;            // Minimum exons per segment for calling
    int max_iterations = 1000;    // Maximum EM iterations
    double convergence = 1e-5;    // Convergence threshold for EM
};

/**
 * Result for a single segment's CNV call.
 */
struct SegmentCall {
    int cn_call;              // Copy number call: -2, -1, 0, +1, +2 (relative to diploid)
    int absolute_cn;          // Absolute copy number: 0, 1, 2, 3, 4+
    double probability;       // Posterior probability of the call
    int state_index;          // Index of the most likely state (0-4)
    double segment_mean;      // Mean log2 ratio of the segment
};

/**
 * Result of FastCall on multiple segments.
 */
struct FastCallResult {
    std::vector<SegmentCall> calls;        // Calls for each segment
    std::vector<double> state_means;       // Fitted means for each state
    std::vector<double> state_sds;         // Fitted standard deviations
    std::vector<double> state_priors;      // Fitted prior probabilities
    int iterations;                        // Number of EM iterations
    bool converged;                        // Whether EM converged
    bool success;                          // Whether the algorithm succeeded
    std::string error_message;             // Error message if failed
};

/**
 * FastCall algorithm for CNV classification.
 *
 * This is a port of the R implementation from LibraryFastCall.R.
 * FastCall uses a 5-state Gaussian mixture model with Expectation-Maximization
 * to classify segments into copy number states.
 *
 * The 5 states correspond to:
 *   - CN=0 (homozygous deletion): mean ~ -3.0
 *   - CN=1 (heterozygous deletion): mean ~ -1.0
 *   - CN=2 (normal): mean ~ 0.0
 *   - CN=3 (single copy gain): mean ~ 0.58
 *   - CN=4+ (amplification): mean ~ 1.0
 *
 * Reference: EXCAVATOR2 paper, D'Aurizio et al., NAR 2016
 */
class FastCall {
public:
    explicit FastCall(const FastCallParameters& params = FastCallParameters());

    /**
     * Call copy number states for segments.
     *
     * @param segment_means Mean log2 ratios for each segment
     * @param segment_sds   Standard deviations for each segment (optional, can be empty)
     * @return FastCallResult containing CN calls and probabilities
     */
    FastCallResult call(
        const std::vector<double>& segment_means,
        const std::vector<double>& segment_sds = {}
    );

    // Getters/setters
    const FastCallParameters& params() const { return params_; }
    void set_params(const FastCallParameters& params) { params_ = params; }

private:
    FastCallParameters params_;

    // EM algorithm state
    std::array<double, N_STATES> means_;
    std::array<double, N_STATES> sds_;
    std::array<double, N_STATES> priors_;

    // Boundary vectors for truncated Gaussians
    std::array<double, N_STATES> lower_bounds_;
    std::array<double, N_STATES> upper_bounds_;

    /**
     * Initialize starting conditions for EM.
     * Ported from StartCond in LibraryFastCall.R
     */
    void initialize_start_conditions(const std::vector<double>& data);

    /**
     * Run the EM algorithm.
     * Ported from EMFastCall in LibraryFastCall.R
     *
     * @param data Segment mean log2 ratios
     * @return Tuple of (iterations, converged)
     */
    std::pair<int, bool> run_em(const std::vector<double>& data);

    /**
     * E-step: compute responsibilities (posterior probabilities).
     * Ported from EStep in LibraryFastCall.R
     *
     * @param data Segment means
     * @return Matrix of responsibilities [n_segments x N_STATES]
     */
    std::vector<std::vector<double>> e_step(const std::vector<double>& data);

    /**
     * M-step: update parameters from responsibilities.
     * Ported from MStep in LibraryFastCall.R
     *
     * @param data Segment means
     * @param tau Responsibilities from E-step
     */
    void m_step(const std::vector<double>& data,
                const std::vector<std::vector<double>>& tau);

    /**
     * Compute posterior probabilities for state assignment.
     * Ported from PosteriorP in LibraryFastCall.R
     *
     * @param data Segment means
     * @return Matrix of posterior probabilities [n_segments x N_STATES]
     */
    std::vector<std::vector<double>> compute_posteriors(const std::vector<double>& data);

    /**
     * Truncated Gaussian PDF.
     * Ported from gfct in LibraryFastCall.R
     *
     * @param x Value to evaluate
     * @param mean Gaussian mean
     * @param sd Gaussian standard deviation
     * @param lower Lower truncation bound
     * @param upper Upper truncation bound
     * @return Truncated Gaussian PDF value
     */
    static double truncated_gaussian_pdf(
        double x, double mean, double sd, double lower, double upper
    );

    /**
     * Assign CN labels from posterior probabilities.
     * Ported from LabelAss in LibraryFastCall.R
     *
     * @param posteriors Posterior probability matrix
     * @param data Original segment means
     * @return Vector of SegmentCall objects
     */
    std::vector<SegmentCall> assign_labels(
        const std::vector<std::vector<double>>& posteriors,
        const std::vector<double>& data
    );

    /**
     * Compute log-likelihood for convergence check.
     */
    double compute_log_likelihood(
        const std::vector<double>& data,
        const std::vector<std::vector<double>>& posteriors
    );
};

} // namespace fastcall
} // namespace excavator
