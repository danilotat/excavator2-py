#include "excavator/hslm.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <limits>

#ifdef EXCAVATOR_HAS_OPENMP
#include <omp.h>
#endif

namespace excavator {
namespace hslm {

HSLM::HSLM(const HSLMParameters& params) : params_(params) {}

HSLMResult HSLM::segment(
    const std::vector<double>& log2_ratios,
    const std::vector<int64_t>& positions
) {
    // Wrap single sequence as matrix for unified processing
    std::vector<std::vector<double>> data_matrix = {log2_ratios};
    return segment_multi(data_matrix, positions);
}

HSLMResult HSLM::segment_multi(
    const std::vector<std::vector<double>>& data_matrix,
    const std::vector<int64_t>& positions
) {
    HSLMResult result;
    result.success = false;

    // Input validation
    if (data_matrix.empty() || data_matrix[0].empty()) {
        result.error_message = "Empty data matrix";
        return result;
    }

    size_t n_positions = data_matrix[0].size();
    if (positions.size() != n_positions) {
        result.error_message = "Position vector length does not match data";
        return result;
    }

    if (n_positions < 2) {
        result.error_message = "Need at least 2 positions for segmentation";
        return result;
    }

    // Estimate parameters from data
    EstimatedParameters est_params = estimate_parameters(data_matrix);

    // Estimate state means
    std::vector<std::vector<double>> muk = estimate_state_means(data_matrix);
    int n_states = static_cast<int>(muk[0].size());

    // Compute distance-dependent transition probabilities
    std::vector<double> eta_vec = compute_eta_vector(positions);
    int n_cov = static_cast<int>(eta_vec.size());

    // Initialize uniform log-probabilities for states
    double log_prob = std::log(1.0 / n_states);
    std::vector<double> etav(n_states, log_prob);

    // Compute transition and emission matrices
    std::vector<std::vector<double>> P(n_states, std::vector<double>(n_states * n_cov, 0.0));
    std::vector<std::vector<double>> emission(n_states, std::vector<double>(n_positions, 0.0));

    compute_transition_emission(
        muk, est_params.mi, eta_vec, data_matrix,
        est_params.smu, est_params.sepsilon,
        P, emission
    );

    // Run Viterbi algorithm
    std::vector<int> state_path = viterbi(etav, P, emission);

    // Extract breakpoints from state path
    std::vector<int> breakpoints = extract_breakpoints(state_path);

    // Filter short segments if requested
    if (params_.min_segment_size > 1) {
        breakpoints = filter_segments(breakpoints, params_.min_segment_size);
    }

    // Compute segment means for first sequence
    std::vector<double> segment_means = compute_segment_means(data_matrix[0], breakpoints);

    // Store state values
    std::vector<double> state_values(n_states);
    for (int i = 0; i < n_states; ++i) {
        state_values[i] = muk[0][i];
    }

    // Populate result
    result.breakpoints = breakpoints;
    result.segment_means = segment_means;
    result.state_path = state_path;
    result.state_values = state_values;
    result.n_segments = static_cast<int>(breakpoints.size()) - 1;
    result.success = true;

    return result;
}

EstimatedParameters HSLM::estimate_parameters(
    const std::vector<std::vector<double>>& data_matrix
) {
    // Ported from ParamEstSeq in LibraryJSLMIn.R
    EstimatedParameters params;
    size_t n_exp = data_matrix.size();

    params.mi.resize(n_exp);
    params.smu.resize(n_exp);
    params.sepsilon.resize(n_exp);

    for (size_t i = 0; i < n_exp; ++i) {
        const auto& seq = data_matrix[i];

        // Compute trimmed variance (1-99 percentile)
        double var = trimmed_variance(seq, 0.01, 0.99);

        // Mean is always 0 for log2 ratios
        params.mi[i] = 0.0;

        // Split variance between state and noise based on omega
        params.smu[i] = std::sqrt(params_.omega * var);
        params.sepsilon[i] = std::sqrt((1.0 - params_.omega) * var);

        // Ensure minimum values for numerical stability
        if (params.smu[i] < EPSILON) params.smu[i] = 0.01;
        if (params.sepsilon[i] < EPSILON) params.sepsilon[i] = 0.01;
    }

    return params;
}

std::vector<std::vector<double>> HSLM::estimate_state_means(
    const std::vector<std::vector<double>>& data_matrix
) {
    // Ported from MukEst in LibraryJSLMIn.R
    // For single sequence or default: use grid from -1 to 1 by 0.1
    size_t n_exp = data_matrix.size();

    // Create default state means: -1, -0.9, ..., 0, ..., 0.9, 1.0
    std::vector<double> default_states;
    for (double v = -1.0; v <= 1.0 + 1e-9; v += 0.1) {
        default_states.push_back(v);
    }

    // Each sequence gets the same state means in single-sample mode
    std::vector<std::vector<double>> muk(n_exp, default_states);

    return muk;
}

std::vector<double> HSLM::compute_eta_vector(
    const std::vector<int64_t>& positions
) {
    // Compute distance-dependent transition probabilities
    // From JointSegIn in LibraryJSLMIn.R:
    // CovPos <- diff(Pos)
    // CovPosNorm <- CovPos / stepeta
    // etavec <- eta + (1 - eta) * exp(log(eta) / CovPosNorm)

    size_t n = positions.size();
    std::vector<double> eta_vec(n - 1);

    double theta = params_.theta;
    double step_eta = params_.step_eta;
    double log_theta = std::log(theta);

    for (size_t i = 0; i < n - 1; ++i) {
        double cov_pos = static_cast<double>(positions[i + 1] - positions[i]);
        double cov_pos_norm = cov_pos / step_eta;

        // eta + (1 - eta) * exp(log(eta) / cov_pos_norm)
        // This gives higher transition probability for larger distances
        if (cov_pos_norm > 0) {
            eta_vec[i] = theta + (1.0 - theta) * std::exp(log_theta / cov_pos_norm);
        } else {
            eta_vec[i] = theta;  // Same position: use base rate
        }

        // Clamp to valid probability range
        eta_vec[i] = std::max(1e-10, std::min(1.0 - 1e-10, eta_vec[i]));
    }

    return eta_vec;
}

void HSLM::compute_transition_emission(
    const std::vector<std::vector<double>>& muk,
    const std::vector<double>& mi,
    const std::vector<double>& eta_vec,
    const std::vector<std::vector<double>>& data,
    const std::vector<double>& smu,
    const std::vector<double>& sepsilon,
    std::vector<std::vector<double>>& P,
    std::vector<std::vector<double>>& emission
) {
    // Ported from TRANSEMISI subroutine in FastJointSLMLibraryI.f

    size_t n_exp = data.size();
    size_t T = data[0].size();
    int K = static_cast<int>(muk[0].size());
    int n_cov = static_cast<int>(eta_vec.size());

    // Step 1: Compute GVECT - log probability of each state given state means
    // GSUP = sum over sequences of -(muk[j,i] - mi[j])^2 / (2 * smu[j]^2)
    std::vector<double> gvect(K);
    for (int i = 0; i < K; ++i) {
        double gsup = 0.0;
        for (size_t j = 0; j < n_exp; ++j) {
            double diff = muk[j][i] - mi[j];
            gsup += -(diff * diff) / (2.0 * smu[j] * smu[j]);
        }
        gvect[i] = gsup;
    }

    // Normalize GVECT using log-sum-exp
    double norm = gvect[0];
    for (int j = 1; j < K; ++j) {
        norm = elnsum(norm, gvect[j]);
    }

    // G matrix: G[i,j] = gvect[j] - norm (same column for all rows)
    std::vector<std::vector<double>> G(K, std::vector<double>(K));
    for (int i = 0; i < K; ++i) {
        for (int j = 0; j < K; ++j) {
            G[i][j] = gvect[j] - norm;
        }
    }

    // Step 2: Compute transition matrix P
    // P has dimensions (K, K * n_cov) - one K x K block per covariate
    // P[j, k + countp] = log(eta) + G[j,k]  if j != k
    // P[j, k + countp] = log(1-eta) + log(1 + exp(log(eta) + G[j,k] - log(1-eta)))  if j == k

    for (int cov_idx = 0; cov_idx < n_cov; ++cov_idx) {
        double eta = eta_vec[cov_idx];
        double log_eta = std::log(eta);
        double log_one_minus_eta = std::log(1.0 - eta);
        int countp = cov_idx * K;

        for (int j = 0; j < K; ++j) {
            for (int k = 0; k < K; ++k) {
                if (j == k) {
                    // Diagonal: elnsum(log(1-eta), log(eta) + G[j,k])
                    P[j][k + countp] = elnsum(log_one_minus_eta, log_eta + G[j][k]);
                } else {
                    // Off-diagonal: log(eta) + G[j,k]
                    P[j][k + countp] = log_eta + G[j][k];
                }
            }
        }
    }

    // Step 3: Compute emission matrix
    // emission[j,k] = sum over sequences of log(N(data[i,k] | muk[i,j], sepsilon[i]))
    // Initialize to zero (will accumulate log-likelihoods)

#ifdef EXCAVATOR_HAS_OPENMP
    #pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < K; ++j) {
        for (size_t k = 0; k < T; ++k) {
            double em = 0.0;
            for (size_t i = 0; i < n_exp; ++i) {
                em += log_normal_pdf(data[i][k], muk[i][j], sepsilon[i]);
            }
            emission[j][k] = em;
        }
    }
}

std::vector<int> HSLM::viterbi(
    const std::vector<double>& etav,
    const std::vector<std::vector<double>>& P,
    const std::vector<std::vector<double>>& emission
) {
    // Ported from BIOVITERBII subroutine in FastJointSLMLibraryI.f

    int K = static_cast<int>(etav.size());
    size_t T = emission[0].size();

    // Delta: forward probabilities
    std::vector<std::vector<double>> delta(K, std::vector<double>(T));

    // Psi: backtracking indices
    std::vector<std::vector<int>> psi(K, std::vector<int>(T, 0));

    // Initialize: delta[i,0] = etav[i] + emission[i,0]
    for (int i = 0; i < K; ++i) {
        delta[i][0] = etav[i] + emission[i][0];
        psi[i][0] = 0;
    }

    // Forward pass
    int countp = 0;
    for (size_t t = 1; t < T; ++t) {
        for (int j = 0; j < K; ++j) {
            double num_max = delta[0][t - 1] + P[0][j + countp];
            int ind = 0;

            for (int k = 1; k < K; ++k) {
                double val = delta[k][t - 1] + P[k][j + countp];
                if (val > num_max) {
                    num_max = val;
                    ind = k;
                }
            }

            psi[j][t] = ind;
            delta[j][t] = num_max + emission[j][t];
        }
        countp += K;
    }

    // Find best final state
    double num_max = delta[0][T - 1];
    int ind = 0;
    for (int k = 1; k < K; ++k) {
        if (delta[k][T - 1] > num_max) {
            num_max = delta[k][T - 1];
            ind = k;
        }
    }

    // Backtrack to find path
    std::vector<int> path(T);
    path[T - 1] = ind;
    for (int t = static_cast<int>(T) - 2; t >= 0; --t) {
        path[t] = psi[path[t + 1]][t + 1];
    }

    return path;
}

std::vector<int> HSLM::extract_breakpoints(const std::vector<int>& state_path) {
    // Ported from SortState in LibraryJSLMIn.R

    std::vector<int> breakpoints;
    breakpoints.push_back(0);  // Always start with 0

    for (size_t i = 0; i < state_path.size() - 1; ++i) {
        if (state_path[i] != state_path[i + 1]) {
            breakpoints.push_back(static_cast<int>(i + 1));  // Breakpoint at next position
        }
    }

    breakpoints.push_back(static_cast<int>(state_path.size()));  // End of data

    return breakpoints;
}

std::vector<double> HSLM::compute_segment_means(
    const std::vector<double>& data,
    const std::vector<int>& breakpoints
) {
    // Ported from SegResults in LibraryJSLMIn.R

    std::vector<double> segment_means;
    segment_means.reserve(breakpoints.size() - 1);

    for (size_t i = 0; i < breakpoints.size() - 1; ++i) {
        int start = breakpoints[i];
        int end = breakpoints[i + 1];

        // Extract segment data
        std::vector<double> segment(data.begin() + start, data.begin() + end);

        // Compute median of segment
        segment_means.push_back(median(segment));
    }

    return segment_means;
}

std::vector<int> HSLM::filter_segments(
    const std::vector<int>& breakpoints,
    int min_size
) {
    // Ported from FilterSeg in LibraryJSLMIn.R
    //
    // R logic: for each segment i with length <= FW, remove the breakpoint
    // at the START of that segment. Special case: if segment 0 is short,
    // we can't remove bp[0], so remove bp[1] instead.

    if (breakpoints.size() <= 2) {
        return breakpoints;  // Can't filter if only one segment
    }

    // Identify breakpoints to remove
    std::vector<bool> remove(breakpoints.size(), false);

    for (size_t seg = 0; seg < breakpoints.size() - 1; ++seg) {
        int seg_length = breakpoints[seg + 1] - breakpoints[seg];
        if (seg_length <= min_size) {
            // Segment is short, mark its starting breakpoint for removal
            // But if seg == 0, we can't remove bp[0], so mark bp[1] instead
            if (seg == 0) {
                remove[1] = true;
            } else {
                remove[seg] = true;
            }
        }
    }

    // Build filtered result, always keeping first and last
    std::vector<int> filtered;
    filtered.push_back(breakpoints[0]);  // Always keep start

    for (size_t i = 1; i < breakpoints.size() - 1; ++i) {
        if (!remove[i]) {
            filtered.push_back(breakpoints[i]);
        }
    }

    filtered.push_back(breakpoints.back());  // Always keep end

    return filtered;
}

PreEstimatedParams HSLM::estimate_params(
    const std::vector<std::vector<double>>& data_matrix
) {
    PreEstimatedParams result;

    if (data_matrix.empty() || data_matrix[0].empty()) {
        result.valid = false;
        return result;
    }

    // Estimate mi, smu, sepsilon
    EstimatedParameters est = estimate_parameters(data_matrix);
    result.mi = est.mi;
    result.smu = est.smu;
    result.sepsilon = est.sepsilon;

    // Estimate state means
    result.muk = estimate_state_means(data_matrix);

    result.valid = true;
    return result;
}

HSLMResult HSLM::segment_with_params(
    const std::vector<double>& log2_ratios,
    const std::vector<int64_t>& positions,
    const PreEstimatedParams& params
) {
    // Wrap single sequence as matrix for unified processing
    std::vector<std::vector<double>> data_matrix = {log2_ratios};

    HSLMResult result;
    result.success = false;

    // Input validation
    if (data_matrix.empty() || data_matrix[0].empty()) {
        result.error_message = "Empty data matrix";
        return result;
    }

    size_t n_positions = data_matrix[0].size();
    if (positions.size() != n_positions) {
        result.error_message = "Position vector length does not match data";
        return result;
    }

    if (n_positions < 2) {
        result.error_message = "Need at least 2 positions for segmentation";
        return result;
    }

    if (!params.valid) {
        result.error_message = "Invalid pre-estimated parameters";
        return result;
    }

    // Use pre-estimated parameters instead of computing from this data
    // We need to adapt the parameters to single-sample mode:
    // Extract first element for single-sample
    EstimatedParameters est_params;
    est_params.mi = {params.mi.empty() ? 0.0 : params.mi[0]};
    est_params.smu = {params.smu.empty() ? 0.01 : params.smu[0]};
    est_params.sepsilon = {params.sepsilon.empty() ? 0.01 : params.sepsilon[0]};

    // Use the state means from pre-estimated params
    // For single sample, we need a 1 x n_states matrix
    std::vector<std::vector<double>> muk;
    if (!params.muk.empty() && !params.muk[0].empty()) {
        muk = {params.muk[0]};  // Take first row for single sample
    } else {
        // Fallback to default state means
        muk = estimate_state_means(data_matrix);
    }

    int n_states = static_cast<int>(muk[0].size());

    // Compute distance-dependent transition probabilities
    std::vector<double> eta_vec = compute_eta_vector(positions);
    int n_cov = static_cast<int>(eta_vec.size());

    // Initialize uniform log-probabilities for states
    double log_prob = std::log(1.0 / n_states);
    std::vector<double> etav(n_states, log_prob);

    // Compute transition and emission matrices
    std::vector<std::vector<double>> P(n_states, std::vector<double>(n_states * n_cov, 0.0));
    std::vector<std::vector<double>> emission(n_states, std::vector<double>(n_positions, 0.0));

    compute_transition_emission(
        muk, est_params.mi, eta_vec, data_matrix,
        est_params.smu, est_params.sepsilon,
        P, emission
    );

    // Run Viterbi algorithm
    std::vector<int> state_path = viterbi(etav, P, emission);

    // Extract breakpoints from state path
    std::vector<int> breakpoints = extract_breakpoints(state_path);

    // Filter short segments if requested
    if (params_.min_segment_size > 1) {
        breakpoints = filter_segments(breakpoints, params_.min_segment_size);
    }

    // Compute segment means for first sequence
    std::vector<double> segment_means = compute_segment_means(data_matrix[0], breakpoints);

    // Store state values
    std::vector<double> state_values(n_states);
    for (int i = 0; i < n_states; ++i) {
        state_values[i] = muk[0][i];
    }

    // Populate result
    result.breakpoints = breakpoints;
    result.segment_means = segment_means;
    result.state_path = state_path;
    result.state_values = state_values;
    result.n_segments = static_cast<int>(breakpoints.size()) - 1;
    result.success = true;

    return result;
}

} // namespace hslm
} // namespace excavator
