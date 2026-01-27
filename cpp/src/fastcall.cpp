#include "excavator/fastcall.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <limits>

namespace excavator {
namespace fastcall {

// Standard normal CDF using error function
static double normal_cdf(double x, double mean, double sd) {
    return 0.5 * (1.0 + std::erf((x - mean) / (sd * std::sqrt(2.0))));
}

// Standard normal PDF
static double normal_pdf(double x, double mean, double sd) {
    double z = (x - mean) / sd;
    return std::exp(-0.5 * z * z) / (sd * std::sqrt(2.0 * PI));
}

FastCall::FastCall(const FastCallParameters& params)
    : params_(params)
{
    // Initialize with default values
    for (int i = 0; i < N_STATES; ++i) {
        means_[i] = DEFAULT_MEANS[i];
        sds_[i] = 0.01;  // Start with small SD
        priors_[i] = DEFAULT_PRIORS[i];
    }
}

FastCallResult FastCall::call(
    const std::vector<double>& segment_means,
    const std::vector<double>& segment_sds
) {
    FastCallResult result;
    result.success = false;

    // Input validation
    if (segment_means.empty()) {
        result.error_message = "Empty segment data";
        return result;
    }

    // Initialize starting conditions
    initialize_start_conditions(segment_means);

    // Run EM algorithm
    auto [iterations, converged] = run_em(segment_means);
    result.iterations = iterations;
    result.converged = converged;

    // Compute final posteriors
    auto posteriors = compute_posteriors(segment_means);

    // Assign labels
    result.calls = assign_labels(posteriors, segment_means);

    // Store fitted parameters
    result.state_means.assign(means_.begin(), means_.end());
    result.state_sds.assign(sds_.begin(), sds_.end());
    result.state_priors.assign(priors_.begin(), priors_.end());

    result.success = true;
    return result;
}

void FastCall::initialize_start_conditions(const std::vector<double>& data) {
    // Ported from StartCond in LibraryFastCall.R
    // muvec <- c(-3, -1, 0, 0.58, 1)
    // lvec <- c(-50, -1.5, -thrd, thru, 0.9)
    // uvec <- c(-1.5, -thrd, thru, 0.9, 50)

    double thrd = params_.thrd;
    double thru = params_.thru;

    // Set means (fixed)
    means_ = DEFAULT_MEANS;

    // Set boundaries for each state
    lower_bounds_ = {-50.0, -1.5, -thrd, thru, 0.9};
    upper_bounds_ = {-1.5, -thrd, thru, 0.9, 50.0};

    // Initialize standard deviations from data in each region
    sds_.fill(0.01);  // Default small SD

    for (int i = 0; i < N_STATES; ++i) {
        double lower = lower_bounds_[i];
        double upper = upper_bounds_[i];

        // Find data points in this region
        std::vector<double> region_data;
        for (double x : data) {
            if (x >= lower && x <= upper) {
                region_data.push_back(x);
            }
        }

        // Compute SD if we have enough data points
        if (region_data.size() > 1) {
            double sum = std::accumulate(region_data.begin(), region_data.end(), 0.0);
            double mean = sum / region_data.size();
            double sq_sum = 0.0;
            for (double x : region_data) {
                sq_sum += (x - mean) * (x - mean);
            }
            double sd = std::sqrt(sq_sum / (region_data.size() - 1));
            if (sd > 0.001) {
                sds_[i] = sd;
            }
        }
    }

    // Ensure minimum SD values
    for (int i = 0; i < N_STATES; ++i) {
        if (sds_[i] < 0.001) {
            sds_[i] = 0.001;
        }
    }

    // Reset priors to defaults
    priors_ = DEFAULT_PRIORS;
}

std::pair<int, bool> FastCall::run_em(const std::vector<double>& data) {
    // Ported from EMFastCall in LibraryFastCall.R

    double threshold = params_.convergence;
    int max_iter = params_.max_iterations;

    // Compute initial log-likelihood
    auto posteriors = compute_posteriors(data);
    double log_likelihood_old = compute_log_likelihood(data, posteriors);

    int iter;
    bool converged = false;

    for (iter = 0; iter < max_iter; ++iter) {
        // E-step
        auto tau = e_step(data);

        // M-step
        m_step(data, tau);

        // Compute new log-likelihood
        posteriors = compute_posteriors(data);
        double log_likelihood_new = compute_log_likelihood(data, posteriors);

        // Check convergence
        if (std::abs(log_likelihood_new - log_likelihood_old) < threshold) {
            converged = true;
            break;
        }

        log_likelihood_old = log_likelihood_new;
    }

    return {iter + 1, converged};
}

std::vector<std::vector<double>> FastCall::e_step(const std::vector<double>& data) {
    // Ported from EStep in LibraryFastCall.R
    // Uses truncated Gaussian for each state

    size_t n = data.size();
    std::vector<std::vector<double>> tau(n, std::vector<double>(N_STATES, 0.0));

    for (size_t i = 0; i < n; ++i) {
        double x = data[i];
        std::vector<double> probs(N_STATES);
        double sum_prob = 0.0;

        for (int j = 0; j < N_STATES; ++j) {
            double lower = lower_bounds_[j];
            double upper = upper_bounds_[j];

            // Check if data point is in this state's range
            bool in_range = (x >= lower && x <= upper);

            if (in_range) {
                double pdf = truncated_gaussian_pdf(x, means_[j], sds_[j], lower, upper);
                // Handle infinity
                if (std::isinf(pdf)) {
                    pdf = 100.0;  // Cap at 100 as in R code
                }
                probs[j] = priors_[j] * pdf;
            } else {
                probs[j] = 0.0;
            }

            sum_prob += probs[j];
        }

        // Handle case where no state matches (all zeros)
        if (sum_prob == 0.0) {
            // Assign to closest mean state
            int closest = 0;
            double min_dist = std::abs(means_[0] - x);
            for (int j = 1; j < N_STATES; ++j) {
                double dist = std::abs(means_[j] - x);
                if (dist < min_dist) {
                    min_dist = dist;
                    closest = j;
                }
            }
            tau[i][closest] = 1.0;
        } else {
            // Normalize
            for (int j = 0; j < N_STATES; ++j) {
                tau[i][j] = probs[j] / sum_prob;
            }
        }
    }

    return tau;
}

void FastCall::m_step(
    const std::vector<double>& data,
    const std::vector<std::vector<double>>& tau
) {
    // Ported from MStep in LibraryFastCall.R
    // Note: means are fixed, only update SD and priors

    size_t n = data.size();

    // Compute column sums (total responsibility for each state)
    std::array<double, N_STATES> col_sums;
    col_sums.fill(0.0);
    for (size_t i = 0; i < n; ++i) {
        for (int j = 0; j < N_STATES; ++j) {
            col_sums[j] += tau[i][j];
        }
    }

    for (int j = 0; j < N_STATES; ++j) {
        if (col_sums[j] > 0) {
            // Update standard deviation (mean is fixed)
            double weighted_sq_sum = 0.0;
            for (size_t i = 0; i < n; ++i) {
                double diff = data[i] - means_[j];
                weighted_sq_sum += tau[i][j] * diff * diff;
            }
            double new_sd = std::sqrt(weighted_sq_sum / col_sums[j]);

            // Apply minimum SD threshold
            if (new_sd < 1e-100) {
                // Keep old SD
            } else {
                sds_[j] = new_sd;
            }

            // Update prior
            priors_[j] = col_sums[j] / n;
        } else {
            // No data assigned to this state
            priors_[j] = 1e-6;
        }
    }

    // Normalize priors
    double prior_sum = std::accumulate(priors_.begin(), priors_.end(), 0.0);
    for (int j = 0; j < N_STATES; ++j) {
        priors_[j] /= prior_sum;
    }
}

std::vector<std::vector<double>> FastCall::compute_posteriors(
    const std::vector<double>& data
) {
    // Ported from PosteriorP in LibraryFastCall.R
    // Uses standard (non-truncated) Gaussian for posterior computation

    size_t n = data.size();
    std::vector<std::vector<double>> posteriors(n, std::vector<double>(N_STATES, 0.0));

    for (size_t i = 0; i < n; ++i) {
        double x = data[i];
        std::vector<double> probs(N_STATES);
        double sum_prob = 0.0;

        for (int j = 0; j < N_STATES; ++j) {
            double pdf = normal_pdf(x, means_[j], sds_[j]);
            probs[j] = priors_[j] * pdf;
            sum_prob += probs[j];
        }

        // Handle case where all probabilities are zero
        if (sum_prob == 0.0) {
            // Assign to closest mean
            int closest = 0;
            double min_dist = std::abs(means_[0] - x);
            for (int j = 1; j < N_STATES; ++j) {
                double dist = std::abs(means_[j] - x);
                if (dist < min_dist) {
                    min_dist = dist;
                    closest = j;
                }
            }
            posteriors[i][closest] = 1.0;
        } else {
            for (int j = 0; j < N_STATES; ++j) {
                posteriors[i][j] = probs[j] / sum_prob;
            }
        }
    }

    return posteriors;
}

double FastCall::truncated_gaussian_pdf(
    double x, double mean, double sd, double lower, double upper
) {
    // Ported from gfct in LibraryFastCall.R:
    // gfct <- function(x, moy, sdev, l, u) {
    //     (dnorm(x, mean=moy, sd=sdev) * (x<=u) * (x>=l)) /
    //     (pnorm(u, mean=moy, sd=sdev) - pnorm(l, mean=moy, sd=sdev))
    // }

    // Check if x is in bounds
    if (x < lower || x > upper) {
        return 0.0;
    }

    // Compute normalization constant
    double cdf_upper = normal_cdf(upper, mean, sd);
    double cdf_lower = normal_cdf(lower, mean, sd);
    double normalization = cdf_upper - cdf_lower;

    if (normalization < 1e-300) {
        // Truncation region has negligible probability
        return 0.0;
    }

    return normal_pdf(x, mean, sd) / normalization;
}

std::vector<SegmentCall> FastCall::assign_labels(
    const std::vector<std::vector<double>>& posteriors,
    const std::vector<double>& data
) {
    // Ported from LabelAss in LibraryFastCall.R
    // CallResults[indcall==1] <- -2  (CN=0, homozygous deletion)
    // CallResults[indcall==2] <- -1  (CN=1, heterozygous deletion)
    // CallResults[indcall==3] <-  0  (CN=2, normal)
    // CallResults[indcall==4] <-  1  (CN=3, single copy gain)
    // CallResults[indcall==5] <-  2  (CN=4+, amplification)

    static const std::array<int, N_STATES> CN_CALLS = {-2, -1, 0, 1, 2};
    static const std::array<int, N_STATES> ABSOLUTE_CN = {0, 1, 2, 3, 4};

    size_t n = posteriors.size();
    std::vector<SegmentCall> calls(n);

    for (size_t i = 0; i < n; ++i) {
        // Find state with maximum posterior probability
        int max_state = 0;
        double max_prob = posteriors[i][0];

        for (int j = 1; j < N_STATES; ++j) {
            if (posteriors[i][j] > max_prob) {
                max_prob = posteriors[i][j];
                max_state = j;
            }
        }

        calls[i].cn_call = CN_CALLS[max_state];
        calls[i].absolute_cn = ABSOLUTE_CN[max_state];
        calls[i].probability = max_prob;
        calls[i].state_index = max_state;
        calls[i].segment_mean = data[i];
    }

    return calls;
}

double FastCall::compute_log_likelihood(
    const std::vector<double>& data,
    const std::vector<std::vector<double>>& posteriors
) {
    // Compute weighted log-likelihood
    double ll = 0.0;
    for (size_t i = 0; i < data.size(); ++i) {
        double sum = 0.0;
        for (int j = 0; j < N_STATES; ++j) {
            double pdf = normal_pdf(data[i], means_[j], sds_[j]);
            sum += priors_[j] * pdf;
        }
        if (sum > 0) {
            ll += std::log(sum);
        }
    }
    return ll;
}

} // namespace fastcall
} // namespace excavator
