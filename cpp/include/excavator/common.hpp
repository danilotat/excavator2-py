#pragma once

#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>
#include <numeric>






namespace excavator {

// Numerical constants
constexpr double PI = 3.14159265358979323846;
constexpr double LOG_ZERO = -std::numeric_limits<double>::infinity();
constexpr double EPSILON = 1e-300;

/**
 * Log-sum-exp: compute log(exp(x) + exp(y)) in a numerically stable way.
 * This is the core utility for log-space arithmetic in HMMs.
 *
 * Ported from Fortran ELNSUM function in FastJointSLMLibraryI.f
 */
inline double elnsum(double x, double y) {
    if (x > y) {
        return x + std::log(1.0 + std::exp(y - x));
    } else {
        return y + std::log(1.0 + std::exp(x - y));
    }
}

/**
 * Log-sum-exp for a vector: compute log(sum(exp(v[i]))) in a numerically stable way.
 */
inline double logsumexp(const std::vector<double>& v) {
    if (v.empty()) {
        return LOG_ZERO;
    }

    double max_val = *std::max_element(v.begin(), v.end());
    if (std::isinf(max_val)) {
        return max_val;
    }

    double sum = 0.0;
    for (double x : v) {
        sum += std::exp(x - max_val);
    }
    return max_val + std::log(sum);
}

/**
 * Compute log of normal PDF.
 * log(N(x | mu, sigma)) = -0.5*log(2*pi*sigma^2) - (x-mu)^2/(2*sigma^2)
 */
inline double log_normal_pdf(double x, double mu, double sigma) {
    double z = (x - mu) / sigma;
    return -0.5 * std::log(2.0 * PI) - std::log(sigma) - 0.5 * z * z;
}

/**
 * Compute variance with optional trimming (quantile-based).
 * Used for robust variance estimation.
 */
double trimmed_variance(const std::vector<double>& data, double lower_q = 0.01, double upper_q = 0.99);

/**
 * Compute median of a vector.
 */
double median(std::vector<double> data);

/**
 * Compute quantile of a vector.
 */
double quantile(std::vector<double> data, double p);

} // namespace excavator
