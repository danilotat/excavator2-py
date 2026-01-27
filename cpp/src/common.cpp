#include "excavator/common.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace excavator {

double trimmed_variance(const std::vector<double>& data, double lower_q, double upper_q) {
    if (data.empty()) {
        return 0.0;
    }

    if (data.size() == 1) {
        return 0.0;
    }

    // Get quantile bounds
    double lower_bound = quantile(data, lower_q);
    double upper_bound = quantile(data, upper_q);

    // Filter data within bounds
    std::vector<double> filtered;
    filtered.reserve(data.size());
    for (double x : data) {
        if (x >= lower_bound && x <= upper_bound) {
            filtered.push_back(x);
        }
    }

    if (filtered.empty()) {
        // Fall back to full data if filtering removes everything
        filtered = data;
    }

    // Compute variance
    double n = static_cast<double>(filtered.size());
    double mean = std::accumulate(filtered.begin(), filtered.end(), 0.0) / n;
    double sum_sq = 0.0;
    for (double x : filtered) {
        double diff = x - mean;
        sum_sq += diff * diff;
    }

    return sum_sq / (n - 1.0);  // Sample variance
}

double median(std::vector<double> data) {
    if (data.empty()) {
        throw std::invalid_argument("Cannot compute median of empty vector");
    }

    size_t n = data.size();
    std::sort(data.begin(), data.end());

    if (n % 2 == 0) {
        return (data[n / 2 - 1] + data[n / 2]) / 2.0;
    } else {
        return data[n / 2];
    }
}

double quantile(std::vector<double> data, double p) {
    if (data.empty()) {
        throw std::invalid_argument("Cannot compute quantile of empty vector");
    }

    if (p < 0.0 || p > 1.0) {
        throw std::invalid_argument("Quantile p must be in [0, 1]");
    }

    std::sort(data.begin(), data.end());

    size_t n = data.size();
    double index = p * (n - 1);
    size_t lower_idx = static_cast<size_t>(std::floor(index));
    size_t upper_idx = static_cast<size_t>(std::ceil(index));

    if (lower_idx == upper_idx) {
        return data[lower_idx];
    }

    double frac = index - lower_idx;
    return data[lower_idx] * (1.0 - frac) + data[upper_idx] * frac;
}

} // namespace excavator
