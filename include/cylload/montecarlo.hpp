// SPDX-License-Identifier: MIT
#pragma once
#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "cylload/model.hpp"

namespace cylload {

/// One uniformly distributed uncertain input.
struct McVariable {
    enum class Kind { Scale, Mass, ShiftX, ShiftY, Length } kind;
    std::string label;
    std::vector<std::string> selectors;  ///< "all", "group:G", "item:NAME", "body:B" (or actuator name for Length)
    double lo = 0, hi = 0;
};

/// Output: force and pressure of an actuator in a pose given by actuator lengths.
struct McOutput {
    std::string label;
    std::string actuator;
    std::map<std::string, double> lengths;  ///< unspecified actuators stay at reference length
};

struct McSpec {
    std::uint64_t seed = 2026;
    int samples = 1000;
    int bootstrap = 1000;  ///< resamples for the confidence intervals (0 = none)
    std::vector<McVariable> vars;
    std::vector<McOutput> outputs;
    static McSpec from_file(const std::string& path);
    static McSpec from_string(const std::string& text);
};

struct McStats {
    double mean = 0, min = 0, p5 = 0, p50 = 0, p95 = 0, max = 0;
    double p5_lo = 0, p5_hi = 0, p95_lo = 0, p95_hi = 0;  ///< 95 % percentile bootstrap intervals
};

struct McResult {
    std::vector<std::string> input_labels, output_labels;
    std::vector<std::vector<double>> inputs;   ///< [sample][variable]
    std::vector<std::vector<double>> force;    ///< [sample][output] N
    std::vector<std::vector<double>> pressure; ///< [sample][output] bar (0 without cylinder data)
    std::vector<McStats> force_stats, pressure_stats;
    std::vector<std::vector<double>> spearman; ///< [output][variable], on the pressure (force without cylinder)
    /// first-order variance-based indices S_i = Var(E[Y|X_i]) / Var(Y), estimated from the
    /// sample with 20 equal-probability bins of each input
    std::vector<std::vector<double>> first_order;
    std::vector<std::vector<double>> first_order_lo, first_order_hi;  ///< 95 % basic bootstrap intervals
    std::vector<std::vector<double>> src;      ///< standardised regression coefficients
    std::vector<double> r2;                    ///< coefficient of determination of that regression
};

/// S_i of y with respect to x from a sample (equal-probability bins).
double first_order_index(const std::vector<double>& x, const std::vector<double>& y, int bins = 20);
/// Standardised regression coefficients of y on the columns X[var][sample]; r2 receives R^2.
std::vector<double> standardized_regression(const std::vector<std::vector<double>>& X, const std::vector<double>& y, double& r2);

/// Small, portable random source: mt19937_64 bits mapped to [0,1) with 53-bit resolution,
/// so a seed gives the same samples with every standard library.
class Uniform {
public:
    explicit Uniform(std::uint64_t seed);
    double operator()();
    double operator()(double lo, double hi) { return lo + (hi - lo) * (*this)(); }

private:
    std::uint64_t mt_[312];
    int idx_;
    std::uint64_t next();
};

McResult run_monte_carlo(const Model& base, const McSpec& spec);
double percentile(std::vector<double> v, double q);  ///< linear interpolation, q in [0,100]
double spearman(const std::vector<double>& x, const std::vector<double>& y);

}  // namespace cylload
