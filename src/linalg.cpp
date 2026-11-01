// SPDX-License-Identifier: MIT
#include "cylload/linalg.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cylload {

LU::LU(const Matrix& A) : lu_(A), piv_(A.rows()) {
    if (A.rows() != A.cols()) throw std::invalid_argument("LU: matrix is not square");
    const std::size_t n = A.rows();
    for (std::size_t i = 0; i < n; ++i) piv_[i] = i;
    double pmin = 1e300, pmax = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t p = k;
        double best = std::fabs(lu_(k, k));
        for (std::size_t i = k + 1; i < n; ++i)
            if (std::fabs(lu_(i, k)) > best) { best = std::fabs(lu_(i, k)); p = i; }
        pmin = std::min(pmin, best);
        pmax = std::max(pmax, best);
        if (best < 1e-300) { singular_ = true; pivot_ratio_ = 0.0; return; }
        if (p != k) {
            for (std::size_t j = 0; j < n; ++j) std::swap(lu_(k, j), lu_(p, j));
            std::swap(piv_[k], piv_[p]);
        }
        for (std::size_t i = k + 1; i < n; ++i) {
            const double f = lu_(i, k) / lu_(k, k);
            lu_(i, k) = f;
            for (std::size_t j = k + 1; j < n; ++j) lu_(i, j) -= f * lu_(k, j);
        }
    }
    pivot_ratio_ = pmax > 0 ? pmin / pmax : 0.0;
    if (pivot_ratio_ < 1e-13) singular_ = true;
}

std::vector<double> LU::solve(const std::vector<double>& b) const {
    if (singular_) throw std::runtime_error("LU: singular matrix (mechanism in a singular or locked position?)");
    const std::size_t n = lu_.rows();
    std::vector<double> x(n);
    for (std::size_t i = 0; i < n; ++i) x[i] = b[piv_[i]];
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < i; ++j) x[i] -= lu_(i, j) * x[j];
    for (std::size_t ii = n; ii-- > 0;) {
        for (std::size_t j = ii + 1; j < n; ++j) x[ii] -= lu_(ii, j) * x[j];
        x[ii] /= lu_(ii, ii);
    }
    return x;
}

std::vector<double> solve(const Matrix& A, const std::vector<double>& b) { return LU(A).solve(b); }

Matrix transpose(const Matrix& A) {
    Matrix T(A.cols(), A.rows());
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) T(j, i) = A(i, j);
    return T;
}

double norm_inf(const std::vector<double>& v) {
    double m = 0.0;
    for (double x : v) m = std::max(m, std::fabs(x));
    return m;
}

}  // namespace cylload
