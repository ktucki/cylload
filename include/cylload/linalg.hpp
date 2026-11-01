// cylload - hydraulic cylinder loads from a planar rigid-body model of a machine
// SPDX-License-Identifier: MIT
#pragma once
#include <cstddef>
#include <vector>

namespace cylload {

/// Dense row-major square/rectangular matrix, sufficient for the small systems
/// (3 x number of bodies) that arise in planar mechanisms.
class Matrix {
public:
    Matrix() = default;
    Matrix(std::size_t rows, std::size_t cols, double v = 0.0) : r_(rows), c_(cols), a_(rows * cols, v) {}
    double& operator()(std::size_t i, std::size_t j) { return a_[i * c_ + j]; }
    double operator()(std::size_t i, std::size_t j) const { return a_[i * c_ + j]; }
    std::size_t rows() const { return r_; }
    std::size_t cols() const { return c_; }

private:
    std::size_t r_ = 0, c_ = 0;
    std::vector<double> a_;
};

/// LU factorisation with partial pivoting of a square matrix.
class LU {
public:
    explicit LU(const Matrix& A);
    bool singular() const { return singular_; }
    /// smallest |pivot| / largest |pivot|: a cheap conditioning indicator
    double pivot_ratio() const { return pivot_ratio_; }
    std::vector<double> solve(const std::vector<double>& b) const;

private:
    Matrix lu_;
    std::vector<std::size_t> piv_;
    bool singular_ = false;
    double pivot_ratio_ = 0.0;
};

std::vector<double> solve(const Matrix& A, const std::vector<double>& b);
Matrix transpose(const Matrix& A);
double norm_inf(const std::vector<double>& v);

}  // namespace cylload
