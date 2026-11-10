// SPDX-License-Identifier: MIT
#include "cylload/montecarlo.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include "cylload/analysis.hpp"
#include "cylload/linalg.hpp"

namespace cylload {

// ------------------------------------------------------------ mt19937_64
Uniform::Uniform(std::uint64_t seed) : idx_(312) {
    mt_[0] = seed;
    for (int i = 1; i < 312; ++i)
        mt_[i] = 6364136223846793005ULL * (mt_[i - 1] ^ (mt_[i - 1] >> 62)) + static_cast<std::uint64_t>(i);
}

std::uint64_t Uniform::next() {
    constexpr std::uint64_t UM = 0xFFFFFFFF80000000ULL, LM = 0x7FFFFFFFULL, A = 0xB5026F5AA96619E9ULL;
    if (idx_ >= 312) {
        for (int i = 0; i < 312; ++i) {
            std::uint64_t x = (mt_[i] & UM) | (mt_[(i + 1) % 312] & LM);
            std::uint64_t xa = x >> 1;
            if (x & 1ULL) xa ^= A;
            mt_[i] = mt_[(i + 156) % 312] ^ xa;
        }
        idx_ = 0;
    }
    std::uint64_t x = mt_[idx_++];
    x ^= (x >> 29) & 0x5555555555555555ULL;
    x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
    x ^= (x << 37) & 0xFFF7EEE000000000ULL;
    x ^= (x >> 43);
    return x;
}

double Uniform::operator()() { return static_cast<double>(next() >> 11) * (1.0 / 9007199254740992.0); }

// ------------------------------------------------------------ spec parser
namespace {
std::vector<std::string> words(const std::string& line) {
    std::string s = line.substr(0, line.find('#'));
    std::istringstream is(s);
    std::vector<std::string> t;
    std::string w;
    while (is >> w) t.push_back(w);
    return t;
}
std::vector<std::string> split_commas(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == ',') { out.push_back(cur); cur.clear(); } else cur += c;
    }
    out.push_back(cur);
    return out;
}
double d(const std::string& s) {
    std::size_t p = 0;
    double v = std::stod(s, &p);
    if (p != s.size()) throw std::runtime_error("bad number '" + s + "'");
    return v;
}
bool selected(const MassItem& mi, const std::string& body_name, const std::vector<std::string>& sel) {
    for (const auto& s : sel) {
        if (s == "all") return true;
        auto c = s.find(':');
        if (c == std::string::npos) continue;
        const std::string k = s.substr(0, c), v = s.substr(c + 1);
        if (k == "group" && mi.group == v) return true;
        if (k == "item" && mi.name == v) return true;
        if (k == "body" && body_name == v) return true;
    }
    return false;
}
bool cylinder_selected(const Actuator& a, const std::string& body_a, const std::vector<std::string>& sel) {
    for (const auto& s : sel) {
        if (s == "all") return true;
        if (s == "group:" + a.name || s == "item:" + a.name) return true;
        (void)body_a;
    }
    return false;
}
}  // namespace

McSpec McSpec::from_string(const std::string& text) {
    McSpec sp;
    std::istringstream is(text);
    std::string line;
    int ln = 0;
    while (std::getline(is, line)) {
        ++ln;
        auto t = words(line);
        if (t.empty()) continue;
        try {
            const auto& k = t[0];
            if (k == "seed") sp.seed = std::stoull(t.at(1));
            else if (k == "samples") sp.samples = std::stoi(t.at(1));
            else if (k == "bootstrap") sp.bootstrap = std::stoi(t.at(1));
            else if (k == "scale" || k == "mass" || k == "length") {
                // scale LABEL SELECTORS LO HI | mass LABEL SELECTORS LO HI | length LABEL ACTUATOR LO HI
                if (t.size() != 5) throw std::runtime_error("usage: " + k + " LABEL SELECTORS LO HI");
                McVariable v;
                v.kind = k == "scale" ? McVariable::Kind::Scale : k == "mass" ? McVariable::Kind::Mass : McVariable::Kind::Length;
                v.label = t[1];
                v.selectors = split_commas(t[2]);
                v.lo = d(t[3]);
                v.hi = d(t[4]);
                sp.vars.push_back(v);
            } else if (k == "shift") {  // shift LABEL SELECTORS DXLO DXHI DYLO DYHI -> LABEL_dx, LABEL_dy
                if (t.size() != 7) throw std::runtime_error("usage: shift LABEL SELECTORS DXLO DXHI DYLO DYHI");
                McVariable vx{McVariable::Kind::ShiftX, t[1] + "_dx", split_commas(t[2]), d(t[3]), d(t[4])};
                McVariable vy{McVariable::Kind::ShiftY, t[1] + "_dy", split_commas(t[2]), d(t[5]), d(t[6])};
                sp.vars.push_back(vx);
                sp.vars.push_back(vy);
            } else if (k == "output") {  // output LABEL ACTUATOR [NAME=L ...]
                if (t.size() < 3) throw std::runtime_error("usage: output LABEL ACTUATOR [NAME=L ...]");
                McOutput o;
                o.label = t[1];
                o.actuator = t[2];
                o.lengths = parse_assignments(std::vector<std::string>(t.begin() + 3, t.end()));
                sp.outputs.push_back(o);
            } else throw std::runtime_error("unknown keyword '" + k + "'");
        } catch (const std::exception& e) {
            throw std::runtime_error("Monte Carlo spec line " + std::to_string(ln) + ": " + e.what());
        }
    }
    if (sp.outputs.empty()) throw std::runtime_error("Monte Carlo spec: no output defined");
    if (sp.samples < 2) throw std::runtime_error("Monte Carlo spec: samples must be >= 2");
    return sp;
}

McSpec McSpec::from_file(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open Monte Carlo spec '" + path + "'");
    std::stringstream ss;
    ss << f.rdbuf();
    return from_string(ss.str());
}

double percentile(std::vector<double> v, double q) {
    if (v.empty()) return NAN;
    std::sort(v.begin(), v.end());
    const double pos = q / 100.0 * static_cast<double>(v.size() - 1);
    const std::size_t i = static_cast<std::size_t>(std::floor(pos));
    if (i + 1 >= v.size()) return v.back();
    return v[i] + (pos - static_cast<double>(i)) * (v[i + 1] - v[i]);
}

namespace {
std::vector<double> ranks(const std::vector<double>& x) {
    std::vector<std::size_t> idx(x.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return x[a] < x[b]; });
    std::vector<double> r(x.size());
    for (std::size_t i = 0; i < idx.size();) {
        std::size_t j = i;
        while (j + 1 < idx.size() && x[idx[j + 1]] == x[idx[i]]) ++j;
        const double avg = 0.5 * static_cast<double>(i + j) + 1.0;
        for (std::size_t k = i; k <= j; ++k) r[idx[k]] = avg;
        i = j + 1;
    }
    return r;
}
}  // namespace

double spearman(const std::vector<double>& x, const std::vector<double>& y) {
    auto rx = ranks(x), ry = ranks(y);
    const double n = static_cast<double>(x.size());
    const double mx = std::accumulate(rx.begin(), rx.end(), 0.0) / n, my = std::accumulate(ry.begin(), ry.end(), 0.0) / n;
    double sxy = 0, sxx = 0, syy = 0;
    for (std::size_t i = 0; i < rx.size(); ++i) {
        sxy += (rx[i] - mx) * (ry[i] - my);
        sxx += (rx[i] - mx) * (rx[i] - mx);
        syy += (ry[i] - my) * (ry[i] - my);
    }
    return (sxx > 0 && syy > 0) ? sxy / std::sqrt(sxx * syy) : 0.0;
}

double first_order_index(const std::vector<double>& x, const std::vector<double>& y, int bins) {
    const std::size_t n = x.size();
    if (n < static_cast<std::size_t>(2 * bins)) return NAN;
    std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return x[a] < x[b]; });
    const double mean = std::accumulate(y.begin(), y.end(), 0.0) / static_cast<double>(n);
    double var = 0;
    for (double v : y) var += (v - mean) * (v - mean);
    var /= static_cast<double>(n);
    if (var <= 0) return 0.0;
    double between = 0;
    for (int b = 0; b < bins; ++b) {
        const std::size_t lo = n * static_cast<std::size_t>(b) / static_cast<std::size_t>(bins);
        const std::size_t hi = n * static_cast<std::size_t>(b + 1) / static_cast<std::size_t>(bins);
        if (hi <= lo) continue;
        double s = 0;
        for (std::size_t k = lo; k < hi; ++k) s += y[idx[k]];
        const double mb = s / static_cast<double>(hi - lo);
        between += static_cast<double>(hi - lo) * (mb - mean) * (mb - mean);
    }
    return between / static_cast<double>(n) / var;
}

std::vector<double> standardized_regression(const std::vector<std::vector<double>>& X, const std::vector<double>& y, double& r2) {
    const std::size_t p = X.size(), n = y.size();
    auto stdz = [n](const std::vector<double>& v, double& sd) {
        const double m = std::accumulate(v.begin(), v.end(), 0.0) / static_cast<double>(n);
        double ss = 0;
        for (double a : v) ss += (a - m) * (a - m);
        sd = std::sqrt(ss / static_cast<double>(n));
        std::vector<double> z(n, 0.0);
        if (sd > 0)
            for (std::size_t i = 0; i < n; ++i) z[i] = (v[i] - m) / sd;
        return z;
    };
    double sy;
    auto zy = stdz(y, sy);
    std::vector<std::vector<double>> Z;
    std::vector<std::size_t> used;
    for (std::size_t j = 0; j < p; ++j) {
        double sx;
        auto z = stdz(X[j], sx);
        if (sx > 0) { Z.push_back(z); used.push_back(j); }
    }
    std::vector<double> beta(p, 0.0);
    r2 = 0;
    if (Z.empty() || sy <= 0) return beta;
    const std::size_t q = Z.size();
    Matrix A(q, q);
    std::vector<double> b(q, 0.0);
    for (std::size_t a = 0; a < q; ++a) {
        for (std::size_t c = 0; c < q; ++c) {
            double s = 0;
            for (std::size_t i = 0; i < n; ++i) s += Z[a][i] * Z[c][i];
            A(a, c) = s;
        }
        for (std::size_t i = 0; i < n; ++i) b[a] += Z[a][i] * zy[i];
    }
    auto x = solve(A, b);
    double ssr = 0;
    for (std::size_t i = 0; i < n; ++i) {
        double f = 0;
        for (std::size_t a = 0; a < q; ++a) f += x[a] * Z[a][i];
        ssr += (zy[i] - f) * (zy[i] - f);
    }
    r2 = 1.0 - ssr / static_cast<double>(n);
    for (std::size_t a = 0; a < q; ++a) beta[used[a]] = x[a];
    return beta;
}

McResult run_monte_carlo(const Model& base, const McSpec& spec) {
    McResult res;
    for (const auto& v : spec.vars) res.input_labels.push_back(v.label);
    for (const auto& o : spec.outputs) res.output_labels.push_back(o.label);
    std::vector<int> out_act;
    for (const auto& o : spec.outputs) {
        out_act.push_back(base.actuator_index(o.actuator));
        for (const auto& kv : o.lengths) base.actuator_index(kv.first);
    }
    for (const auto& v : spec.vars)
        if (v.kind == McVariable::Kind::Length) base.actuator_index(v.selectors.at(0));
    Uniform rnd(spec.seed);
    for (int n = 0; n < spec.samples; ++n) {
        Model m = base;
        std::vector<double> x;
        std::map<std::string, double> len_override;
        for (const auto& v : spec.vars) {
            const double u = rnd(v.lo, v.hi);
            x.push_back(u);
            for (auto& mi : m.masses) {
                if (!selected(mi, m.bodies[mi.body], v.selectors)) continue;
                switch (v.kind) {
                    case McVariable::Kind::Scale: mi.m *= u; break;
                    case McVariable::Kind::Mass: mi.m = u; break;
                    case McVariable::Kind::ShiftX: mi.c[0] += u; break;
                    case McVariable::Kind::ShiftY: mi.c[1] += u; break;
                    default: break;
                }
            }
            if (v.kind == McVariable::Kind::Scale)
                for (auto& a : m.actuators)
                    if (a.has_cylinder && cylinder_selected(a, m.bodies[a.body_a], v.selectors)) a.cyl.mass *= u;
            if (v.kind == McVariable::Kind::Length) len_override[v.selectors[0]] = u;
        }
        Kinematics kin(m);
        std::vector<double> f, p;
        for (std::size_t i = 0; i < spec.outputs.size(); ++i) {
            const auto& o = spec.outputs[i];
            Pose s = kin.reference();
            auto L = s.lengths;
            for (const auto& kv : len_override) L[m.actuator_index(kv.first)] = kv.second;
            for (const auto& kv : o.lengths) L[m.actuator_index(kv.first)] = kv.second;
            s = kin.solve(s, L);
            const double F = forces_virtual_work(kin, s)[out_act[i]];
            const auto& a = m.actuators[out_act[i]];
            f.push_back(F);
            p.push_back(a.has_cylinder ? pressure_bar(a.cyl, F) : 0.0);
        }
        res.inputs.push_back(x);
        res.force.push_back(f);
        res.pressure.push_back(p);
    }
    auto stats = [](const std::vector<double>& v) {
        McStats s;
        s.mean = std::accumulate(v.begin(), v.end(), 0.0) / static_cast<double>(v.size());
        s.min = *std::min_element(v.begin(), v.end());
        s.max = *std::max_element(v.begin(), v.end());
        s.p5 = percentile(v, 5);
        s.p50 = percentile(v, 50);
        s.p95 = percentile(v, 95);
        return s;
    };
    for (std::size_t i = 0; i < spec.outputs.size(); ++i) {
        std::vector<double> f, p;
        for (int n = 0; n < spec.samples; ++n) {
            f.push_back(res.force[n][i]);
            p.push_back(res.pressure[n][i]);
        }
        res.force_stats.push_back(stats(f));
        res.pressure_stats.push_back(stats(p));
        // rank correlation on the pressure (what the hydraulic design depends on);
        // on the force for actuators without cylinder data
        const bool cyl = base.actuators[out_act[i]].has_cylinder;
        std::vector<double> rho;
        for (std::size_t j = 0; j < spec.vars.size(); ++j) {
            std::vector<double> xj;
            for (int n = 0; n < spec.samples; ++n) xj.push_back(res.inputs[n][j]);
            rho.push_back(spearman(xj, cyl ? p : f));
        }
        res.spearman.push_back(rho);
        std::vector<std::vector<double>> X;
        std::vector<double> s1;
        for (std::size_t j = 0; j < spec.vars.size(); ++j) {
            std::vector<double> xj;
            for (int n = 0; n < spec.samples; ++n) xj.push_back(res.inputs[static_cast<std::size_t>(n)][j]);
            s1.push_back(first_order_index(xj, cyl ? p : f));
            X.push_back(xj);
        }
        res.first_order.push_back(s1);
        // percentile bootstrap: resample the N runs with replacement and recompute the estimators
        {
            const std::size_t N = static_cast<std::size_t>(spec.samples), nv = spec.vars.size();
            const int B = spec.bootstrap;
            std::vector<double> lo(nv, NAN), hi(nv, NAN);
            if (B > 0) {
                Uniform br(spec.seed ^ 0x9E3779B97F4A7C15ULL ^ static_cast<std::uint64_t>(i));
                const auto& yv = cyl ? p : f;
                std::vector<std::vector<double>> sb(nv);
                std::vector<double> b5, b95;
                std::vector<double> yb(N), xb(N);
                std::vector<std::size_t> pick(N);
                for (int b = 0; b < B; ++b) {
                    for (auto& k : pick) k = std::min(N - 1, static_cast<std::size_t>(br() * static_cast<double>(N)));
                    for (std::size_t k = 0; k < N; ++k) yb[k] = yv[pick[k]];
                    b5.push_back(percentile(yb, 5));
                    b95.push_back(percentile(yb, 95));
                    for (std::size_t j = 0; j < nv; ++j) {
                        for (std::size_t k = 0; k < N; ++k) xb[k] = X[j][pick[k]];
                        sb[j].push_back(first_order_index(xb, yb));
                    }
                }
                auto& st = cyl ? res.pressure_stats.back() : res.force_stats.back();
                st.p5_lo = percentile(b5, 2.5);
                st.p5_hi = percentile(b5, 97.5);
                st.p95_lo = percentile(b95, 2.5);
                st.p95_hi = percentile(b95, 97.5);
                // S_i: basic (reverse-percentile) interval. Resampling with replacement repeats
                // points, which raises the binning bias of the estimator; reflecting the bootstrap
                // distribution about the estimate removes that shift to first order.
                for (std::size_t j = 0; j < nv; ++j) {
                    lo[j] = std::max(0.0, 2.0 * s1[j] - percentile(sb[j], 97.5));
                    hi[j] = std::max(0.0, 2.0 * s1[j] - percentile(sb[j], 2.5));
                }
            }
            res.first_order_lo.push_back(lo);
            res.first_order_hi.push_back(hi);
        }
        double r2 = 0;
        res.src.push_back(standardized_regression(X, cyl ? p : f, r2));
        res.r2.push_back(r2);
    }
    return res;
}

}  // namespace cylload
