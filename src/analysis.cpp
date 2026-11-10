// SPDX-License-Identifier: MIT
#include "cylload/analysis.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cylload {

namespace {
constexpr double kPi = 3.14159265358979323846;
double cross(const Vec2& r, const Vec2& f) { return r[0] * f[1] - r[1] * f[0]; }
}  // namespace

double potential(const Kinematics& kin, const Pose& s) {
    const auto& m = kin.model();
    double V = 0.0;
    for (const auto& mi : m.all_masses()) {
        Vec2 r = kin.point(s, mi.body, mi.c);
        V -= mi.m * (m.gravity[0] * r[0] + m.gravity[1] * r[1]) / 1000.0;  // J (r in mm)
    }
    for (const auto& ef : m.forces) {  // a constant force has the potential -f.r
        Vec2 r = kin.point(s, ef.body, ef.p);
        V -= (ef.f[0] * r[0] + ef.f[1] * r[1]) / 1000.0;
    }
    return V;
}

std::vector<double> forces_virtual_work(const Kinematics& kin, const Pose& s) {
    const auto& m = kin.model();
    Matrix dq = kin.sensitivity(s);
    std::vector<double> F(m.actuators.size(), 0.0);
    const auto masses = m.all_masses();
    for (std::size_t k = 0; k < F.size(); ++k)
        for (const auto& mi : masses) {
            Vec2 v = kin.point_rate(s, dq, static_cast<int>(k), mi.body, mi.c);
            F[k] -= mi.m * (m.gravity[0] * v[0] + m.gravity[1] * v[1]);  // N
        }
    for (std::size_t k = 0; k < F.size(); ++k)
        for (const auto& ef : m.forces) {
            Vec2 v = kin.point_rate(s, dq, static_cast<int>(k), ef.body, ef.p);
            F[k] -= ef.f[0] * v[0] + ef.f[1] * v[1];
        }
    return F;
}

std::vector<double> forces_finite_difference(const Kinematics& kin, const Pose& s, double h) {
    std::vector<double> F(s.lengths.size());
    for (std::size_t k = 0; k < F.size(); ++k) {
        auto lp = s.lengths, lm = s.lengths;
        lp[k] += h;
        lm[k] -= h;
        Pose sp = kin.solve(s, lp), sm = kin.solve(s, lm);
        F[k] = (potential(kin, sp) - potential(kin, sm)) / (2.0 * h) * 1000.0;
    }
    return F;
}

StaticResult statics(const Kinematics& kin, const Pose& s) {
    const auto& m = kin.model();
    const int nb = m.n_moving();
    const int nr = static_cast<int>(m.revolutes.size()), np = static_cast<int>(m.prismatics.size());
    const int na = static_cast<int>(m.actuators.size());
    const int nu = 2 * nr + 2 * np + na;
    if (nu != 3 * nb) throw std::runtime_error("statics: model is not statically determinate");
    Matrix A(3 * nb, nu);
    std::vector<double> rhs(3 * nb, 0.0);
    // force f (per unit unknown) applied to body b at point r: rows 3(b-1) .. +2
    auto apply = [&](int body, int col, const Vec2& f, const Vec2& r) {
        if (body == 0) return;
        const int row = 3 * (body - 1);
        A(row, col) += f[0];
        A(row + 1, col) += f[1];
        A(row + 2, col) += cross(r, f) / 1000.0;  // N m
    };
    int col = 0;
    for (const auto& j : m.revolutes) {
        Vec2 r = kin.point(s, j.body_b, j.p);
        for (int c = 0; c < 2; ++c) {
            Vec2 f{c == 0 ? 1.0 : 0.0, c == 1 ? 1.0 : 0.0};
            apply(j.body_b, col, f, r);
            apply(j.body_a, col, {-f[0], -f[1]}, r);
            ++col;
        }
    }
    for (const auto& j : m.prismatics) {
        Vec2 r = kin.point(s, j.body_b, j.p);
        const double ta = kin.angle(s, j.body_a);
        Vec2 n{-(std::cos(ta) * j.u[1] + std::sin(ta) * j.u[0]), std::cos(ta) * j.u[0] - std::sin(ta) * j.u[1]};
        apply(j.body_b, col, n, r);
        apply(j.body_a, col, {-n[0], -n[1]}, r);
        ++col;
        if (j.body_b) A(3 * (j.body_b - 1) + 2, col) += 1.0;  // couple [N m]
        if (j.body_a) A(3 * (j.body_a - 1) + 2, col) -= 1.0;
        ++col;
    }
    for (int k = 0; k < na; ++k) {
        const auto& a = m.actuators[k];
        Vec2 u = kin.actuator_dir(s, k);
        apply(a.body_b, col, u, kin.point(s, a.body_b, a.pb));
        apply(a.body_a, col, {-u[0], -u[1]}, kin.point(s, a.body_a, a.pa));
        ++col;
    }
    for (const auto& mi : m.all_masses()) {
        if (mi.body == 0) continue;
        Vec2 r = kin.point(s, mi.body, mi.c);
        Vec2 w{mi.m * m.gravity[0], mi.m * m.gravity[1]};
        const int row = 3 * (mi.body - 1);
        rhs[row] -= w[0];
        rhs[row + 1] -= w[1];
        rhs[row + 2] -= cross(r, w) / 1000.0;
    }
    for (const auto& ef : m.forces) {
        if (ef.body == 0) continue;
        Vec2 r = kin.point(s, ef.body, ef.p);
        const int row = 3 * (ef.body - 1);
        rhs[row] -= ef.f[0];
        rhs[row + 1] -= ef.f[1];
        rhs[row + 2] -= cross(r, ef.f) / 1000.0;
    }
    auto x = solve(A, rhs);
    StaticResult out;
    col = 0;
    for (const auto& j : m.revolutes) {
        out.reactions.push_back({j.name, x[col], x[col + 1], 0.0});
        col += 2;
    }
    for (const auto& j : m.prismatics) {
        const double ta = kin.angle(s, j.body_a);
        Vec2 n{-(std::cos(ta) * j.u[1] + std::sin(ta) * j.u[0]), std::cos(ta) * j.u[0] - std::sin(ta) * j.u[1]};
        out.reactions.push_back({j.name, x[col] * n[0], x[col] * n[1], x[col + 1]});
        col += 2;
    }
    for (int k = 0; k < na; ++k) out.force_eq.push_back(x[col++]);
    out.force_vw = forces_virtual_work(kin, s);
    out.potential = potential(kin, s);
    return out;
}

double pressure_bar(const Cylinder& c, double F) {
    // loaded chamber; a back-pressure in the other chamber adds p_back times the area ratio
    const double Ap = c.area_piston(), Ar = c.area_annulus();
    return F >= 0 ? F / Ap * 10.0 + c.p_back * Ar / Ap : -F / Ar * 10.0 + c.p_back * Ap / Ar;
}

double reduced_mass(const Kinematics& kin, const Pose& s, int k) {
    const auto& m = kin.model();
    Matrix dq = kin.sensitivity(s);
    double mr = 0.0;
    for (const auto& mi : m.all_masses()) {
        if (mi.body == 0) continue;
        if (!m.hydraulics.pin_masses_in_reduced_mass && mi.name.size() > 6 &&
            (mi.name.compare(mi.name.size() - 6, 6, "/pin_a") == 0 || mi.name.compare(mi.name.size() - 6, 6, "/pin_b") == 0))
            continue;
        Vec2 v = kin.point_rate(s, dq, k, mi.body, mi.c);
        const double w = dq(3 * (mi.body - 1) + 2, k) * 1000.0;  // rad per m of stroke
        mr += mi.m * (v[0] * v[0] + v[1] * v[1]) + mi.j * w * w;
    }
    return mr;
}

double hydraulic_stiffness(const Cylinder& c, const Hydraulics& h, double L) {
    const double Ap = c.area_piston() * 1e-6, Ar = c.area_annulus() * 1e-6;
    const double x = (L - c.lmin) / 1000.0, s = c.stroke() / 1000.0, v0 = h.dead_volume * 1e-3, v0r = h.dead_volume_rod * 1e-3;
    return h.beta * (Ap * Ap / (v0 + Ap * x) + Ar * Ar / (v0r + Ar * (s - x)));
}

double ramp_factor(double tau) {
    if (tau <= 0) return 1.0;
    return tau < kPi ? 2.0 * std::fabs(std::sin(tau / 2.0)) / tau : 2.0 / tau;
}

double euler_load(const Cylinder& c, double L, double E) {
    const double I = kPi * std::pow(c.rod, 4) / 64.0;
    return kPi * kPi * E * I / (L * L);
}

std::vector<SweepPoint> sweep(const Kinematics& kin, const SweepOptions& o) {
    const auto& m = kin.model();
    if (o.n < 2) throw std::invalid_argument("sweep: n must be >= 2");
    if (o.fixed.size() != m.actuators.size()) throw std::invalid_argument("sweep: wrong number of actuator lengths");
    const int rk = o.report < 0 ? o.actuator : o.report;
    const auto& act = m.actuators[rk];
    std::vector<SweepPoint> out;
    auto L = o.fixed;
    L[o.actuator] = o.from;
    Pose s = kin.solve(kin.reference(), L);
    for (int i = 0; i < o.n; ++i) {
        L[o.actuator] = o.from + (o.to - o.from) * i / (o.n - 1);
        s = kin.solve(s, L);
        SweepPoint p;
        p.L = L[o.actuator];
        for (int b = 1; b <= m.n_moving(); ++b) p.angle_deg.push_back(kin.angle(s, b) * 180.0 / kPi);
        auto st = statics(kin, s);
        p.F = st.force_vw[rk];
        p.F_eq = st.force_eq[rk];
        p.F_fd = forces_finite_difference(kin, s, o.fd_step)[rk];
        p.V = st.potential;
        p.m_red = reduced_mass(kin, s, rk);
        if (act.has_cylinder) {
            p.p = pressure_bar(act.cyl, p.F);
            p.k_hyd = hydraulic_stiffness(act.cyl, m.hydraulics, L[rk]);  // at the reported cylinder's length
            p.f_n = p.m_red > 0 ? std::sqrt(p.k_hyd / p.m_red) / (2.0 * kPi) : 0.0;
        }
        out.push_back(p);
    }
    return out;
}

SweepSummary summarize(const Kinematics& kin, const SweepOptions& o, const std::vector<SweepPoint>& pts) {
    const int rk = o.report < 0 ? o.actuator : o.report;
    const auto& act = kin.model().actuators[rk];
    SweepSummary s;
    s.F_min = s.F_max = pts[0].F;
    s.L_at_F_min = s.L_at_F_max = pts[0].L;
    double fabsmax = 0, dfd = 0, deq = 0;
    s.f_min = s.f_max = pts[0].f_n;
    s.m_red_min = s.m_red_max = pts[0].m_red;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const auto& p = pts[i];
        if (p.F < s.F_min) { s.F_min = p.F; s.L_at_F_min = p.L; }
        if (p.F > s.F_max) { s.F_max = p.F; s.L_at_F_max = p.L; }
        if (p.p > s.p_max) { s.p_max = p.p; s.L_at_p_max = p.L; }
        fabsmax = std::max(fabsmax, std::fabs(p.F));
        dfd = std::max(dfd, std::fabs(p.F - p.F_fd));
        deq = std::max(deq, std::fabs(p.F - p.F_eq));
        s.f_min = std::min(s.f_min, p.f_n);
        s.f_max = std::max(s.f_max, p.f_n);
        s.m_red_min = std::min(s.m_red_min, p.m_red);
        s.m_red_max = std::max(s.m_red_max, p.m_red);
    }
    // energy balance: composite Simpson rule on the uniform grid (trapezoid for an even count)
    const std::size_t n = pts.size();
    const double h = (pts.back().L - pts.front().L) / static_cast<double>(n - 1) / 1000.0;
    if (n >= 3 && (n - 1) % 2 == 0) {
        double acc = pts.front().F + pts.back().F;
        for (std::size_t i = 1; i + 1 < n; ++i) acc += (i % 2 ? 4.0 : 2.0) * pts[i].F;
        s.work_J = acc * h / 3.0;
    } else {
        for (std::size_t i = 1; i < n; ++i) s.work_J += 0.5 * (pts[i].F + pts[i - 1].F) * h;
    }
    s.dV_J = pts.back().V - pts.front().V;
    if (rk != o.actuator) s.work_J = NAN;  // the energy balance holds for the swept actuator only
    s.dev_fd = fabsmax > 0 ? dfd / fabsmax : 0;
    s.dev_eq = fabsmax > 0 ? deq / fabsmax : 0;
    if (act.has_cylinder) {
        s.buckling_Pcr_Lmax = euler_load(act.cyl, act.cyl.lmax);
        s.buckling_n = s.F_max > 0 ? s.buckling_Pcr_Lmax / s.F_max : INFINITY;  // never pushed: no buckling
    }
    return s;
}

double peak_pressure(const Cylinder& c, const std::vector<SweepPoint>& pts, double flow_lpm, int dir, double t_r) {
    const double Q = flow_lpm / 60000.0;  // m^3/s
    const double A = (dir >= 0 ? c.area_piston() : c.area_annulus()) * 1e-6;
    const double v = Q / A;
    double pk = 0.0;
    for (const auto& p : pts) {
        const double w = 2.0 * kPi * p.f_n;
        const double dF = p.m_red * v * w * ramp_factor(w * t_r);
        pk = std::max({pk, pressure_bar(c, p.F + dF), pressure_bar(c, p.F - dF)});
    }
    return pk;
}

double min_ramp_time(const Cylinder& c, const std::vector<SweepPoint>& pts, double flow_lpm, int dir, double p_lim) {
    double lo = 0.0, hi = 1e-3;
    if (peak_pressure(c, pts, flow_lpm, dir, 1e300) > p_lim) return INFINITY;  // static load alone too high
    if (peak_pressure(c, pts, flow_lpm, dir, lo) <= p_lim) return 0.0;
    while (peak_pressure(c, pts, flow_lpm, dir, hi) > p_lim) hi *= 2.0;
    for (int i = 0; i < 60; ++i) {
        const double mid = 0.5 * (lo + hi);
        (peak_pressure(c, pts, flow_lpm, dir, mid) > p_lim ? lo : hi) = mid;
    }
    return hi;
}

}  // namespace cylload
