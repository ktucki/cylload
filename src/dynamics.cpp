// SPDX-License-Identifier: MIT
#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "cylload/analysis.hpp"

namespace cylload {

namespace {
constexpr double kPi = 3.14159265358979323846;

struct Table {
    std::vector<double> L, F, m;
    double lo = 0, h = 1;
    void at(double x, double& Fx, double& mx, double& dm) const {  // x in mm
        const double u = std::min(std::max((x - lo) / h, 0.0), static_cast<double>(L.size() - 1) - 1e-9);
        const std::size_t i = static_cast<std::size_t>(u);
        const double w = u - static_cast<double>(i);
        Fx = F[i] + w * (F[i + 1] - F[i]);
        mx = m[i] + w * (m[i + 1] - m[i]);
        dm = (m[i + 1] - m[i]) / h * 1000.0;  // kg per m
    }
};
}  // namespace

SimResult simulate(const Kinematics& kin, const SimOptions& o, int stride) {
    const auto& m = kin.model();
    const auto& a = m.actuators.at(static_cast<std::size_t>(o.actuator));
    if (!a.has_cylinder) throw std::runtime_error("simulate: actuator " + a.name + " has no cylinder data");
    const Cylinder& c = a.cyl;
    const double L0 = o.fixed.at(static_cast<std::size_t>(o.actuator));
    if (L0 < c.lmin || L0 > c.lmax) throw std::runtime_error("simulate: start length outside the stroke");

    // tabulate static force and reduced mass over the stroke
    Table tb;
    tb.lo = c.lmin;
    tb.h = c.stroke() / (o.table_n - 1);
    auto Ls = o.fixed;
    Pose s = kin.reference();
    for (int i = 0; i < o.table_n; ++i) {
        Ls[static_cast<std::size_t>(o.actuator)] = c.lmin + tb.h * i;
        s = kin.solve(s, Ls);
        tb.L.push_back(Ls[static_cast<std::size_t>(o.actuator)]);
        tb.F.push_back(forces_virtual_work(kin, s)[static_cast<std::size_t>(o.actuator)]);
        tb.m.push_back(reduced_mass(kin, s, o.actuator));
    }
    const double A1 = c.area_piston() * 1e-6, A2 = c.area_annulus() * 1e-6, stroke = c.stroke() / 1000.0;
    const double beta = m.hydraulics.beta, V0 = m.hydraulics.dead_volume * 1e-3, lmin = c.lmin / 1000.0;
    const double Q = o.flow_lpm / 60000.0;
    const double v = o.dir >= 0 ? Q / A1 : Q / A2;

    SimResult r;
    double F0, m0, dm0;
    tb.at(L0, F0, m0, dm0);
    const double k0 = hydraulic_stiffness(c, m.hydraulics, L0);
    const double w0 = std::sqrt(k0 / m0);
    r.m0 = m0;
    r.f0 = w0 / (2 * kPi);
    r.v = v;
    r.dF_formula = m0 * v * w0 * ramp_factor(w0 * o.t_ramp);
    const double cdamp = 2.0 * o.damping_ratio * std::sqrt(k0 * m0);
    const double t_end = o.t_end > 0 ? o.t_end : o.t_ramp + 6.0 * 2.0 * kPi / w0;

    // state: x [m], xd [m/s], p1, p2 [Pa]
    double y[4] = {L0 / 1000.0, 0.0, 0.0, 0.0};
    if (F0 >= 0) { y[3] = o.p_low; y[2] = (F0 + y[3] * A2) / A1; }
    else { y[2] = o.p_low; y[3] = (y[2] * A1 - F0) / A2; }

    auto deriv = [&](double t, const double* q, double* d) {
        const double ramp = o.t_ramp > 0 ? std::min(t / o.t_ramp, 1.0) : 1.0;
        const double qv = o.dir * v * ramp;
        double Fs, mr, dm;
        tb.at(q[0] * 1000.0, Fs, mr, dm);
        const double z = q[0] - lmin;
        const double V1 = V0 + A1 * z, V2 = V0 + A2 * (stroke - z);
        const double Fh = q[2] * A1 - q[3] * A2;
        d[0] = q[1];
        d[1] = (Fh - Fs - cdamp * q[1] - 0.5 * dm * q[1] * q[1]) / mr;
        d[2] = beta / V1 * (A1 * qv - A1 * q[1]);
        d[3] = beta / V2 * (A2 * q[1] - A2 * qv);
    };
    const int nsteps = static_cast<int>(std::ceil(t_end / o.dt));
    for (int i = 0; i <= nsteps; ++i) {
        const double t = i * o.dt;
        double Fs, mr, dm;
        tb.at(y[0] * 1000.0, Fs, mr, dm);
        const double dF = y[2] * A1 - y[3] * A2 - Fs;
        r.dF_max = std::max(r.dF_max, std::fabs(dF));
        r.p_peak = std::max({r.p_peak, y[2] / 1e5, y[3] / 1e5});
        if (i % stride == 0) r.history.push_back({t, y[0] * 1000.0, y[1] * 1000.0, y[2] / 1e5, y[3] / 1e5, dF});
        if (y[0] * 1000.0 <= c.lmin || y[0] * 1000.0 >= c.lmax) { r.hit_end = true; break; }
        double k1[4], k2[4], k3[4], k4[4], tmp[4];
        deriv(t, y, k1);
        for (int j = 0; j < 4; ++j) tmp[j] = y[j] + 0.5 * o.dt * k1[j];
        deriv(t + 0.5 * o.dt, tmp, k2);
        for (int j = 0; j < 4; ++j) tmp[j] = y[j] + 0.5 * o.dt * k2[j];
        deriv(t + 0.5 * o.dt, tmp, k3);
        for (int j = 0; j < 4; ++j) tmp[j] = y[j] + o.dt * k3[j];
        deriv(t + o.dt, tmp, k4);
        for (int j = 0; j < 4; ++j) y[j] += o.dt / 6.0 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
    }
    return r;
}

}  // namespace cylload

namespace cylload {

std::vector<MapCell> load_map(const Kinematics& kin, const MapOptions& o) {
    const auto& m = kin.model();
    const auto& rep = m.actuators.at(static_cast<std::size_t>(o.report));
    if (o.nx < 1 || o.ny < 1) throw std::invalid_argument("map: nx and ny must be positive");
    std::vector<MapCell> out;
    auto L = o.fixed;
    Pose row_start = kin.reference();
    for (int j = 0; j < o.ny; ++j) {
        L[static_cast<std::size_t>(o.ay)] = o.ny > 1 ? o.y0 + (o.y1 - o.y0) * j / (o.ny - 1) : o.y0;
        L[static_cast<std::size_t>(o.ax)] = o.x0;
        row_start = kin.solve(row_start, L);
        Pose s = row_start;
        for (int i = 0; i < o.nx; ++i) {
            L[static_cast<std::size_t>(o.ax)] = o.nx > 1 ? o.x0 + (o.x1 - o.x0) * i / (o.nx - 1) : o.x0;
            s = kin.solve(s, L);
            MapCell c;
            c.Lx = L[static_cast<std::size_t>(o.ax)];
            c.Ly = L[static_cast<std::size_t>(o.ay)];
            for (int b = 1; b <= m.n_moving(); ++b) c.angle_deg.push_back(kin.angle(s, b) * 180.0 / 3.14159265358979323846);
            c.F = forces_virtual_work(kin, s)[static_cast<std::size_t>(o.report)];
            c.m_red = reduced_mass(kin, s, o.report);
            if (rep.has_cylinder) {
                const double Lr = L[static_cast<std::size_t>(o.report)];
                c.p = pressure_bar(rep.cyl, c.F);
                c.f_n = std::sqrt(hydraulic_stiffness(rep.cyl, m.hydraulics, Lr) / c.m_red) / (2 * 3.14159265358979323846);
                if (o.flow_lpm > 0) {
                    SweepPoint sp;
                    sp.L = Lr;
                    sp.F = c.F;
                    sp.m_red = c.m_red;
                    sp.f_n = c.f_n;
                    const std::vector<SweepPoint> one{sp};
                    c.tr_ext = min_ramp_time(rep.cyl, one, o.flow_lpm, +1, o.p_lim);
                    c.tr_ret = min_ramp_time(rep.cyl, one, o.flow_lpm, -1, o.p_lim);
                }
            }
            out.push_back(c);
        }
    }
    return out;
}

}  // namespace cylload
