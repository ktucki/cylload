// SPDX-License-Identifier: MIT
#pragma once
#include <string>
#include <vector>

#include "cylload/model.hpp"

namespace cylload {

/// Joint reaction from the link-equilibrium solution.
struct Reaction {
    std::string joint;
    double fx = 0, fy = 0;  ///< force on body_b from body_a [N] (prismatic: normal force vector)
    double moment = 0;      ///< prismatic joints only: locking couple on body_b [N m]
};

struct StaticResult {
    std::vector<double> force_vw;   ///< virtual-work forces dV/dL [N], positive = pushing
    std::vector<double> force_eq;   ///< forces from the equilibrium of every link [N]
    std::vector<Reaction> reactions;
    double potential = 0;           ///< gravitational potential energy [J]
};

/// Gravitational potential energy [J] of all mass items in pose s.
double potential(const Kinematics& kin, const Pose& s);
/// Cylinder forces from the principle of virtual work, F_k = dV/dL_k (exact sensitivities).
std::vector<double> forces_virtual_work(const Kinematics& kin, const Pose& s);
/// Same forces by central differences of V with repeated position solutions (step h in mm).
std::vector<double> forces_finite_difference(const Kinematics& kin, const Pose& s, double h = 1e-3);
/// Independent route: force and moment balance of every body, unknown joint reactions
/// and actuator forces solved together. Does not use the constraint Jacobian.
StaticResult statics(const Kinematics& kin, const Pose& s);

/// Pressure [bar] needed for force F [N]: push -> piston area, pull -> annulus.
double pressure_bar(const Cylinder& c, double F);
/// Reduced (equivalent) mass seen by actuator k with the others locked [kg].
double reduced_mass(const Kinematics& kin, const Pose& s, int k);
/// Hydraulic stiffness of a cylinder with both chambers closed [N/m].
double hydraulic_stiffness(const Cylinder& c, const Hydraulics& h, double L_mm);
/// Valve-ramp dynamic factor kappa(tau), tau = omega * t_ramp.
double ramp_factor(double tau);
/// Euler critical load [N], pinned-pinned column of rod diameter and length L (conservative).
double euler_load(const Cylinder& c, double L_mm, double E_MPa = 210000.0);

struct SweepPoint {
    double L = 0;            ///< driven actuator length [mm]
    std::vector<double> angle_deg;  ///< body rotations from the reference configuration
    double F = 0, F_fd = 0, F_eq = 0;  ///< [N]
    double p = 0;            ///< [bar]
    double m_red = 0;        ///< [kg]
    double k_hyd = 0;        ///< [N/m]
    double f_n = 0;          ///< natural frequency [Hz]
    double V = 0;            ///< potential energy [J]
};

struct SweepOptions {
    int actuator = 0;           ///< actuator whose length is swept
    int report = -1;            ///< actuator whose loads are reported (-1: the swept one)
    std::vector<double> fixed;  ///< lengths of all actuators; entry `actuator` is overwritten
    double from = 0, to = 0;
    int n = 101;
    double fd_step = 1e-3;
};

struct SweepSummary {
    double F_min = 0, F_max = 0, L_at_F_min = 0, L_at_F_max = 0;
    double p_max = 0, L_at_p_max = 0;
    double work_J = 0, dV_J = 0;          ///< energy balance: integral F dL (Simpson) vs. Delta V (swept = reported only)
    double dev_fd = 0, dev_eq = 0;        ///< max |F - F_fd| / max|F|, max |F - F_eq| / max|F|
    double f_min = 0, f_max = 0, m_red_min = 0, m_red_max = 0;
    double buckling_Pcr_Lmax = 0, buckling_n = 0;  ///< Euler load at full extension, factor vs. max push force
};

std::vector<SweepPoint> sweep(const Kinematics& kin, const SweepOptions& o);
SweepSummary summarize(const Kinematics& kin, const SweepOptions& o, const std::vector<SweepPoint>& pts);

/// Peak pressure over a sweep for a valve opening ramp t_r [s] and flow Q [L/min]
/// (dir = +1 extension, -1 retraction); static force +/- m v omega kappa.
double peak_pressure(const Cylinder& c, const std::vector<SweepPoint>& pts, double flow_lpm, int dir, double t_r);
/// Shortest ramp time keeping the peak pressure below p_lim (bisection, kappa is monotone).
double min_ramp_time(const Cylinder& c, const std::vector<SweepPoint>& pts, double flow_lpm, int dir, double p_lim);

}  // namespace cylload

namespace cylload {

/// Time-domain start-up simulation of one cylinder (the other actuators locked):
/// Lagrange equation with the configuration-dependent reduced mass,
///   m(q) q'' + m'(q) q'^2 / 2 + dV/dq + c q' = p1 A1 - p2 A2,
/// and oil compressibility in both chambers fed by an ideal flow valve whose opening
/// rises linearly over the ramp time:
///   p1' = beta / V1 (A1 v r(t) - A1 q'),   p2' = beta / V2 (A2 q' - A2 v r(t)).
struct SimOptions {
    int actuator = 0;
    std::vector<double> fixed;   ///< lengths of all actuators; entry `actuator` = start length
    int dir = +1;                ///< +1 extension, -1 retraction
    double flow_lpm = 10.0;      ///< valve flow at full opening [L/min]
    double t_ramp = 0.1;         ///< opening ramp [s]
    double damping_ratio = 0.0;  ///< viscous damping as a fraction of critical (at the start)
    double dt = 2e-5;            ///< time step (RK4) [s]
    double t_end = 0.0;          ///< 0: ramp + 6 natural periods
    double p_low = 2e5;          ///< pressure of the unloaded chamber at rest [Pa]
    int table_n = 201;           ///< tabulation of F_st and m_red over the stroke
};

struct SimSample {
    double t, L, v, p1, p2, dF;  ///< s, mm, mm/s, bar, bar, N (hydraulic minus static force)
};

struct SimResult {
    std::vector<SimSample> history;  ///< every `stride`-th step
    double dF_max = 0;               ///< max |hydraulic - static force| [N]
    double dF_formula = 0;           ///< m v omega kappa(omega t_r) at the start position [N]
    double m0 = 0, f0 = 0, v = 0;    ///< reduced mass [kg], natural frequency [Hz], velocity [m/s] at the start
    double p_peak = 0;               ///< highest chamber pressure [bar]
    bool hit_end = false;            ///< the piston reached the end of the stroke
};

SimResult simulate(const Kinematics& kin, const SimOptions& o, int stride = 25);

}  // namespace cylload

namespace cylload {

/// Grid over the lengths of two actuators with the loads of a reported cylinder at every
/// node, including the local minimum valve-ramp time (condition applied at that node
/// only), e.g. tilt angle x hammer position for a ramp-time map in the controller.
struct MapOptions {
    int ax = 0, ay = 1;          ///< swept actuators (x fastest)
    int report = 0;              ///< actuator whose loads are mapped
    std::vector<double> fixed;   ///< lengths of all actuators
    double x0 = 0, x1 = 0, y0 = 0, y1 = 0;
    int nx = 21, ny = 21;
    double flow_lpm = 0.0;       ///< 0: no ramp times
    double p_lim = 140.0;
};

struct MapCell {
    double Lx = 0, Ly = 0;
    std::vector<double> angle_deg;
    double F = 0, p = 0, m_red = 0, f_n = 0;
    double tr_ext = 0, tr_ret = 0;  ///< local minimum ramp times [s] (inf: static load above the limit)
};

std::vector<MapCell> load_map(const Kinematics& kin, const MapOptions& o);

}  // namespace cylload
