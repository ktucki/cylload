// SPDX-License-Identifier: MIT
// cylload command-line interface
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "cylload/analysis.hpp"
#include "cylload/montecarlo.hpp"
#include "cylload/step.hpp"
#include "cylload/svg.hpp"
#include "cylload/version.hpp"

using namespace cylload;

namespace {

const char* kUsage =
    "cylload " CYLLOAD_VERSION " - hydraulic cylinder loads from a planar rigid-body model\n"
    "\n"
    "usage:\n"
    "  cylload info  MODEL\n"
    "  cylload pose  MODEL [ACT=L ...]\n"
    "  cylload sweep MODEL ACT [ACT=L ...] [--from L] [--to L] [--n N] [--report ACT] [--out FILE.csv]\n"
    "  cylload ramp  MODEL ACT [ACT=L ...] --flow LPM [--plim BAR] [--n N] [--out FILE.csv]\n"
    "  cylload simulate MODEL ACT [ACT=L ...] --flow LPM [--ramp S] [--dir ext|ret] [--start L]\n"
    "                [--damping ZETA] [--dt S] [--tend S] [--stride N] [--ramps T1,T2,...] [--out FILE.csv]\n"
    "  cylload map   MODEL ACT_X ACT_Y [ACT=L ...] [--report ACT] [--nx N] [--ny N] [--xfrom L --xto L]\n"
    "                [--yfrom L --yto L] [--flow LPM] [--plim BAR] [--out FILE.csv]\n"
    "  cylload mc    MODEL SPEC [--samples N] [--seed S] [--out FILE.csv] [--sensitivity FILE.csv] [--bootstrap B]\n"
    "  cylload draw  MODEL [ACT=L ...] [--overlay ACT=L,ACT=L ...] [--width MM] [--height MM] [--scale MM_PER_M [--ybottom MM]] [--color yes] [--cg yes] [--items no] [--labels no] [--font_pt PT] [--title TEXT]\n                [--step FILE --parts CSV [--level 2] [--plane xy] [--root NAME]] [--sketch FILE] [--view XMIN,YMIN,XMAX,YMAX]\n"
    "                [--out FILE.svg]\n"
    "  cylload plot  SPEC.txt [--data DIR] [--out FILE.svg]\n"
    "  cylload step  FILE.step [--root NAME] [--level K] [--prefix P ...] [--mass yes|no] [--bspline parametric|fan]\n"
    "                [--density KG_M3] [--out FILE.csv]\n"
    "  cylload masstable FILE.step ASSIGN.csv [--root NAME] [--level K] [--plane xy]\n"
    "                [--properties solid|outline] [--remainder KG] [--bspline parametric|fan] [--out MASSES.csv]\n"
    "\n"
    "Lengths in mm, masses in kg, forces in N, pressures in bar. Actuators that are not\n"
    "given (ACT=L) keep their length from the reference configuration of the model.\n";

struct Args {
    std::vector<std::string> pos;
    std::vector<std::pair<std::string, std::string>> opt;
    bool has(const std::string& k) const {
        for (const auto& o : opt) if (o.first == k) return true;
        return false;
    }
    std::string get(const std::string& k, const std::string& d = "") const {
        for (const auto& o : opt) if (o.first == k) return o.second;
        return d;
    }
    std::vector<std::string> all(const std::string& k) const {
        std::vector<std::string> v;
        for (const auto& o : opt) if (o.first == k) v.push_back(o.second);
        return v;
    }
};

Args parse_args(int argc, char** argv, int from) {
    Args a;
    for (int i = from; i < argc; ++i) {
        std::string s = argv[i];
        if (s.rfind("--", 0) == 0) {
            if (i + 1 >= argc) throw std::runtime_error("option " + s + " needs a value");
            a.opt.push_back({s.substr(2), argv[++i]});
        } else {
            a.pos.push_back(s);
        }
    }
    return a;
}

double to_d(const std::string& s) { return std::stod(s); }

/// actuator lengths: reference values overridden by NAME=L arguments
std::vector<double> lengths_from(const Model& m, const Kinematics& kin, const std::vector<std::string>& assigns) {
    auto L = kin.reference().lengths;
    for (const auto& [k, v] : parse_assignments(assigns)) L[m.actuator_index(k)] = v;
    return L;
}

std::ostream& outstream(const Args& a, std::ofstream& f) {
    if (!a.has("out")) return std::cout;
    f.open(a.get("out"));
    if (!f) throw std::runtime_error("cannot write " + a.get("out"));
    return f;
}

int cmd_info(const Args& a) {
    Model m = Model::from_file(a.pos.at(0));
    Kinematics kin(m);
    std::printf("model        %s\n", m.title.c_str());
    std::printf("bodies       %d moving (+ground), DOF %d, actuators %zu\n", m.n_moving(), m.dof(), m.actuators.size());
    std::printf("joints       %zu revolute, %zu prismatic\n", m.revolutes.size(), m.prismatics.size());
    std::printf("mass items   %zu, total %.3f kg (incl. cylinders)\n", m.all_masses().size(), m.total_mass());
    for (int b = 1; b <= m.n_moving(); ++b) {
        double mb = 0;
        for (const auto& mi : m.all_masses()) if (mi.body == b) mb += mi.m;
        std::printf("  body %-10s %10.3f kg\n", m.bodies[b].c_str(), mb);
    }
    for (const auto& ef : m.forces)
        std::printf("force    %-10s on %-10s at (%.1f, %.1f) mm: (%.1f, %.1f) N\n", ef.name.c_str(), m.bodies[ef.body].c_str(), ef.p[0], ef.p[1], ef.f[0], ef.f[1]);
    for (const auto& ac : m.actuators) {
        std::printf("actuator %-6s %s -> %s  L_ref %.3f mm", ac.name.c_str(), m.bodies[ac.body_a].c_str(), m.bodies[ac.body_b].c_str(), ac.l_ref);
        if (ac.has_cylinder)
            std::printf("  cylinder %g/%g, %g..%g mm, %g kg", ac.cyl.bore, ac.cyl.rod, ac.cyl.lmin, ac.cyl.lmax, ac.cyl.mass);
        std::printf("\n");
    }
    return 0;
}

int cmd_pose(const Args& a) {
    Model m = Model::from_file(a.pos.at(0));
    Kinematics kin(m);
    auto L = lengths_from(m, kin, std::vector<std::string>(a.pos.begin() + 1, a.pos.end()));
    Pose s = kin.solve(kin.reference(), L);
    auto st = statics(kin, s);
    auto fd = forces_finite_difference(kin, s);
    std::printf("%-8s %10s %13s %13s %13s %9s\n", "actuator", "L_mm", "F_vw_N", "F_eq_N", "F_fd_N", "p_bar");
    for (std::size_t k = 0; k < m.actuators.size(); ++k) {
        const auto& ac = m.actuators[k];
        std::printf("%-8s %10.3f %13.3f %13.3f %13.3f", ac.name.c_str(), L[k], st.force_vw[k], st.force_eq[k], fd[k]);
        if (ac.has_cylinder) std::printf(" %9.2f", pressure_bar(ac.cyl, st.force_vw[k]));
        std::printf("\n");
    }
    std::printf("\n%-14s %12s %12s %12s %10s\n", "joint", "Fx_N", "Fy_N", "|F|_N", "M_Nm");
    for (const auto& r : st.reactions)
        std::printf("%-14s %12.1f %12.1f %12.1f %10.1f\n", r.joint.c_str(), r.fx, r.fy, std::hypot(r.fx, r.fy), r.moment);
    std::printf("\nbody rotations [deg]:");
    for (int b = 1; b <= m.n_moving(); ++b) std::printf(" %s=%.3f", m.bodies[b].c_str(), kin.angle(s, b) * 180.0 / 3.14159265358979323846);
    std::printf("\npotential energy %.3f J\n", st.potential);
    return 0;
}

SweepOptions sweep_options(const Model& m, const Kinematics& kin, const Args& a) {
    SweepOptions o;
    o.actuator = m.actuator_index(a.pos.at(1));
    o.fixed = lengths_from(m, kin, std::vector<std::string>(a.pos.begin() + 2, a.pos.end()));
    const auto& ac = m.actuators[o.actuator];
    o.from = a.has("from") ? to_d(a.get("from")) : (ac.has_cylinder ? ac.cyl.lmin : ac.l_ref);
    o.to = a.has("to") ? to_d(a.get("to")) : (ac.has_cylinder ? ac.cyl.lmax : ac.l_ref);
    o.n = std::stoi(a.get("n", "101"));
    if (a.has("report")) o.report = m.actuator_index(a.get("report"));
    if (o.from == o.to) throw std::runtime_error("sweep: give --from and --to for an actuator without cylinder data");
    return o;
}

void print_summary(const Model& m, const SweepOptions& o, const SweepSummary& s) {
    const auto& ac = m.actuators[o.report < 0 ? o.actuator : o.report];
    std::fprintf(stderr, "actuator %s, L %.2f..%.2f mm, %d points\n", m.actuators[o.actuator].name.c_str(), o.from, o.to, o.n);
    if (o.report >= 0 && o.report != o.actuator) std::fprintf(stderr, "  reported loads: actuator %s\n", ac.name.c_str());
    std::fprintf(stderr, "  F_max %.1f N at L=%.2f   F_min %.1f N at L=%.2f\n", s.F_max, s.L_at_F_max, s.F_min, s.L_at_F_min);
    if (ac.has_cylinder) std::fprintf(stderr, "  p_max %.2f bar at L=%.2f\n", s.p_max, s.L_at_p_max);
    std::fprintf(stderr, "  checks: |F_vw-F_fd|/max|F| = %.2e, |F_vw-F_eq|/max|F| = %.2e\n", s.dev_fd, s.dev_eq);
    if (!std::isnan(s.work_J))
        std::fprintf(stderr, "  energy: integral F dL = %.3f J, Delta V = %.3f J (rel. diff %.2e)\n", s.work_J, s.dV_J,
                     s.dV_J != 0 ? std::fabs(s.work_J - s.dV_J) / std::fabs(s.dV_J) : 0.0);
    std::fprintf(stderr, "  reduced mass %.1f..%.1f kg", s.m_red_min, s.m_red_max);
    if (ac.has_cylinder)
        std::isfinite(s.buckling_n)
            ? std::fprintf(stderr, ", natural frequency %.2f..%.2f Hz\n  Euler load at L_max %.0f N, safety factor %.2f\n", s.f_min, s.f_max,
                           s.buckling_Pcr_Lmax, s.buckling_n)
            : std::fprintf(stderr, ", natural frequency %.2f..%.2f Hz\n  Euler load at L_max %.0f N, not loaded in compression\n", s.f_min,
                           s.f_max, s.buckling_Pcr_Lmax);
    else
        std::fprintf(stderr, "\n");
}

void write_sweep_csv(std::ostream& os, const Model& m, const std::vector<SweepPoint>& pts) {
    os << "L_mm";
    for (int b = 1; b <= m.n_moving(); ++b) os << ",rot_" << m.bodies[b] << "_deg";
    os << ",F_N,F_fd_N,F_eq_N,p_bar,m_red_kg,k_hyd_N_per_m,f_n_Hz,V_J\n";
    os.precision(10);
    for (const auto& p : pts) {
        os << p.L;
        for (double x : p.angle_deg) os << ',' << x;
        os << ',' << p.F << ',' << p.F_fd << ',' << p.F_eq << ',' << p.p << ',' << p.m_red << ',' << p.k_hyd << ',' << p.f_n << ','
           << p.V << '\n';
    }
}

int cmd_sweep(const Args& a) {
    Model m = Model::from_file(a.pos.at(0));
    Kinematics kin(m);
    auto o = sweep_options(m, kin, a);
    auto pts = sweep(kin, o);
    std::ofstream f;
    write_sweep_csv(outstream(a, f), m, pts);
    print_summary(m, o, summarize(kin, o, pts));
    return 0;
}

int cmd_ramp(const Args& a) {
    Model m = Model::from_file(a.pos.at(0));
    Kinematics kin(m);
    auto o = sweep_options(m, kin, a);
    if (o.report >= 0 && o.report != o.actuator) throw std::runtime_error("ramp: --report is not supported");
    const auto& ac = m.actuators[o.actuator];
    if (!ac.has_cylinder) throw std::runtime_error("ramp: actuator " + ac.name + " has no cylinder data");
    if (!a.has("flow")) throw std::runtime_error("ramp: --flow LPM is required");
    const double Q = to_d(a.get("flow")), plim = to_d(a.get("plim", "140"));
    auto pts = sweep(kin, o);
    const double s = ac.cyl.stroke() / 1000.0;
    const double vext = Q / 60000.0 / (ac.cyl.area_piston() * 1e-6), vret = Q / 60000.0 / (ac.cyl.area_annulus() * 1e-6);
    std::fprintf(stderr, "actuator %s, flow %.2f L/min: v_ext %.1f mm/s (%.2f s), v_ret %.1f mm/s (%.2f s)\n", ac.name.c_str(), Q,
                 vext * 1000, s / vext, vret * 1000, s / vret);
    const double t_ext = min_ramp_time(ac.cyl, pts, Q, +1, plim), t_ret = min_ramp_time(ac.cyl, pts, Q, -1, plim);
    std::fprintf(stderr, "minimum valve ramp for p <= %.1f bar: extension %.3f s, retraction %.3f s\n", plim, t_ext, t_ret);
    std::ofstream f;
    std::ostream& os = outstream(a, f);
    os << "t_ramp_s,p_peak_ext_bar,p_peak_ret_bar\n";
    for (double t : {0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0})
        os << t << ',' << peak_pressure(ac.cyl, pts, Q, +1, t) << ',' << peak_pressure(ac.cyl, pts, Q, -1, t) << '\n';
    return 0;
}

int cmd_mc(const Args& a) {
    Model m = Model::from_file(a.pos.at(0));
    McSpec sp = McSpec::from_file(a.pos.at(1));
    if (a.has("samples")) sp.samples = std::stoi(a.get("samples"));
    if (a.has("seed")) sp.seed = std::stoull(a.get("seed"));
    if (a.has("bootstrap")) sp.bootstrap = std::stoi(a.get("bootstrap"));
    auto r = run_monte_carlo(m, sp);
    std::fprintf(stderr, "Monte Carlo: %d samples, seed %llu\n", sp.samples, static_cast<unsigned long long>(sp.seed));
    std::fprintf(stderr, "%-14s %9s %9s %9s %9s %9s   (pressure, bar)\n", "output", "min", "p5", "p50", "p95", "max");
    for (std::size_t i = 0; i < r.output_labels.size(); ++i) {
        const auto& s = r.pressure_stats[i];
        std::fprintf(stderr, "%-14s %9.2f %9.2f %9.2f %9.2f %9.2f", r.output_labels[i].c_str(), s.min, s.p5, s.p50, s.p95, s.max);
        if (sp.bootstrap > 0)
            std::fprintf(stderr, "   95%% CI: p5 %.2f..%.2f, p95 %.2f..%.2f", s.p5_lo, s.p5_hi, s.p95_lo, s.p95_hi);
        std::fprintf(stderr, "\n");
    }
    if (sp.bootstrap > 0) std::fprintf(stderr, "(confidence intervals from %d bootstrap resamples)\n", sp.bootstrap);
    std::fprintf(stderr, "\nSpearman rank correlation (pressure vs. input):\n%-14s", "input");
    for (const auto& o : r.output_labels) std::fprintf(stderr, " %12s", o.c_str());
    std::fprintf(stderr, "\n");
    for (std::size_t j = 0; j < r.input_labels.size(); ++j) {
        std::fprintf(stderr, "%-14s", r.input_labels[j].c_str());
        for (std::size_t i = 0; i < r.output_labels.size(); ++i) std::fprintf(stderr, " %12.3f", r.spearman[i][j]);
        std::fprintf(stderr, "\n");
    }
    std::fprintf(stderr, "\nFirst-order indices S_i / standardised regression coefficients (R^2 in the last row):\n%-14s", "input");
    for (const auto& o : r.output_labels) std::fprintf(stderr, " %17s", o.c_str());
    std::fprintf(stderr, "\n");
    for (std::size_t j = 0; j < r.input_labels.size(); ++j) {
        std::fprintf(stderr, "%-14s", r.input_labels[j].c_str());
        for (std::size_t i = 0; i < r.output_labels.size(); ++i) std::fprintf(stderr, "    %6.3f / %6.3f", r.first_order[i][j], r.src[i][j]);
        std::fprintf(stderr, "\n");
    }
    if (sp.bootstrap > 0) {
        std::fprintf(stderr, "\n95%% bootstrap confidence intervals of S_i:\n%-14s", "input");
        for (const auto& o : r.output_labels) std::fprintf(stderr, " %17s", o.c_str());
        std::fprintf(stderr, "\n");
        for (std::size_t j = 0; j < r.input_labels.size(); ++j) {
            std::fprintf(stderr, "%-14s", r.input_labels[j].c_str());
            for (std::size_t i = 0; i < r.output_labels.size(); ++i)
                std::fprintf(stderr, "     %5.3f..%5.3f", r.first_order_lo[i][j], r.first_order_hi[i][j]);
            std::fprintf(stderr, "\n");
        }
    }
    std::fprintf(stderr, "%-14s", "R^2");
    for (std::size_t i = 0; i < r.output_labels.size(); ++i) std::fprintf(stderr, " %17.3f", r.r2[i]);
    std::fprintf(stderr, "\n");
    if (a.has("sensitivity")) {
        std::ofstream sf(a.get("sensitivity"));
        if (!sf) throw std::runtime_error("cannot write " + a.get("sensitivity"));
        sf.precision(8);
        sf << "output,input,spearman,first_order,first_order_ci_lo,first_order_ci_hi,src,r2\n";
        for (std::size_t i = 0; i < r.output_labels.size(); ++i)
            for (std::size_t j = 0; j < r.input_labels.size(); ++j)
                sf << r.output_labels[i] << ',' << r.input_labels[j] << ',' << r.spearman[i][j] << ',' << r.first_order[i][j] << ','
                   << r.first_order_lo[i][j] << ',' << r.first_order_hi[i][j] << ',' << r.src[i][j] << ',' << r.r2[i] << '\n';
        sf << "#stats,output,min,p5,p50,p95,max,mean,p5_ci_lo,p5_ci_hi,p95_ci_lo,p95_ci_hi\n";
        for (std::size_t i = 0; i < r.output_labels.size(); ++i) {
            const auto& st = r.pressure_stats[i];
            sf << "#stats," << r.output_labels[i] << ',' << st.min << ',' << st.p5 << ',' << st.p50 << ',' << st.p95 << ',' << st.max << ','
               << st.mean << ',' << st.p5_lo << ',' << st.p5_hi << ',' << st.p95_lo << ',' << st.p95_hi << '\n';
        }
    }
    std::ofstream f;
    std::ostream& os = outstream(a, f);
    os.precision(10);
    os << "sample";
    for (const auto& l : r.input_labels) os << ',' << l;
    for (const auto& l : r.output_labels) os << ',' << l << "_F_N," << l << "_p_bar";
    os << '\n';
    for (std::size_t n = 0; n < r.inputs.size(); ++n) {
        os << n;
        for (double x : r.inputs[n]) os << ',' << x;
        for (std::size_t i = 0; i < r.output_labels.size(); ++i) os << ',' << r.force[n][i] << ',' << r.pressure[n][i];
        os << '\n';
    }
    return 0;
}

int cmd_step(const Args& a) {
    const std::size_t level = std::stoul(a.get("level", "1"));
    const bool mass = a.get("mass", "yes") != "no";
    const double rho = to_d(a.get("density", "7850"));  // kg/m^3, only for the reference steel mass column
    set_bspline_integration(a.get("bspline", "parametric") != "fan");
    StepModel sm = read_step(a.pos.at(0), a.get("root"), mass);
    std::size_t nf = 0, nt = 0, np = 0;
    double area = 0, area_t = 0;
    for (const auto& o : sm.occurrences) { nf += o.faces; nt += o.faces_triangulated; np += o.faces_parametric; area += o.area; area_t += o.area_triangulated; }
    std::fprintf(stderr, "root '%s': %zu entities, %zu tree nodes, %zu with geometry, length unit %s", sm.root.c_str(), sm.n_entities,
                 sm.n_tree_nodes, sm.occurrences.size(), sm.length_unit.c_str());
    if (mass)
        std::fprintf(stderr, "; %zu faces integrated (%zu B-spline faces in their parameter domain, %zu triangulated = %.3g%% of the area)",
                     nf, np, nt, area > 0 ? 100.0 * area_t / area : 0.0);
    std::fprintf(stderr, "\n");
    auto groups = group_by_level(sm, level);
    std::ofstream f;
    std::ostream& os = outstream(a, f);
    os.precision(10);
    auto row = [&](const StepGroup& g) {
        os << '"' << g.key << "\"," << g.n_occurrences << ',' << g.n_points;
        for (int i = 0; i < 3; ++i) os << ',' << g.lo[i];
        for (int i = 0; i < 3; ++i) os << ',' << g.hi[i];
        for (int i = 0; i < 3; ++i) os << ',' << (g.lo[i] + g.hi[i]) / 2;
        if (mass) {
            os << ',' << g.mass.volume << ',' << g.mass.volume * rho * 1e-9;
            for (int i = 0; i < 3; ++i) os << ',' << g.mass.centroid[i];
            os << ',' << g.mass.k2(0, 1) << ',' << g.faces << ',' << g.faces_parametric << ',' << g.faces_triangulated << ',' << g.area << ','
               << g.area_triangulated << ',' << g.volume_check << ',' << g.moment_check;
        }
        os << '\n';
    };
    os << "name,occurrences,points,xmin,ymin,zmin,xmax,ymax,zmax,xc,yc,zc";
    if (mass) os << ",volume_mm3,mass_at_density_kg,xg,yg,zg,kz2_m2,faces,faces_bspline_parametric,faces_triangulated,area_mm2,area_triangulated_mm2,volume_check,moment_check";
    os << '\n';
    auto prefixes = a.all("prefix");
    if (prefixes.empty()) {
        for (const auto& g : groups) row(g);
    } else {
        for (const auto& p : prefixes) row(union_by_prefix(groups, p));
    }
    return 0;
}

int cmd_simulate(const Args& a) {
    Model m = Model::from_file(a.pos.at(0));
    Kinematics kin(m);
    SimOptions o;
    o.actuator = m.actuator_index(a.pos.at(1));
    o.fixed = lengths_from(m, kin, std::vector<std::string>(a.pos.begin() + 2, a.pos.end()));
    const auto& c = m.actuators[static_cast<std::size_t>(o.actuator)].cyl;
    if (a.has("start")) o.fixed[static_cast<std::size_t>(o.actuator)] = to_d(a.get("start"));
    else if (a.pos.size() <= 2 || parse_assignments(std::vector<std::string>(a.pos.begin() + 2, a.pos.end())).count(a.pos.at(1)) == 0)
        o.fixed[static_cast<std::size_t>(o.actuator)] = 0.5 * (c.lmin + c.lmax);
    if (!a.has("flow")) throw std::runtime_error("simulate: --flow LPM is required");
    o.flow_lpm = to_d(a.get("flow"));
    o.t_ramp = to_d(a.get("ramp", "0.1"));
    o.dir = a.get("dir", "ext") == "ret" ? -1 : +1;
    o.damping_ratio = to_d(a.get("damping", "0"));
    o.dt = to_d(a.get("dt", "2e-5"));
    o.t_end = to_d(a.get("tend", "0"));
    if (a.has("ramps")) {  // table of excess forces for a list of ramp times
        std::vector<double> trs;
        std::string cur;
        for (char ch : a.get("ramps") + ",") {
            if (ch == ',') { if (!cur.empty()) trs.push_back(to_d(cur)); cur.clear(); } else cur += ch;
        }
        std::ofstream f;
        std::ostream& os = outstream(a, f);
        os.precision(10);
        os << "t_ramp_s,tau,kappa_formula,kappa_sim,dF_sim_N,dF_formula_N,dF_sim_less_damping_drag_N\n";
        for (double tr : trs) {
            o.t_ramp = tr;
            auto r = simulate(kin, o, 1000000);
            const double w = 2 * 3.14159265358979323846 * r.f0, base = r.m0 * r.v * w;
            os << tr << ',' << w * tr << ',' << ramp_factor(w * tr) << ',' << r.dF_max / base << ',' << r.dF_max << ',' << r.dF_formula << ','
               << r.dF_max_osc << '\n';
            std::fprintf(stderr, "t_r %.4f s: simulation %.1f N, formula %.1f N (%+.2f %%)\n", tr, r.dF_max, r.dF_formula,
                         100.0 * (r.dF_max / r.dF_formula - 1));
        }
        return 0;
    }
    auto r = simulate(kin, o, std::stoi(a.get("stride", "25")));
    std::fprintf(stderr, "start L=%.2f mm, %s at %.2f L/min (v=%.1f mm/s), ramp %.3f s, damping %.3f\n",
                 o.fixed[static_cast<std::size_t>(o.actuator)], o.dir > 0 ? "extension" : "retraction", o.flow_lpm, r.v * 1000,
                 o.t_ramp, o.damping_ratio);
    std::fprintf(stderr, "reduced mass %.1f kg, natural frequency %.3f Hz\n", r.m0, r.f0);
    std::fprintf(stderr, "max excess force: simulation %.1f N, formula m v w kappa %.1f N (%+.2f %%)\n", r.dF_max, r.dF_formula,
                 r.dF_formula > 0 ? 100.0 * (r.dF_max / r.dF_formula - 1) : 0.0);
    std::fprintf(stderr, "peak chamber pressure %.1f bar%s\n", r.p_peak, r.hit_end ? "; the piston reached the end of the stroke" : "");
    std::ofstream f;
    std::ostream& os = outstream(a, f);
    os.precision(10);
    os << "t_s,L_mm,v_mm_s,p1_bar,p2_bar,dF_N\n";
    for (const auto& h : r.history) os << h.t << ',' << h.L << ',' << h.v << ',' << h.p1 << ',' << h.p2 << ',' << h.dF << '\n';
    return 0;
}

int cmd_map(const Args& a) {
    Model m = Model::from_file(a.pos.at(0));
    Kinematics kin(m);
    MapOptions o;
    o.ax = m.actuator_index(a.pos.at(1));
    o.ay = m.actuator_index(a.pos.at(2));
    o.report = m.actuator_index(a.get("report", a.pos.at(1)));
    o.fixed = lengths_from(m, kin, std::vector<std::string>(a.pos.begin() + 3, a.pos.end()));
    auto range = [&](int k, const std::string& lo, const std::string& hi, double& r0, double& r1) {
        const auto& ac = m.actuators[static_cast<std::size_t>(k)];
        r0 = a.has(lo) ? to_d(a.get(lo)) : (ac.has_cylinder ? ac.cyl.lmin : ac.l_ref);
        r1 = a.has(hi) ? to_d(a.get(hi)) : (ac.has_cylinder ? ac.cyl.lmax : ac.l_ref);
    };
    range(o.ax, "xfrom", "xto", o.x0, o.x1);
    range(o.ay, "yfrom", "yto", o.y0, o.y1);
    o.nx = std::stoi(a.get("nx", "41"));
    o.ny = std::stoi(a.get("ny", "21"));
    o.flow_lpm = to_d(a.get("flow", "0"));
    o.p_lim = to_d(a.get("plim", "140"));
    auto cells = load_map(kin, o);
    double pmax = 0, tmax = 0;
    for (const auto& c : cells) {
        pmax = std::max(pmax, c.p);
        if (std::isfinite(c.tr_ret)) tmax = std::max(tmax, c.tr_ret);
    }
    std::fprintf(stderr, "map %s x %s (%d x %d), loads of %s: p_max %.2f bar", m.actuators[static_cast<std::size_t>(o.ax)].name.c_str(),
                 m.actuators[static_cast<std::size_t>(o.ay)].name.c_str(), o.nx, o.ny, m.actuators[static_cast<std::size_t>(o.report)].name.c_str(), pmax);
    if (o.flow_lpm > 0) std::fprintf(stderr, ", largest local ramp time (retraction) %.3f s", tmax);
    std::fprintf(stderr, "\n");
    std::ofstream f;
    std::ostream& os = outstream(a, f);
    os.precision(10);
    os << "Lx_mm,Ly_mm";
    for (int b = 1; b <= m.n_moving(); ++b) os << ",rot_" << m.bodies[static_cast<std::size_t>(b)] << "_deg";
    os << ",F_N,p_bar,m_red_kg,f_n_Hz,tr_ext_s,tr_ret_s\n";
    for (const auto& c : cells) {
        os << c.Lx << ',' << c.Ly;
        for (double x : c.angle_deg) os << ',' << x;
        os << ',' << c.F << ',' << c.p << ',' << c.m_red << ',' << c.f_n << ',' << c.tr_ext << ',' << c.tr_ret << '\n';
    }
    return 0;
}

std::string dir_of_path(const std::string& p) {
    auto k = p.find_last_of("/\\");
    return k == std::string::npos ? std::string(".") : p.substr(0, k);
}

int cmd_plot(const Args& a) {
    std::ifstream f(a.pos.at(0));
    if (!f) throw std::runtime_error("cannot open " + a.pos.at(0));
    std::stringstream ss;
    ss << f.rdbuf();
    Figure fig = figure_from_spec(ss.str(), a.has("data") ? a.get("data") : dir_of_path(a.pos.at(0)));
    std::ofstream o;
    outstream(a, o) << render_svg(fig);
    return 0;
}

std::vector<std::map<std::string, std::string>> read_csv(const std::string& path);

bool ends_with(const std::string& s, const std::string& suf) {
    return s.size() > suf.size() && s.compare(s.size() - suf.size(), suf.size(), suf) == 0;
}

int cmd_draw(const Args& a) {
    Model m = Model::from_file(a.pos.at(0));
    Kinematics kin(m);
    std::vector<Pose> poses;
    poses.push_back(kin.solve(kin.reference(), lengths_from(m, kin, std::vector<std::string>(a.pos.begin() + 1, a.pos.end()))));
    for (const auto& ov : a.all("overlay")) {  // --overlay "S8=698,S10=758"
        std::vector<std::string> as;
        std::string cur;
        for (char c : ov + ",") {
            if (c == ',') { if (!cur.empty()) as.push_back(cur); cur.clear(); } else cur += c;
        }
        auto L = poses.front().lengths;
        for (const auto& [k, v] : parse_assignments(as)) L[static_cast<std::size_t>(m.actuator_index(k))] = v;
        poses.push_back(kin.solve(poses.front(), L));
    }
    std::ofstream o;
    DrawOptions d;
    d.scale_mm_per_m = to_d(a.get("scale", "0"));
    d.width_mm = to_d(a.get("width", "120"));
    d.height_mm = to_d(a.get("height", d.scale_mm_per_m > 0 ? "0" : "100"));
    if (a.has("ybottom")) d.y_bottom = to_d(a.get("ybottom", "0"));
    d.labels = a.get("labels", "yes") != "no";
    d.colour = a.get("color", "no") == "yes";
    d.cg = a.get("cg", "no") == "yes";
    d.items = a.get("items", "yes") != "no";
    d.font_pt = to_d(a.get("font_pt", "7"));
    d.title = a.get("title", "");
    if (a.has("step")) {  // solid-model parts: --step FILE --parts CSV (prefix,body) [--level 2] [--plane xy]
        const std::size_t level = std::stoul(a.get("level", "2"));
        const std::string plane = a.get("plane", "xy");
        if (plane.size() != 2 || plane.find_first_not_of("xyz") != std::string::npos || plane[0] == plane[1])
            throw std::runtime_error("--plane must be xy, yz, zx, ...");
        const int ci = plane[0] - 'x', cj = plane[1] - 'x', ck = 3 - ci - cj;
        const double dsign = ((cj - ci + 3) % 3 == 1) ? 1.0 : -1.0;  // viewer on the side of i x j
        struct Rule { std::string prefix; int body = 0; int act = -1; bool rod = false; bool skip = false; };
        std::vector<Rule> rules;
        for (const auto& r : read_csv(a.get("parts", ""))) {
            Rule ru;
            ru.prefix = r.at("prefix");
            const std::string b = r.at("body");
            if (b == "-") ru.skip = true;
            else if (ends_with(b, ".rod") || ends_with(b, ".barrel")) {
                ru.rod = ends_with(b, ".rod");
                ru.act = m.actuator_index(b.substr(0, b.size() - (ru.rod ? 4 : 7)));
            } else if (b == "ground") ru.body = 0;
            else ru.body = m.body_index(b);
            rules.push_back(ru);
        }
        const StepModel sm = read_step(a.get("step", ""), a.get("root"), false, true);
        std::size_t used = 0, skipped = 0;
        for (const auto& oc : sm.occurrences) {
            const std::string& key = oc.path[std::min(level, oc.path.size() - 1)];
            const Rule* ru = nullptr;
            for (const auto& r : rules)
                if (key.rfind(r.prefix, 0) == 0) { ru = &r; break; }
            if (!ru || ru->skip) { ++skipped; continue; }
            for (const auto& so : oc.solids) {
                DrawPart dp;
                dp.body = ru->body;
                dp.actuator = ru->act;
                dp.rod = ru->rod;
                double zs = 0; std::size_t nz = 0;
                dp.curved = so.curved;
                for (const auto& f : so.faces) {
                    std::vector<std::vector<Vec2>> loops;
                    for (const auto& lp : f) {
                        std::vector<Vec2> l2;
                        for (const auto& q : lp) { l2.push_back({q[static_cast<std::size_t>(ci)], q[static_cast<std::size_t>(cj)]}); zs += q[static_cast<std::size_t>(ck)]; ++nz; }
                        loops.push_back(std::move(l2));
                    }
                    dp.faces.push_back(std::move(loops));
                }
                for (const auto& e : so.edges) {
                    std::vector<Vec2> e2;
                    for (const auto& q : e) e2.push_back({q[static_cast<std::size_t>(ci)], q[static_cast<std::size_t>(cj)]});
                    dp.edges.push_back(std::move(e2));
                }
                dp.depth = nz ? dsign * zs / static_cast<double>(nz) : 0.0;
                d.parts.push_back(std::move(dp));
            }
            ++used;
        }
        std::fprintf(stderr, "solid model: %zu occurrences drawn, %zu left out, %zu solids\n", used, skipped, d.parts.size());
    }
    if (a.has("sketch")) {  // 2D line work in the reference configuration: lines "BODY [fill] X1 Y1 X2 Y2 ..." (mm)
        std::ifstream f(a.get("sketch", ""));
        if (!f) throw std::runtime_error("cannot open " + a.get("sketch", ""));
        std::string line;
        std::size_t n = 0, ln = 0;
        while (std::getline(f, line)) {
            ++ln;
            const auto hsh = line.find('#');
            if (hsh != std::string::npos) line.erase(hsh);
            std::istringstream is(line);
            std::string b;
            if (!(is >> b)) continue;
            std::string tok;
            bool fill = false;
            std::vector<double> v;
            while (is >> tok) {
                if (tok == "fill") { fill = true; continue; }
                v.push_back(to_d(tok));
            }
            if (v.size() < 4 || v.size() % 2) throw std::runtime_error("sketch line " + std::to_string(ln) + ": BODY [fill] X1 Y1 X2 Y2 ...");
            DrawPart dp;
            if (ends_with(b, ".rod") || ends_with(b, ".barrel")) {
                dp.rod = ends_with(b, ".rod");
                dp.actuator = m.actuator_index(b.substr(0, b.size() - (dp.rod ? 4 : 7)));
            } else dp.body = b == "ground" ? 0 : m.body_index(b);
            std::vector<Vec2> pts;
            for (std::size_t i = 0; i + 1 < v.size(); i += 2) pts.push_back({v[i], v[i + 1]});
            if (fill) { dp.faces.push_back({pts}); dp.curved.push_back(0); }
            dp.edges.push_back(std::move(pts));
            dp.depth = static_cast<double>(ln);  // later lines are drawn on top
            d.parts.push_back(std::move(dp));
            ++n;
        }
        std::fprintf(stderr, "sketch: %zu lines\n", n);
    }
    if (a.has("view")) {  // --view XMIN,YMIN,XMAX,YMAX
        std::string v = a.get("view", ""), cur;
        for (char c : v + ",") {
            if (c == ',') { if (!cur.empty()) d.view.push_back(to_d(cur)); cur.clear(); } else cur += c;
        }
        if (d.view.size() != 4) throw std::runtime_error("--view needs XMIN,YMIN,XMAX,YMAX");
    }
    outstream(a, o) << draw_mechanism(kin, poses, d);
    return 0;
}

/// read a small CSV with a header into rows of name -> value
std::vector<std::map<std::string, std::string>> read_csv(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::vector<std::map<std::string, std::string>> rows;
    std::vector<std::string> head;
    std::string line;
    auto split = [](const std::string& l) {
        std::vector<std::string> out;
        std::string cur;
        bool q = false;
        for (char c : l) {
            if (c == '"') { q = !q; continue; }
            if (c == ',' && !q) { out.push_back(cur); cur.clear(); continue; }
            if (c != '\r') cur += c;
        }
        out.push_back(cur);
        for (auto& x : out) {
            auto a = x.find_first_not_of(" \t"), b = x.find_last_not_of(" \t");
            x = a == std::string::npos ? std::string() : x.substr(a, b - a + 1);
        }
        return out;
    };
    while (std::getline(f, line)) {
        if (!line.empty() && static_cast<unsigned char>(line[0]) == 0xEF && line.size() >= 3) line = line.substr(3);
        auto t = line.find_first_not_of(" \t\r");
        if (t == std::string::npos || line[t] == '#') continue;
        auto v = split(line);
        if (head.empty()) { head = v; continue; }
        std::map<std::string, std::string> r;
        for (std::size_t i = 0; i < head.size(); ++i) r[head[i]] = i < v.size() ? v[i] : std::string();
        rows.push_back(r);
    }
    return rows;
}

int cmd_masstable(const Args& a) {
    // assignment: name,body,group,mass_kg,prefix[,remainder][,x_mm,y_mm,J_kgm2 when prefix is empty]
    const std::size_t level = std::stoul(a.get("level", "2"));
    const std::string plane = a.get("plane", "xy");
    if (plane.size() != 2 || plane.find_first_not_of("xyz") != std::string::npos || plane[0] == plane[1])
        throw std::runtime_error("--plane must be xy, yz, zx, ...");
    const int ci = plane[0] - 'x', cj = plane[1] - 'x';
    auto rows = read_csv(a.pos.at(1));
    const std::string props = a.get("properties", "solid");
    if (props != "solid" && props != "outline") throw std::runtime_error("--properties must be solid or outline");
    const bool outline = props == "outline";
    set_bspline_integration(a.get("bspline", "parametric") != "fan");
    StepModel sm = read_step(a.pos.at(0), a.get("root"), !outline);
    auto groups = group_by_level(sm, level);
    struct Out { std::string name, body, group, src; double m, x, y, k2; bool rem; };
    std::vector<Out> out;
    double m_rem = 0;
    for (const auto& r : rows) {
        Out o;
        o.name = r.at("name");
        o.body = r.at("body");
        o.group = r.count("group") && !r.at("group").empty() ? r.at("group") : o.body;
        o.m = std::stod(r.at("mass_kg"));
        o.rem = !(r.count("remainder") && (r.at("remainder") == "no" || r.at("remainder") == "0"));
        const std::string pre = r.count("prefix") ? r.at("prefix") : std::string();
        if (o.body == "-") {  // takes part in the remainder only (e.g. a cylinder whose mass is set in the model file)
            o.x = o.y = o.k2 = 0;
            o.src = "remainder only";
        } else if (!pre.empty()) {
            auto u = union_by_prefix(groups, pre);
            char buf[200];
            if (outline) {  // centre of the bounding outline, moment of inertia of a uniform rectangular prism
                const double da = (u.hi[ci] - u.lo[ci]) / 1000.0, db = (u.hi[cj] - u.lo[cj]) / 1000.0;
                o.x = 0.5 * (u.lo[ci] + u.hi[ci]);
                o.y = 0.5 * (u.lo[cj] + u.hi[cj]);
                o.k2 = (da * da + db * db) / 12.0;
                std::snprintf(buf, sizeof buf, "STEP %s: outline %.0f x %.0f mm, %zu occurrences", pre.c_str(), da * 1000, db * 1000, u.n_occurrences);
            } else {
                if (u.mass.volume <= 0) throw std::runtime_error("no solid volume for prefix '" + pre + "'");
                o.x = u.mass.centroid[ci];
                o.y = u.mass.centroid[cj];
                o.k2 = u.mass.k2(ci, cj);
                std::snprintf(buf, sizeof buf, "STEP %s: V=%.0f mm3 (%.2f kg at 7850 kg/m3), %zu occurrences, %zu faces triangulated",
                              pre.c_str(), u.mass.volume, u.mass.volume * 7.85e-6, u.n_occurrences, u.faces_triangulated);
            }
            o.src = buf;
            // a bought-in component modelled as a shell or a simplified body: given centre of gravity / inertia
            if (r.count("x_mm") && !r.at("x_mm").empty() && r.count("y_mm") && !r.at("y_mm").empty()) {
                o.x = std::stod(r.at("x_mm"));
                o.y = std::stod(r.at("y_mm"));
                o.src += ", centre of gravity given";
            }
            if (r.count("J_kgm2") && !r.at("J_kgm2").empty()) {
                o.k2 = o.m > 0 ? std::stod(r.at("J_kgm2")) / o.m : 0.0;
                o.src += ", J given";
            }
        } else {
            o.x = std::stod(r.at("x_mm"));
            o.y = std::stod(r.at("y_mm"));
            const double J = r.count("J_kgm2") && !r.at("J_kgm2").empty() ? std::stod(r.at("J_kgm2")) : 0.0;
            o.k2 = o.m > 0 ? J / o.m : 0.0;
            o.src = "given";
        }
        if (o.rem) m_rem += o.m;
        out.push_back(o);
    }
    const double remainder = to_d(a.get("remainder", "0"));
    const double k = m_rem > 0 ? 1.0 + remainder / m_rem : 1.0;
    std::fprintf(stderr, "%zu items, %.3f kg (%.3f kg take part in the remainder), remainder %.3f kg -> scale factor %.10f\n", out.size(),
                 [&] { double s = 0; for (auto& o : out) s += o.m; return s; }(), m_rem, remainder, k);
    std::ofstream f;
    std::ostream& os = outstream(a, f);
    os.precision(10);
    os << "name,body,group,mass_kg,x_mm,y_mm,J_kgm2,source\n";
    for (const auto& o : out) {
        const double m = o.rem ? o.m * k : o.m;
        if (o.body == "-") {
            std::fprintf(stderr, "  %s: scaled mass %.10g kg (set it in the model file)\n", o.name.c_str(), m);
            continue;
        }
        os << '"' << o.name << "\"," << o.body << ',' << o.group << ',' << m << ',' << o.x << ',' << o.y << ',' << m * o.k2 << ",\""
           << o.src << "\"\n";
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2 || !std::strcmp(argv[1], "-h") || !std::strcmp(argv[1], "--help")) {
        std::fputs(kUsage, argc < 2 ? stderr : stdout);
        return argc < 2 ? 1 : 0;
    }
    if (!std::strcmp(argv[1], "--version")) { std::puts(CYLLOAD_VERSION); return 0; }
    try {
        const std::string cmd = argv[1];
        Args a = parse_args(argc, argv, 2);
        if (a.pos.empty()) throw std::runtime_error("missing arguments; see cylload --help");
        if (cmd == "info") return cmd_info(a);
        if (cmd == "pose") return cmd_pose(a);
        if (cmd == "sweep") return cmd_sweep(a);
        if (cmd == "ramp") return cmd_ramp(a);
        if (cmd == "mc") return cmd_mc(a);
        if (cmd == "step") return cmd_step(a);
        if (cmd == "masstable") return cmd_masstable(a);
        if (cmd == "simulate") return cmd_simulate(a);
        if (cmd == "map") return cmd_map(a);
        if (cmd == "plot") return cmd_plot(a);
        if (cmd == "draw") return cmd_draw(a);
        throw std::runtime_error("unknown command '" + cmd + "'; see cylload --help");
    } catch (const std::exception& e) {
        std::fprintf(stderr, "cylload: error: %s\n", e.what());
        return 2;
    }
}
