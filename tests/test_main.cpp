// SPDX-License-Identifier: MIT
// Unit, verification and regression tests for cylload (no external framework).
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <sstream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "cylload/analysis.hpp"
#include "cylload/montecarlo.hpp"
#include "cylload/step.hpp"
#include "cylload/svg.hpp"

using namespace cylload;

namespace {
int g_fail = 0, g_checks = 0;
const std::string SRC = CYLLOAD_SOURCE_DIR;
constexpr double kPi = 3.14159265358979323846;

void check(bool ok, const std::string& what) {
    ++g_checks;
    if (!ok) { ++g_fail; std::printf("  FAIL: %s\n", what.c_str()); }
}
void near(double a, double b, double tol, const std::string& what) {
    const bool ok = std::fabs(a - b) <= tol;
    ++g_checks;
    if (!ok) { ++g_fail; std::printf("  FAIL: %s: got %.10g expected %.10g (tol %.3g)\n", what.c_str(), a, b, tol); }
}
void rel(double a, double b, double rtol, const std::string& what) { near(a, b, rtol * std::fabs(b), what); }

void run(const char* name, const std::function<void()>& f) {
    const int before = g_fail, checks_before = g_checks;
    try {
        f();
    } catch (const std::exception& e) {
        ++g_fail;
        std::printf("  FAIL: exception: %s\n", e.what());
    }
    std::printf("[%s] %s (%d checks)\n", g_fail == before ? " ok " : "FAIL", name, g_checks - checks_before);
}

// ------------------------------------------------------------------------------------------
// 1. Synthetic boom: closed form from the law of cosines, independent of the solver.
void test_boom_closed_form() {
    Model m = Model::from_file(SRC + "/examples/boom/boom.mech");
    Kinematics kin(m);
    const double Kx = 300, Ky = -400, Jx = 900, Jy = 150;
    const double OK = std::hypot(Kx, Ky), OJ = std::hypot(Jx, Jy);
    const double a0 = std::atan2(Jy, Jx) - std::atan2(Ky, Kx);
    struct It { double m, x, y, j; };
    const std::vector<It> items{{250, 1500, 200, 40}, {120, 2900, 300, 5}, {10, Jx, Jy, 0}};  // cylinder half at J
    for (double L = 700; L <= 1200; L += 25) {
        Pose s = kin.solve(kin.reference(), {L});
        const double al = std::acos((OK * OK + OJ * OJ - L * L) / (2 * OK * OJ));
        const double th = al - a0, dth = L / (OK * OJ * std::sin(al));  // rad, rad/mm
        double M = 0, mr = 0;
        for (const auto& it : items) {
            const double xr = std::cos(th) * it.x - std::sin(th) * it.y;
            M += it.m * xr;
            mr += it.m * (it.x * it.x + it.y * it.y) * dth * dth + it.j * (dth * 1000) * (dth * 1000);
        }
        const double Fc = 9.81 * M * dth;
        auto st = statics(kin, s);
        rel(kin.angle(s, 1), th, 1e-12, "boom angle");
        rel(st.force_vw[0], Fc, 1e-12, "boom virtual work vs closed form");
        rel(st.force_eq[0], Fc, 1e-12, "boom link equilibrium vs closed form");
        rel(forces_finite_difference(kin, s)[0], Fc, 1e-8, "boom finite difference vs closed form");
        rel(reduced_mass(kin, s, 0), mr, 1e-10, "boom reduced mass vs closed form");
        // pin reaction = total weight minus vertical cylinder force component (vertical balance)
        Vec2 u = kin.actuator_dir(s, 0);
        rel(st.reactions[0].fy + Fc * u[1], 9.81 * 380.0, 1e-10, "boom vertical balance");
    }
}

// 1b. The same boom with a constant external force (e.g. a hitch or digging load) at a boom point:
// its contribution to the cylinder force is -f . dr/dL from the closed form.
void test_external_force() {
    std::ifstream fin(SRC + "/examples/boom/boom.mech");
    std::stringstream ss;
    ss << fin.rdbuf() << "\nforce tip boom 2900 300 1500 -4000\n";
    Model m = Model::from_string(ss.str(), SRC + "/examples/boom");
    Model m0 = Model::from_file(SRC + "/examples/boom/boom.mech");
    Kinematics kin(m), kin0(m0);
    const double Kx = 300, Ky = -400, Jx = 900, Jy = 150;
    const double OK = std::hypot(Kx, Ky), OJ = std::hypot(Jx, Jy);
    const double a0 = std::atan2(Jy, Jx) - std::atan2(Ky, Kx);
    for (double L = 700; L <= 1200; L += 50) {
        Pose s = kin.solve(kin.reference(), {L}), s0 = kin0.solve(kin0.reference(), {L});
        const double al = std::acos((OK * OK + OJ * OJ - L * L) / (2 * OK * OJ));
        const double th = al - a0, dth = L / (OK * OJ * std::sin(al));
        const double drx = -std::sin(th) * 2900 - std::cos(th) * 300, dry = std::cos(th) * 2900 - std::sin(th) * 300;
        const double Fext = -(1500 * drx - 4000 * dry) * dth;
        const double F0 = statics(kin0, s0).force_vw[0];
        auto st = statics(kin, s);
        rel(st.force_vw[0] - F0, Fext, 1e-10, "external force: virtual work vs closed form");
        rel(st.force_eq[0] - F0, Fext, 1e-10, "external force: link equilibrium vs closed form");
        rel(forces_finite_difference(kin, s)[0], st.force_vw[0], 1e-8, "external force: energy differences");
        rel(st.reactions[0].fy + st.force_vw[0] * kin.actuator_dir(s, 0)[1], 9.81 * 380.0 + 4000.0, 1e-10, "external force: vertical balance");
    }
}

// 2. Block on an inclined slide pushed by a cylinder along the slide: F = m g sin(beta).
void test_incline() {
    const double beta = 30.0 * kPi / 180.0;
    const double ux = std::cos(beta), uy = std::sin(beta);
    char txt[1024];
    std::snprintf(txt, sizeof txt,
                  "body block\nprismatic slide ground block 0 0 %.15g %.15g\n"
                  "actuator cyl ground %.15g %.15g block 0 0 bore=50 rod=25 min=100 max=900\n"
                  "mass load block 0 0 200\n",
                  ux, uy, -500 * ux, -500 * uy);
    Model m = Model::from_string(txt);
    Kinematics kin(m);
    for (double L : {300.0, 500.0, 800.0}) {
        Pose s = kin.solve(kin.reference(), {L});
        auto st = statics(kin, s);
        rel(st.force_vw[0], 200 * 9.81 * std::sin(beta), 1e-12, "incline virtual work");
        rel(st.force_eq[0], 200 * 9.81 * std::sin(beta), 1e-12, "incline equilibrium");
        rel(std::hypot(st.reactions[0].fx, st.reactions[0].fy), 200 * 9.81 * std::cos(beta), 1e-12, "incline normal force");
        rel(reduced_mass(kin, s, 0), 200.0, 1e-12, "incline reduced mass");
    }
}

// 3. Formulas of the hydraulic part.
void test_hydraulics() {
    Cylinder c;
    c.bore = 80; c.rod = 40; c.lmin = 758; c.lmax = 1198;
    rel(c.area_piston(), kPi * 1600, 1e-14, "piston area");
    rel(c.area_annulus(), kPi * 1200, 1e-14, "annulus area");
    rel(pressure_bar(c, 50000), 50000 / (kPi * 1600) * 10, 1e-14, "push pressure");
    rel(pressure_bar(c, -20000), 20000 / (kPi * 1200) * 10, 1e-14, "pull pressure");
    Hydraulics h;
    const double Ap = kPi * 0.04 * 0.04, Ar = Ap - kPi * 0.02 * 0.02;
    rel(hydraulic_stiffness(c, h, 900), 1e9 * (Ap * Ap / (3e-4 + Ap * 0.142) + Ar * Ar / (3e-4 + Ar * 0.298)), 1e-12, "stiffness");
    near(ramp_factor(0), 1.0, 1e-15, "kappa(0)");
    rel(ramp_factor(1e-6), 1.0, 1e-9, "kappa -> 1");
    rel(ramp_factor(kPi), 2 / kPi, 1e-12, "kappa(pi)");
    rel(ramp_factor(10), 0.2, 1e-14, "kappa(10)");
    double prev = 2;
    bool mono = true;  // one check for the whole grid of 3000 points
    for (double t = 0; t < 30; t += 0.01) { mono &= ramp_factor(t) <= prev + 1e-15; prev = ramp_factor(t); }
    check(mono, "kappa non-increasing on a grid of 3000 points");
    rel(euler_load(c, 1000), kPi * kPi * 210000 * kPi * std::pow(40, 4) / 64 / 1e6, 1e-12, "Euler load");
    // version 1.2: back-pressure in the unloaded chamber and separate dead volumes
    c.p_back = 20;
    rel(pressure_bar(c, 50000), 50000 / (kPi * 1600) * 10 + 20 * 1200.0 / 1600.0, 1e-14, "push pressure with back-pressure");
    rel(pressure_bar(c, -20000), 20000 / (kPi * 1200) * 10 + 20 * 1600.0 / 1200.0, 1e-14, "pull pressure with back-pressure");
    Hydraulics h2;
    h2.dead_volume = 0.5;
    h2.dead_volume_rod = 1.5;
    rel(hydraulic_stiffness(c, h2, 900), 1e9 * (Ap * Ap / (5e-4 + Ap * 0.142) + Ar * Ar / (1.5e-3 + Ar * 0.298)), 1e-12,
        "stiffness with different dead volumes");
    Model hm = Model::from_string(
        "body b\nrevolute O ground b 0 0\nactuator c ground 0 -500 b 300 0 bore=80 rod=40 min=500 max=900\n"
        "mass m b 600 0 100 0\nhydraulics beta=1.2e9 dead_volume=0.4 dead_volume_rod=0.9 back_pressure=15\n");
    check(hm.hydraulics.dead_volume == 0.4 && hm.hydraulics.dead_volume_rod == 0.9, "dead volumes parsed (dead_volume sets both, then rod side)");
    check(hm.actuators[0].cyl.p_back == 15, "back_pressure copied to the cylinders");
}

// 4. Horty tower, simplified variant (outline centres, prism inertia): regression against the
//    independent Python implementation used for applsci-4655781 (Table A2, first column).
void test_horty_outline() {
    Model m = Model::from_file(SRC + "/examples/horty/tower_outline.mech");
    Kinematics kin(m);
    rel(m.total_mass(), 1150.0 - 39.4 + 56.0, 1e-8, "tower mass incl. post");
    auto L = kin.reference().lengths;
    L[m.actuator_index("S8")] = 973.0;
    L[m.actuator_index("S10")] = 1183.0;
    Pose s = kin.solve(kin.reference(), L);
    auto st = statics(kin, s);
    rel(st.force_vw[1], -13867.778194, 2e-8, "F10, mast vertical");
    rel(st.force_vw[2], 7120.218206, 2e-8, "F19");
    for (std::size_t k = 0; k < 4; ++k) rel(st.force_eq[k], st.force_vw[k], 1e-10, "Horty equilibrium vs virtual work");

    SweepOptions o8;
    o8.actuator = 0; o8.fixed = L; o8.from = 698; o8.to = 1048; o8.n = 141;
    auto p8 = sweep(kin, o8);
    auto s8 = summarize(kin, o8, p8);
    rel(s8.F_min, 9505.819157, 2e-8, "F8 min");
    rel(s8.F_max, 48473.955202, 2e-8, "F8 max");
    near(s8.L_at_F_max, 698.0, 1e-9, "F8 max at retracted");
    rel(s8.p_max, 61.718957, 2e-8, "p8 max");
    rel(s8.work_J, s8.dV_J, 1e-8, "energy balance, lift (Simpson)");
    rel(s8.dV_J, 9961.508577, 1e-8, "potential change, lift");
    check(s8.dev_fd < 1e-8 && s8.dev_eq < 1e-12, "lift: three force routes agree");
    rel(s8.m_red_min, 15925.25, 0.01, "lift reduced mass (cylinder pins included)");

    SweepOptions o10 = o8;
    o10.actuator = 1; o10.from = 758; o10.to = 1198; o10.n = 177;
    auto p10 = sweep(kin, o10);
    auto s10 = summarize(kin, o10, p10);
    rel(s10.F_max, 50838.925584, 2e-8, "F10 max");
    near(s10.L_at_F_max, 758.0, 1e-9, "F10 max at retracted");
    near(p10.front().angle_deg[m.body_index("W") - 1], -90.6184469, 1e-6, "tower tilt at L10 min");
    near(p10.back().angle_deg[m.body_index("W") - 1], 4.0154317, 1e-6, "tower tilt at L10 max");
    rel(s10.F_min, -17958.141964, 2e-8, "F10 min");
    rel(s10.p_max, 101.140829, 2e-8, "p10 max");
    rel(s10.work_J, s10.dV_J, 1e-7, "energy balance, tilt (Simpson)");
    rel(s10.dV_J, 6686.314399, 1e-8, "potential change, tilt");
    check(s10.dev_fd < 1e-8 && s10.dev_eq < 1e-12, "tilt: three force routes agree");
    rel(s10.f_min, 3.7088, 0.005, "tilt natural frequency min");
    rel(s10.f_max, 4.4650, 0.005, "tilt natural frequency max");
    // valve ramp: the published grid search gave 0.445 s (step 1 ms) at Q = 11.25 L/min
    const double tr = min_ramp_time(m.actuators[1].cyl, p10, 11.25, +1, 140.0);
    near(tr, 0.445, 0.002, "minimum ramp time, tilt, extension");
    check(peak_pressure(m.actuators[1].cyl, p10, 11.25, +1, tr) <= 140.0 + 1e-9, "ramp result respects the limit");

    // hammer position along the mast changes the tilt load but not the lift load
    auto Lh = L;
    Lh[m.actuator_index("winch")] += 1500;
    Lh[m.actuator_index("S10")] = 758;
    Pose sh = kin.solve(s, Lh);
    auto Lb = Lh;
    Lb[m.actuator_index("winch")] = L[m.actuator_index("winch")];
    Pose sb = kin.solve(s, Lb);
    check(forces_virtual_work(kin, sh)[1] > forces_virtual_work(kin, sb)[1], "raised hammer increases the tilt load");
    rel(forces_virtual_work(kin, sh)[0], forces_virtual_work(kin, sb)[0], 1e-12, "lift load independent of hammer position");
}

void test_post_lift_outline() {
    Model m = Model::from_file(SRC + "/examples/horty/post_lift_outline.mech");
    Kinematics kin(m);
    SweepOptions o;
    o.actuator = 0; o.fixed = {630}; o.from = 630; o.to = 1030; o.n = 121;
    auto p = sweep(kin, o);
    auto s = summarize(kin, o, p);
    rel(s.F_max, 12585.543808, 2e-8, "post lift F max");
    rel(s.p_max, 64.097648, 2e-8, "post lift p max");
    rel(s.F_min, -4230.021817, 2e-8, "post lift F min");
    rel(s.f_min, 2.2488, 0.005, "post lift natural frequency min");
}

// Stone picker (examples/usarya): kinematics against the closed-form triangle model of the machine
// (Mathcad) and pressures against the published values (USARYA 2, Agriculture).
void test_usarya() {
    const std::string dir = SRC + "/examples/usarya/";
    {   // folding of the working module: Mathcad gives angle Sp-Lp1-PSRRp = 2.4037490335 rad at S = 1657 mm
        Model m = Model::from_file(dir + "rotor.mech");
        Kinematics kin(m);
        const Pose s = kin.solve(kin.reference(), {1657.0});
        const Vec2 lp{890.0, 0.0}, sp{-325.26912, 325.26912};
        const Vec2 p = kin.point(s, 1, {964.49487, 488.08454});
        const double v0 = p[0] - lp[0], v1 = p[1] - lp[1], w0 = sp[0] - lp[0], w1 = sp[1] - lp[1];
        near(std::fabs(std::atan2(v0 * w1 - v1 * w0, v0 * w0 + v1 * w1)), 2.4037490335, 5e-7, "stone picker folding angle vs Mathcad");
        // closed form: chi(S) from |PSRRp(chi) - Sp| = S, chi measured as in the Mathcad model
        auto chi_of = [&](double S) {
            double a = -0.3, b = 1.9;
            auto f = [&](double c) {
                const double e1x = std::cos(c), e1y = std::sin(c), e2x = std::sin(c), e2y = -std::cos(c);
                const double x = 890.0 + 440 * e1x - 224 * e2x, y = 440 * e1y - 224 * e2y;
                return std::hypot(x - sp[0], y - sp[1]) - S;
            };
            for (int i = 0; i < 200; ++i) { const double c = 0.5 * (a + b); (f(a) * f(c) <= 0 ? b : a) = c; }
            return 0.5 * (a + b);
        };
        double dmax = 0;
        for (double S = 1023; S <= 1657; S += 31.7) {
            const Pose q = kin.solve(kin.reference(), {S});
            dmax = std::max(dmax, std::fabs(0.9484443460058755 + kin.angle(q, 1) - chi_of(S)));
        }
        near(dmax, 0.0, 1e-7, "stone picker folding: rotation vs closed form over the stroke");
        SweepOptions o; o.actuator = 0; o.fixed = {1013}; o.from = 1013; o.to = 1583; o.n = 58;
        auto pts = sweep(kin, o);
        auto sm = summarize(kin, o, pts);
        near(sm.p_max, 65.03, 0.01, "stone picker folding: maximum pressure (published 65.0 bar)");
    }
    {   // drawbar: frame tilt 3.1325305 deg at ED = 737 mm (closed form), 10.1929247 deg at the reference 935 mm
        Model m = Model::from_file(dir + "drawbar.mech");
        Kinematics kin(m);
        const Pose s = kin.solve(kin.reference(), {737.0});
        near(10.192924729 + kin.angle(s, 1) * 180.0 / M_PI, 3.1325305278, 1e-6, "stone picker drawbar: frame tilt vs closed form");
    }
    struct Case { const char* file; double from, to, pmax; };
    for (const Case& c : {Case{"cases/lift_m6000_b0.mech", 1150.388, 1890.388, 145.61}, Case{"cases/lift_m6000_b3.1325.mech", 1150.388, 1890.388, 150.15},
                          Case{"cases/lift_m6000_b10.1929.mech", 1150.388, 1890.388, 159.22}, Case{"cases/lift_m0_b0.mech", 1150.388, 1890.388, 23.49},
                          Case{"cases/tip_m6000_b0.mech", 988.435, 1588.435, 143.76}, Case{"cases/tip_m6000_b3.1325.mech", 988.435, 1588.435, 146.69},
                          Case{"cases/tip_m6000_b10.1929.mech", 988.435, 1588.435, 151.70}}) {
        Model m = Model::from_file(dir + c.file);
        Kinematics kin(m);
        SweepOptions o; o.actuator = 0; o.fixed = {c.from}; o.from = c.from; o.to = c.to; o.n = 75;
        auto pts = sweep(kin, o);
        auto sm = summarize(kin, o, pts);
        near(sm.p_max, c.pmax, 0.011, std::string("stone picker ") + c.file + ": maximum pressure vs published");
        check(sm.dev_fd < 1e-8 && sm.dev_eq < 1e-8, std::string("stone picker ") + c.file + ": three routes agree");
    }
}

// 5. Monte Carlo helpers.
void test_monte_carlo() {
    std::mt19937_64 ref(5489);
    Uniform u(5489);
    bool same = true;
    for (int i = 0; i < 20000; ++i) same &= (u() == static_cast<double>(ref() >> 11) * (1.0 / 9007199254740992.0));
    check(same, "Uniform reproduces std::mt19937_64");
    std::vector<double> x, y, z;
    for (int i = 0; i < 100; ++i) { x.push_back(i); y.push_back(std::exp(0.05 * i)); z.push_back(-i * i); }
    near(spearman(x, y), 1.0, 1e-14, "Spearman monotone increasing");
    near(spearman(x, z), -1.0, 1e-14, "Spearman monotone decreasing");
    near(spearman({1, 2, 2, 3}, {1, 2, 2, 3}), 1.0, 1e-14, "Spearman with ties");
    near(percentile({1, 2, 3, 4, 5}, 50), 3, 1e-14, "median");
    near(percentile({1, 2, 3, 4}, 5), 1.15, 1e-14, "5th percentile, linear");

    Model m = Model::from_file(SRC + "/examples/boom/boom.mech");
    McSpec sp = McSpec::from_string(
        "seed 7\nsamples 400\nscale k_arm item:arm 0.8 1.2\nmass m_tool item:tool 50 150\n"
        "length L cyl 750 1150\noutput F cyl\n");
    auto r = run_monte_carlo(m, sp);
    check(r.inputs.size() == 400 && r.force_stats.size() == 1, "MC sizes");
    check(r.spearman[0][1] > 0.5 && r.spearman[0][0] > 0.3, "MC: tool mass dominates, arm mass matters");
    check(r.force_stats[0].p5 < r.force_stats[0].p50 && r.force_stats[0].p50 < r.force_stats[0].p95, "MC percentiles ordered");
    // the same seed gives the same samples
    auto r2 = run_monte_carlo(m, sp);
    check(r2.force[123][0] == r.force[123][0], "MC reproducible");
}

// 6. STEP reader on a hand-written two-level assembly.
void test_step() {
    StepModel s = read_step(SRC + "/tests/data/mini_assembly.step");
    check(s.root == "TOP", "STEP root found");
    check(s.occurrences.size() == 2, "STEP: two occurrences with geometry");
    auto g = group_by_level(s, 2);
    for (const auto& o : s.occurrences) {
        if (o.path.back() == "Bracket") {
            near(o.lo[0], 980, 1e-9, "bracket xmin"); near(o.hi[0], 1000, 1e-9, "bracket xmax");
            near(o.lo[1], 500, 1e-9, "bracket ymin"); near(o.hi[1], 600, 1e-9, "bracket ymax");
            near(o.hi[2], 10, 1e-9, "bracket zmax");
            check(o.n_points == 4, "bracket points (placements excluded)");
        } else {
            check(o.path.back() == "Pin \xC5\xBC\xC3\xB3\xC5\x82ty", "STEP \\X2\\ name decoding");
            check(o.path.size() == 3, "pin path depth");
            near(o.lo[0], 985, 1e-9, "pin xmin"); near(o.hi[0], 995, 1e-9, "pin xmax");
            near(o.lo[1], 545, 1e-9, "pin ymin"); near(o.hi[1], 555, 1e-9, "pin ymax");
            near(o.lo[2], 10, 1e-9, "pin zmin"); near(o.hi[2], 50, 1e-9, "pin zmax");
        }
    }
    auto u = union_by_prefix(group_by_level(s, 1), "Bracket");
    near(u.lo[2], 0, 1e-9, "union zmin"); near(u.hi[2], 50, 1e-9, "union zmax (bracket + pin)");
    check(decode_step_string("it''s \\X\\E9") == "it's \xC3\xA9", "STEP \\X\\ and quote decoding");
}


// 5. Horty, mass properties from the solid model (adopted in applsci-4655781): published values
//    (Tables 4, 5 and A2; curve data of Figures 4a, 4b and 5a). Published values are rounded,
//    and the reduced mass here also contains the cylinder masses lumped at their pins (0.3-1.1 %).
void test_horty_solid() {
    Model m = Model::from_file(SRC + "/examples/horty/tower.mech");
    Kinematics kin(m);
    auto L = kin.reference().lengths;
    L[m.actuator_index("S8")] = 973.0;
    L[m.actuator_index("S10")] = 1183.0;
    SweepOptions o8;
    o8.actuator = 0; o8.fixed = L; o8.from = 698; o8.to = 1048; o8.n = 141;
    auto p8 = sweep(kin, o8);
    rel(p8.front().p, 61.84077, 1e-4, "p8 tower lowered (Fig. 4a data)");
    near(p8.back().p, 12.2, 0.05, "p8 tower raised (Table 4)");
    SweepOptions o10 = o8;
    o10.actuator = 1; o10.from = 758; o10.to = 1198; o10.n = 177;
    auto p10 = sweep(kin, o10);
    rel(p10.front().p, 108.41555, 1e-4, "p10 mast horizontal (Fig. 4b data)");
    near(p10.back().F / 1000, -20.0, 0.05, "F10 end of stroke (Table 4)");
    near(p10.back().p, 53.1, 0.05, "p10 end of stroke (Table 4)");
    const int iw = m.body_index("W") - 1;
    double th0 = 0;
    for (std::size_t i = 0; i + 1 < p10.size(); ++i)
        if (p10[i].F * p10[i + 1].F < 0) th0 = p10[i].angle_deg[iw];
    near(th0, -18.2, 0.6, "load reversal angle (Table 4)");
    Pose s = kin.solve(kin.reference(), L);
    auto st = statics(kin, s);
    near(st.force_vw[1] / 1000, -15.6, 0.05, "F10 mast vertical (Table 4)");
    near(pressure_bar(m.actuators[1].cyl, st.force_vw[1]), 41.4, 0.05, "p10 mast vertical (Table 4)");
    near(st.force_vw[2] / 1000, 7.1, 0.05, "F19 (Table 4)");
    const auto& c10 = m.actuators[1].cyl;
    rel(min_ramp_time(c10, p10, 56.1, -1, 140.0), 3.80, 0.015, "t_r* tilting, 56.1 L/min (Table 5)");
    rel(min_ramp_time(c10, p10, 56.1, -1, 160.0), 2.33, 0.015, "t_r* tilting, 160 bar (Table 5)");
    near(min_ramp_time(m.actuators[0].cyl, p8, 56.1, -1, 140.0), 0.12, 0.006, "t_r* lifting (Table 5)");
    auto s10 = summarize(kin, o10, p10);
    near(s10.f_min, 3.6, 0.05, "tilt natural frequency min (Table 5)");
    near(s10.f_max, 4.4, 0.08, "tilt natural frequency max (Table 5)");

    Model pm = Model::from_file(SRC + "/examples/horty/post_lift.mech");
    Kinematics pk(pm);
    SweepOptions op;
    op.actuator = 0; op.fixed = {630}; op.from = 630; op.to = 1030; op.n = 121;
    auto pp = sweep(pk, op);
    rel(pp.front().p, 57.95876, 5e-4, "post lift p at L=630 (Fig. 5a data)");
    near(pp.front().F / 1000, 11.4, 0.05, "post lift F (Table 4)");
    rel(min_ramp_time(pm.actuators[0].cyl, pp, 56.1, -1, 140.0), 5.22, 0.01, "t_r* post lifting 56 kg (Table 5)");
    for (auto& mi : pm.masses)
        if (mi.name == "post") { mi.m = 110.0; mi.j = 110.0 * 16.0 / 12.0; }
    Kinematics pk2(pm);
    auto pp2 = sweep(pk2, op);
    rel(min_ramp_time(pm.actuators[0].cyl, pp2, 56.1, -1, 140.0), 7.35, 0.01, "t_r* post lifting 110 kg (Table 5)");
}

// 6. STEP mass properties: solids of revolution with all analytic face types, a cone apex,
//    sphere poles and a lemon torus (negative major radius), placed with rotations.
void test_step_revolved() {
    StepModel s = read_step(SRC + "/tests/data/revolved_solids.step", "", true);
    std::ifstream f(SRC + "/tests/data/revolved_reference.csv");
    check(static_cast<bool>(f), "reference table present");
    std::string line;
    std::getline(f, line);
    int n = 0;
    while (std::getline(f, line)) {
        std::stringstream ss(line);
        std::string name, x;
        std::getline(ss, name, ',');
        std::vector<double> v;
        while (std::getline(ss, x, ',')) v.push_back(std::stod(x));
        for (const auto& o : s.occurrences) {
            if (o.path.back() != name) continue;
            ++n;
            const double V = v[0];
            const double scaleL = std::cbrt(V);
            rel(o.volume, V, 1e-8, name + " volume");
            for (int i = 0; i < 3; ++i) near(o.centroid[i], v[1 + i], 1e-7 * scaleL, name + " centroid");
            const double tr = v[4] + v[7] + v[9];
            const int idx[6][2] = {{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}};
            for (int k = 0; k < 6; ++k) near(o.cov[3 * idx[k][0] + idx[k][1]], v[4 + k], 1e-8 * tr, name + " second moment");
            check(o.faces_triangulated == 0, name + ": all faces integrated analytically");
            check(o.volume_check < 1e-9 && o.moment_check < 1e-9, name + ": field consistency");
        }
    }
    check(n == 5, "all five test solids found");
}

// 6b. Version 1.2: B-spline faces integrated in their parameter domain (doubly curved bicubic face with
//     interior knots, convex and concave rational fillets); version 1.1 replaced them by triangle fans.
struct RefRow { std::string name; std::vector<double> v; };
std::vector<RefRow> read_reference(const std::string& path) {
    std::ifstream f(path);
    std::vector<RefRow> out;
    std::string line;
    std::getline(f, line);
    while (std::getline(f, line)) {
        std::stringstream ss(line);
        RefRow r;
        std::string x;
        std::getline(ss, r.name, ',');
        while (std::getline(ss, x, ',')) r.v.push_back(std::stod(x));
        out.push_back(r);
    }
    return out;
}
void test_step_bspline() {
    const auto ref = read_reference(SRC + "/tests/data/bspline_reference.csv");
    check(ref.size() == 3, "reference table present");
    for (bool parametric : {true, false}) {
        set_bspline_integration(parametric);
        StepModel s = read_step(SRC + "/tests/data/bspline_solids.step", "", true);
        int n = 0;
        for (const auto& r : ref)
            for (const auto& o : s.occurrences) {
                if (o.path.back() != r.name) continue;
                ++n;
                const auto& v = r.v;
                if (parametric) {
                    const double L = std::cbrt(v[0]), tr = v[4] + v[7] + v[9];
                    rel(o.volume, v[0], 1e-9, r.name + " volume");
                    for (int i = 0; i < 3; ++i) near(o.centroid[i], v[1 + i], 1e-8 * L, r.name + " centroid");
                    const int idx[6][2] = {{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}};
                    for (int k = 0; k < 6; ++k) near(o.cov[3 * idx[k][0] + idx[k][1]], v[4 + k], 1e-8 * tr, r.name + " second moment");
                    check(o.faces_parametric == 1 && o.faces_triangulated == 0, r.name + ": B-spline face integrated in (u, v)");
                    check(o.volume_check < 1e-9 && o.moment_check < 1e-9, r.name + ": field consistency");
                } else {
                    // the version 1.1 fans miss the curvature of the face interior: errors of about 1-2 %
                    const double e = std::fabs(o.volume - v[0]) / v[0];
                    check(e > 5e-3 && e < 3e-2, r.name + ": triangle-fan volume error of 0.5-3 %");
                    check(o.faces_triangulated == 1, r.name + ": B-spline face triangulated");
                }
            }
        check(n == 3, "all three B-spline test solids found");
    }
    set_bspline_integration(true);
}

// STEP length units: the same solids written in mm, inch, m and without a unit.
void test_step_units() {
    std::ifstream f(SRC + "/tests/data/bspline_solids.step");
    std::stringstream ss;
    ss << f.rdbuf();
    const std::string mm = ss.str(), unit = "#4 = (LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));";
    check(mm.find(unit) != std::string::npos, "unit entity present in the test file");
    set_bspline_integration(true);
    const StepModel ref = read_step(SRC + "/tests/data/bspline_solids.step", "", true);
    check(ref.length_unit == "mm" && ref.length_unit_mm == 1.0, "millimetre file recognised");
    struct Case { std::string name, entity; double f; };
    const std::vector<Case> cases = {
        {"inch", "#4 = (CONVERSION_BASED_UNIT('INCH',#9001)LENGTH_UNIT()NAMED_UNIT(#9002));\n"
                 "#9000 = (LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));\n"
                 "#9001 = LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#9000);\n"
                 "#9002 = DIMENSIONAL_EXPONENTS(1.,0.,0.,0.,0.,0.,0.);", 1.0 / 25.4},
        {"m", "#4 = (LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT($,.METRE.));", 1e-3},
        {"none", "#4 = DIMENSIONAL_EXPONENTS(1.,0.,0.,0.,0.,0.,0.);", 1.0}};
    for (const auto& c : cases) {
        // same file with the coordinates divided by the unit, so the results must be equal in mm
        std::string txt = mm;
        txt.replace(txt.find(unit), unit.size(), c.entity);
        if (c.f != 1.0) {
            std::string out;
            std::size_t p = 0;
            while (true) {
                const auto q = txt.find("CARTESIAN_POINT('',(", p);
                if (q == std::string::npos) { out += txt.substr(p); break; }
                const auto o = q + 20, e = txt.find(')', o);
                out += txt.substr(p, o - p);
                std::stringstream in(txt.substr(o, e - o));
                std::string x;
                bool first = true;
                while (std::getline(in, x, ',')) {
                    char buf[40];
                    std::snprintf(buf, sizeof buf, "%.17g", std::stod(x) * c.f);
                    out += (first ? "" : ",") + std::string(buf);
                    first = false;
                }
                p = e;
            }
            txt = out;
            // radii of circles and lengths of vectors are lengths too
            for (const char* key : {"VECTOR(", "CIRCLE(", "CYLINDRICAL_SURFACE(", "CONICAL_SURFACE(", "SPHERICAL_SURFACE(", "TOROIDAL_SURFACE("}) {
                std::size_t k = 0;
                while ((k = txt.find(key, k)) != std::string::npos) {
                    const auto e = txt.find(");", k);
                    auto c1 = txt.rfind(',', e);
                    const double v = std::stod(txt.substr(c1 + 1, e - c1 - 1));
                    char buf[40];
                    std::snprintf(buf, sizeof buf, "%.17g", v * c.f);
                    txt.replace(c1 + 1, e - c1 - 1, buf);
                    k = e;
                }
            }
        }
        const std::string path = (std::filesystem::temp_directory_path() / ("cylload_units_" + c.name + ".step")).string();
        std::ofstream(path) << txt;
        const StepModel s = read_step(path, "", true);
        if (c.name == "none") {
            check(s.length_unit == "mm (assumed)", "missing unit reported as assumed mm");
            continue;
        }
        check(std::fabs(s.length_unit_mm - 1.0 / c.f) < 1e-12 * s.length_unit_mm, "unit factor " + c.name);
        for (std::size_t i = 0; i < ref.occurrences.size(); ++i) {
            const auto& a = s.occurrences[i];
            const auto& b = ref.occurrences[i];
            rel(a.volume, b.volume, 1e-9, c.name + ": volume converted to mm3");
            for (int k = 0; k < 3; ++k) near(a.centroid[k], b.centroid[k], 1e-7, c.name + ": centroid converted to mm");
            rel(a.cov[0] + a.cov[4] + a.cov[8], b.cov[0] + b.cov[4] + b.cov[8], 1e-9, c.name + ": second moments converted");
            rel(a.area, b.area, 1e-9, c.name + ": area converted");
            for (int k = 0; k < 3; ++k) near(a.hi[k], b.hi[k], 1e-7, c.name + ": bounding box converted");
        }
    }
}

// Open example: synthetic backhoe attachment, exact assembly properties and three force routes.
void test_excavator() {
    std::ifstream f(SRC + "/examples/excavator/excavator_reference.csv");
    std::string line;
    std::vector<std::pair<std::string, std::vector<double>>> ref;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#' || line.rfind("assembly", 0) == 0) continue;
        std::stringstream ss(line);
        std::string name, x;
        std::getline(ss, name, ',');
        std::vector<double> v;
        while (std::getline(ss, x, ',')) v.push_back(std::stod(x));
        ref.push_back({name, v});
    }
    check(ref.size() == 13, "excavator reference table");
    set_bspline_integration(true);
    const StepModel sm = read_step(SRC + "/examples/excavator/excavator.step", "", true);
    const auto groups = group_by_level(sm, 1);
    std::size_t bs = 0;
    for (const auto& [name, v] : ref) {
        const auto g = union_by_prefix(groups, name);
        const double V = v[1], L = std::cbrt(V), tr = v[5] + v[8] + v[10];
        check(g.n_occurrences == static_cast<std::size_t>(v[0]), name + ": parts");
        rel(g.mass.volume, V, 1e-7, name + ": volume");
        for (int k = 0; k < 3; ++k) near(g.mass.centroid[k], v[2 + k], 1e-7 * L, name + ": centroid");
        const int idx[6][2] = {{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}};
        for (int k = 0; k < 6; ++k) near(g.mass.cov[3 * idx[k][0] + idx[k][1]], v[5 + k], 1e-7 * tr, name + ": second moment");
        check(g.faces_triangulated == 0, name + ": no triangulated face");
        bs += g.faces_parametric;
    }
    check(bs == 53, "53 B-spline faces integrated in their parameter domain");
    Model m = Model::from_file(SRC + "/examples/excavator/excavator.mech");
    Kinematics kin(m);
    for (int k = 0; k < 3; ++k) {
        const auto& a = m.actuators[k];
        for (double t : {0.0, 0.25, 0.5, 0.75, 1.0}) {
            Pose s = kin.reference();
            auto Lv = s.lengths;
            Lv[k] = a.cyl.lmin + t * (a.cyl.lmax - a.cyl.lmin);
            s = kin.solve(s, Lv);
            const auto fv = forces_virtual_work(kin, s);
            const auto fe = statics(kin, s).force_eq;
            for (int j = 0; j < 3; ++j) rel(fe[j], fv[j], 1e-9, "excavator: equilibrium vs virtual work");
        }
    }
}

// Bootstrap confidence intervals of the Monte Carlo estimators.
void test_bootstrap() {
    Model m = Model::from_file(SRC + "/examples/boom/boom.mech");
    McSpec sp = McSpec::from_string(
        "seed 11\nsamples 2000\nbootstrap 300\nscale k_arm item:arm 0.8 1.2\nmass m_tool item:tool 50 150\n"
        "scale k_dummy group:none 0.9 1.1\noutput F cyl\n");
    auto r = run_monte_carlo(m, sp);
    const auto& st = r.pressure_stats[0];
    check(st.p5_lo <= st.p5 && st.p5 <= st.p5_hi && st.p95_lo <= st.p95 && st.p95 <= st.p95_hi, "percentiles inside their intervals");
    check(st.p95_hi - st.p95_lo < 0.05 * (st.max - st.min), "narrow interval for N = 2000");
    for (int j = 0; j < 2; ++j)  // active inputs (the interval of an inactive one is bias-corrected towards zero)
        check(r.first_order_lo[0][j] <= r.first_order[0][j] && r.first_order[0][j] <= r.first_order_hi[0][j], "S_i inside its interval");
    // an input that does not act: S_i is only the estimator bias, about (bins - 1) / N
    check(r.first_order[0][2] < 0.03 && r.first_order_lo[0][2] < 1e-3, "inactive input: S_i near its bias, interval reaches zero");
    check(r.first_order_lo[0][1] > 0.3, "dominant input clearly non-zero");
    sp.bootstrap = 0;
    auto r0 = run_monte_carlo(m, sp);
    check(std::isnan(r0.first_order_lo[0][0]), "no intervals without bootstrap");
}

// 8. Version 1.1: time-domain simulation, maps, sensitivity indices, figures.
void test_simulation() {
    Model m = Model::from_file(SRC + "/examples/horty/tower.mech");
    Kinematics kin(m);
    SimOptions o;
    o.actuator = 1;
    o.fixed = kin.reference().lengths;
    o.fixed[0] = 973;
    o.fixed[1] = 978;  // mid-stroke
    o.flow_lpm = 11.25;
    for (double tr : {0.0001, 0.1, 0.3, 0.5}) {
        o.t_ramp = tr;
        auto r = simulate(kin, o);
        rel(r.dF_max, r.dF_formula, 0.01, "simulation vs ramp formula (paper: 0.2-1.2 %)");
        check(!r.hit_end, "start-up stays inside the stroke");
    }
    // published (applsci-4655781, Sec. 3.4): 47.0, 36.2, 12.9, 7.8 kN from the simulation
    o.t_ramp = 0.0001;
    near(simulate(kin, o).dF_max / 1000, 47.0, 0.4, "simulated excess force, instantaneous opening");
    o.t_ramp = 0.5;
    near(simulate(kin, o).dF_max / 1000, 7.8, 0.1, "simulated excess force, 0.5 s ramp");
    o.flow_lpm = 14.8;
    o.t_ramp = 1.0;
    near(simulate(kin, o).dF_max / 1000, 5.17, 0.06, "simulated excess force, 14.8 L/min and 1 s");
    // damping lowers the excess over the quasi-steady force (inertia and oscillation), so zeta = 0 is the bound
    for (double tr : {0.1, 0.5, 1.0}) {
        o.t_ramp = tr;
        o.damping_ratio = 0.0;
        const double d0 = simulate(kin, o).dF_max_osc;
        o.damping_ratio = 0.1;
        check(simulate(kin, o).dF_max_osc < d0, "damping lowers the excess force");
    }
    o.damping_ratio = 0.0;
    // with damping the oscillation after the ramp decays: compare the swing over the last period
    o.flow_lpm = 11.25;
    o.t_ramp = 0.0001;
    auto swing = [&](double zeta) {
        o.damping_ratio = zeta;
        auto r = simulate(kin, o);
        double lo = 1e300, hi = -1e300;
        for (std::size_t i = r.history.size() * 5 / 6; i < r.history.size(); ++i) {
            lo = std::min(lo, r.history[i].dF);
            hi = std::max(hi, r.history[i].dF);
        }
        return hi - lo;
    };
    check(swing(0.2) < 0.1 * swing(0.0), "damping makes the oscillation decay");
}

void test_map() {
    Model m = Model::from_file(SRC + "/examples/horty/tower.mech");
    Kinematics kin(m);
    MapOptions o;
    o.ax = 1;
    o.ay = 3;
    o.report = 1;
    o.fixed = kin.reference().lengths;
    o.fixed[0] = 973;
    o.x0 = 758; o.x1 = 1198; o.nx = 5;
    o.y0 = 4942.5; o.y1 = 10721.5; o.ny = 3;
    o.flow_lpm = 56.1;
    auto cells = load_map(kin, o);
    check(cells.size() == 15, "map size");
    near(cells.front().p, 27.21713, 0.01, "map corner: mast horizontal, hammer lowest (Fig. 8a)");
    near(cells[10].p, 189.3, 0.05, "map corner: mast horizontal, hammer highest (Table 6)");
    check(std::isinf(cells[10].tr_ret), "static load above 140 bar gives no admissible ramp");
    // the local ramp time at the published nominal point is below the whole-stroke value 3.84 s
    check(cells[5].tr_ret <= 3.85 && cells[5].tr_ret > 0.5, "local ramp time at s_b = 2.89 m, mast horizontal");
}

void test_sensitivity() {
    Uniform u(11);
    std::vector<double> x1, x2, x3, y;
    for (int i = 0; i < 20000; ++i) {
        const double a = u(), b = u(), c = u();
        x1.push_back(a); x2.push_back(b); x3.push_back(c);
        y.push_back(3 * a + b);  // Var: 9/12 and 1/12 -> S1 = 0.9, S2 = 0.1, S3 = 0
    }
    near(first_order_index(x1, y), 0.9, 0.02, "first-order index of a dominant input");
    near(first_order_index(x2, y), 0.1, 0.02, "first-order index of a minor input");
    near(first_order_index(x3, y), 0.0, 0.01, "first-order index of an inert input");
    double r2 = 0;
    auto b = standardized_regression({x1, x2, x3}, y, r2);
    near(r2, 1.0, 1e-9, "R^2 of a linear model");
    near(b[0], 3 / std::sqrt(10.0), 0.01, "SRC of x1");
    near(b[1], 1 / std::sqrt(10.0), 0.01, "SRC of x2");
    near(b[2], 0.0, 0.01, "SRC of x3");
}

void test_figures() {
    double a0, a1, st;
    nice_axis(0.3, 108.42, 6, a0, a1, st);
    check(a0 == 0 && st == 20 && a1 == 120, "axis ends on a labelled tick (0..120 by 20)");
    nice_axis(-90.6, 4.02, 6, a0, a1, st);
    check(a0 == -100 && a1 == 20 && st == 20, "negative axis");
    Figure f;
    Panel p;
    p.x.label = "s_b (m)";
    p.letter = "(a)";
    Series s;
    s.x = {0, 1, 2};
    s.y = {1, 4, 9};
    s.label = "test";
    p.series.push_back(s);
    f.panels.push_back(p);
    const std::string svg = render_svg(f);
    check(svg.find("baseline-shift=\"sub\"") != std::string::npos, "subscript markup");
    check(svg.find("#000") != std::string::npos && svg.find("rgb(") == std::string::npos, "black only");
    Model m = Model::from_file(SRC + "/examples/boom/boom.mech");
    Kinematics kin(m);
    const std::string d = draw_mechanism(kin, {kin.reference()}, 100, 80);
    check(d.find("<polygon") != std::string::npos && d.find("cyl") != std::string::npos, "mechanism drawing");
    // solid-model parts: every integrated face is also sampled for drawing, and drawn with the link
    StepModel sg = read_step(SRC + "/tests/data/revolved_solids.step", "", true, true);
    std::size_t f_int = 0, f_draw = 0, n_edges = 0;
    DrawOptions o;
    o.width_mm = 100; o.height_mm = 80;
    for (const auto& oc : sg.occurrences) {
        f_int += oc.faces;
        for (const auto& so : oc.solids) {
            f_draw += so.faces.size();
            n_edges += so.edges.size();
            check(so.curved.size() == so.faces.size(), "curved flag per face");
            DrawPart dp;
            dp.body = 1;
            for (const auto& fc : so.faces) {
                std::vector<std::vector<Vec2>> loops;
                for (const auto& lp : fc) {
                    std::vector<Vec2> l2;
                    for (const auto& q : lp) l2.push_back({q[0] * 0.01, q[1] * 0.01});
                    loops.push_back(l2);
                }
                dp.faces.push_back(loops);
            }
            dp.curved = so.curved;
            o.parts.push_back(dp);
        }
    }
    check(f_draw == f_int && f_draw > 0, "faces sampled for drawing = faces integrated");
    check(n_edges > 0, "edges sampled for drawing");
    const std::string d2 = draw_mechanism(kin, {kin.reference()}, o);
    check(d2.find("fill-rule=\"evenodd\"") != std::string::npos, "mechanism drawing with solid-model parts");
}

// 7. Input validation.
void test_validation() {
    bool threw = false;
    try { Model::from_string("body a b\nrevolute A ground a 0 0\nrevolute B a b 1 0\n"); } catch (const std::exception&) { threw = true; }
    check(threw, "under-actuated model rejected");
    threw = false;
    try { Model::from_string("body a\nrevolute A ground a 0 0\nactuator c ground 1 0 a 0 1 bore=40 rod=50 min=1 max=2\n"); }
    catch (const std::exception&) { threw = true; }
    check(threw, "rod larger than bore rejected");
    threw = false;
    try {
        Model m = Model::from_string("body a\nrevolute A ground a 0 0\nactuator c ground 100 0 a 0 100\nmass x a 50 50 1\n");
        Kinematics k(m);
        k.solve(k.reference(), {500.0});
    } catch (const std::exception&) { threw = true; }
    check(threw, "unreachable length reported");
}
}  // namespace

int main() {
    run("synthetic boom vs closed form", test_boom_closed_form);
    run("external force on the boom vs closed form", test_external_force);
    run("inclined slide vs closed form", test_incline);
    run("hydraulic formulas", test_hydraulics);
    run("Horty tower, outline variant vs. Python implementation", test_horty_outline);
    run("Horty post lift, outline variant vs. Python implementation", test_post_lift_outline);
    run("Horty, solid-model variant vs. published results", test_horty_solid);
    run("stone picker: Mathcad model and published pressures", test_usarya);
    run("Monte Carlo helpers", test_monte_carlo);
    run("STEP reader", test_step);
    run("STEP mass properties of solids of revolution", test_step_revolved);
    run("STEP mass properties of solids with B-spline faces", test_step_bspline);
    run("STEP length units (mm, inch, m, missing)", test_step_units);
    run("open example: synthetic backhoe from STEP to forces", test_excavator);
    run("Monte Carlo bootstrap confidence intervals", test_bootstrap);
    run("input validation", test_validation);
    run("time-domain simulation vs. ramp formula and published values", test_simulation);
    run("load and ramp-time map", test_map);
    run("sensitivity indices", test_sensitivity);
    run("figures (SVG)", test_figures);
    std::printf("%d checks, %d failed\n", g_checks, g_fail);
    return g_fail == 0 ? 0 : 1;
}
