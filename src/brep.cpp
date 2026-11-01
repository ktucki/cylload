// SPDX-License-Identifier: MIT
// Volume, first and second moments of B-rep solids from a STEP file by the
// divergence theorem:  int_V phi dV = sum_faces int_S (Phi, 0, 0) . n dA,  dPhi/dx = phi.
//
// Analytic faces (plane, cylinder, cone, sphere, torus) are integrated in their
// parameter domain (u, v). With G(u, v) = int_{v0}^{v} Phi(S) (S_u x S_v)_x ds, Green's
// theorem turns the face integral into the loop integral  -oint G du  over the face
// bounds, which are traversed exactly as stored in the file (the orientation of the
// bounds already encodes the face side, so ADVANCED_FACE.same_sense is not needed).
// u is the periodic (angular) parameter, so the seam of a periodic face does not
// contribute (du = 0 along it). For a cone reaching its apex or a sphere face reaching
// a pole, v0 is placed at the degenerate point, where G = 0 (the apex or pole is an
// implicit edge of the parameter domain). Tori with a negative major radius (barrel
// rollers) need no special treatment. Edges (lines, circles, ellipses, B-spline and
// rational B-spline curves) are integrated with Gauss-Legendre quadrature in their own
// parameter. Other faces (B-spline surfaces) are replaced by a fan of triangles spanning
// each bound. Void shells of BREP_WITH_VOIDS always remove material.
// The volume is also obtained with the fields (0, y, 0) and (0, 0, z), and the second
// moments with a second field each; the differences are reported as consistency checks.
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <stdexcept>
#include <string>

#include "step_entities.hpp"

namespace cylload {
namespace step_detail {

namespace {

constexpr double kPi = 3.14159265358979323846;
using V3 = std::array<double, 3>;

V3 add(const V3& a, const V3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
V3 sub(const V3& a, const V3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
V3 scl(const V3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
double dot(const V3& a, const V3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
V3 crs(const V3& a, const V3& b) { return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; }

// 8-point Gauss-Legendre on [-1, 1]
const double GX[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
                      0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363};
const double GW[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
                      0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};

// ------------------------------------------------------------------ argument values
struct Val {
    char t = '$';  // L list, R ref, N number, S string, E enum, $ null, * derived, T typed
    double n = 0;
    long r = 0;
    std::string s;
    std::vector<Val> l;
};

void ws(const std::string& s, std::size_t& p) {
    while (p < s.size() && std::isspace(static_cast<unsigned char>(s[p]))) ++p;
}

Val parse_val(const std::string& s, std::size_t& p);

std::vector<Val> parse_list_body(const std::string& s, std::size_t& p) {  // after '(' up to ')'
    std::vector<Val> out;
    ws(s, p);
    if (p < s.size() && s[p] == ')') { ++p; return out; }
    while (p < s.size()) {
        out.push_back(parse_val(s, p));
        ws(s, p);
        if (p < s.size() && s[p] == ',') { ++p; continue; }
        if (p < s.size() && s[p] == ')') { ++p; break; }
        throw std::runtime_error("STEP: malformed parameter list");
    }
    return out;
}

Val parse_val(const std::string& s, std::size_t& p) {
    ws(s, p);
    Val v;
    if (p >= s.size()) return v;
    const char c = s[p];
    if (c == '(') { ++p; v.t = 'L'; v.l = parse_list_body(s, p); return v; }
    if (c == '#') {
        ++p;
        v.t = 'R';
        while (p < s.size() && std::isdigit(static_cast<unsigned char>(s[p]))) v.r = v.r * 10 + (s[p++] - '0');
        return v;
    }
    if (c == '\'') {
        v.t = 'S';
        ++p;
        while (p < s.size()) {
            if (s[p] == '\'') {
                if (p + 1 < s.size() && s[p + 1] == '\'') { v.s += '\''; p += 2; continue; }
                ++p;
                break;
            }
            v.s += s[p++];
        }
        return v;
    }
    if (c == '.') {  // enumeration .T. .F. .UNSPECIFIED.
        v.t = 'E';
        auto q = s.find('.', p + 1);
        v.s = s.substr(p + 1, q - p - 1);
        p = q + 1;
        return v;
    }
    if (c == '$') { ++p; v.t = '$'; return v; }
    if (c == '*') { ++p; v.t = '*'; return v; }
    if (std::isalpha(static_cast<unsigned char>(c))) {  // typed parameter, e.g. LENGTH_MEASURE(1.)
        std::size_t q = p;
        while (q < s.size() && (std::isalnum(static_cast<unsigned char>(s[q])) || s[q] == '_')) ++q;
        v.t = 'T';
        v.s = s.substr(p, q - p);
        p = q;
        ws(s, p);
        if (p < s.size() && s[p] == '(') { ++p; v.l = parse_list_body(s, p); }
        return v;
    }
    v.t = 'N';
    char* end = nullptr;
    v.n = std::strtod(s.c_str() + p, &end);
    p = static_cast<std::size_t>(end - s.c_str());
    return v;
}

std::vector<Val> args_of(const Entity& e) {
    std::size_t p = 0;
    std::string body = e.args + ")";
    return parse_list_body(body, p);
}

/// complex entity "A(..) B(..)" -> map of partial entity name to its parameters
std::map<std::string, std::vector<Val>> complex_parts(const Entity& e) {
    std::map<std::string, std::vector<Val>> out;
    const std::string& s = e.args;
    std::size_t p = 0;
    while (true) {
        ws(s, p);
        if (p >= s.size()) break;
        std::size_t q = p;
        while (q < s.size() && (std::isalnum(static_cast<unsigned char>(s[q])) || s[q] == '_')) ++q;
        if (q == p) break;
        std::string name = s.substr(p, q - p);
        p = q;
        ws(s, p);
        if (p >= s.size() || s[p] != '(') break;
        ++p;
        out[name] = parse_list_body(s, p);
    }
    return out;
}

double num(const Val& v) {
    if (v.t == 'N') return v.n;
    if (v.t == 'T' && !v.l.empty()) return num(v.l[0]);
    throw std::runtime_error("STEP: expected a number");
}

// ------------------------------------------------------------------ geometry access
struct Geo {
    const Entities& e;
    std::unordered_map<long, V3> pts;

    const Entity& ent(long id) const {
        const Entity* x = get(e, id);
        if (!x) throw std::runtime_error("STEP: missing entity #" + std::to_string(id));
        return *x;
    }
    V3 point(long id) {
        auto it = pts.find(id);
        if (it != pts.end()) return it->second;
        auto v = numbers(ent(id).args);
        V3 p{v.size() > 0 ? v[0] : 0, v.size() > 1 ? v[1] : 0, v.size() > 2 ? v[2] : 0};
        pts.emplace(id, p);
        return p;
    }
    V3 vertex(long id) {  // VERTEX_POINT -> CARTESIAN_POINT
        auto r = refs(ent(id).args);
        return point(r.at(0));
    }
    struct Frame { V3 o, x, y, z; };
    Frame frame(long id) {
        Mat4 m = placement(e, id);
        return {{m[3], m[7], m[11]}, {m[0], m[4], m[8]}, {m[1], m[5], m[9]}, {m[2], m[6], m[10]}};
    }
};

// ------------------------------------------------------------------ curves
struct Nurbs {
    int p = 0;
    std::vector<V3> P;
    std::vector<double> w, U;
    V3 eval(double t) const {
        const int n = static_cast<int>(P.size());
        t = std::min(std::max(t, U[p]), U[n]);
        int k = p;
        while (k < n - 1 && t >= U[k + 1]) ++k;
        std::vector<std::array<double, 4>> d(p + 1);
        for (int j = 0; j <= p; ++j) {
            const auto& c = P[k - p + j];
            const double wj = w[k - p + j];
            d[j] = {c[0] * wj, c[1] * wj, c[2] * wj, wj};
        }
        for (int r = 1; r <= p; ++r)
            for (int j = p; j >= r; --j) {
                const double den = U[j + 1 + k - r] - U[j + k - p];
                const double a = den > 0 ? (t - U[j + k - p]) / den : 0.0;
                for (int c = 0; c < 4; ++c) d[j][c] = (1 - a) * d[j - 1][c] + a * d[j][c];
            }
        return {d[p][0] / d[p][3], d[p][1] / d[p][3], d[p][2] / d[p][3]};
    }
    double t0() const { return U[p]; }
    double t1() const { return U[P.size()]; }
    int spans() const {
        int s = 0;
        for (std::size_t i = p; i < P.size(); ++i) s += U[i + 1] > U[i];
        return std::max(1, s);
    }
    double project(const V3& q) const {  // parameter of the closest point
        const int N = 64 * spans();
        double best = 1e300, tb = t0();
        for (int i = 0; i <= N; ++i) {
            const double t = t0() + (t1() - t0()) * i / N;
            const double d = dot(sub(eval(t), q), sub(eval(t), q));
            if (d < best) { best = d; tb = t; }
        }
        double a = std::max(t0(), tb - (t1() - t0()) / N), b = std::min(t1(), tb + (t1() - t0()) / N);
        for (int it = 0; it < 80; ++it) {
            const double m1 = a + (b - a) * 0.381966, m2 = a + (b - a) * 0.618034;
            const double d1 = dot(sub(eval(m1), q), sub(eval(m1), q)), d2 = dot(sub(eval(m2), q), sub(eval(m2), q));
            (d1 < d2 ? b : a) = d1 < d2 ? m2 : m1;
        }
        return 0.5 * (a + b);
    }
};

Nurbs nurbs_from(const std::vector<Val>& bsc, const std::vector<Val>& knots_part, const std::vector<Val>* weights, Geo& g) {
    // bsc: (degree, control points, form, closed, self_intersect); knots_part: (mults, knots, spec)
    Nurbs c;
    c.p = static_cast<int>(num(bsc.at(0)));
    for (const auto& r : bsc.at(1).l) c.P.push_back(g.point(r.r));
    const auto& mult = knots_part.at(0).l;
    const auto& kn = knots_part.at(1).l;
    for (std::size_t i = 0; i < kn.size(); ++i)
        for (int m = 0; m < static_cast<int>(num(mult.at(i))); ++m) c.U.push_back(num(kn[i]));
    c.w.assign(c.P.size(), 1.0);
    if (weights)
        for (std::size_t i = 0; i < c.P.size() && i < weights->size(); ++i) c.w[i] = num((*weights)[i]);
    if (c.U.size() != c.P.size() + c.p + 1) throw std::runtime_error("STEP: inconsistent B-spline knot vector");
    return c;
}

/// An edge as traversed in a loop: P(s), s in [0, 1] from its start to its end vertex.
struct EdgePath {
    int kind = 0;  // 0 line, 1 conic (circle/ellipse), 2 nurbs
    V3 a{}, b{};   // line end points
    Geo::Frame f{};
    double ra = 0, rb = 0, t0 = 0, t1 = 0;  // conic radii and angles, or nurbs parameters
    Nurbs nb;
    int nsub = 1;
    V3 at(double s) const {
        if (kind == 0) return add(a, scl(sub(b, a), s));
        const double t = t0 + (t1 - t0) * s;
        if (kind == 1) return add(f.o, add(scl(f.x, ra * std::cos(t)), scl(f.y, rb * std::sin(t))));
        return nb.eval(t);
    }
};

double angle_in(const Geo::Frame& f, const V3& q, double ra, double rb) {
    const V3 d = sub(q, f.o);
    return std::atan2(dot(d, f.y) / rb, dot(d, f.x) / ra);
}

EdgePath edge_path(long oriented_edge, Geo& g) {
    const auto oe = args_of(g.ent(oriented_edge));          // (name, *, *, edge_element, orientation)
    const bool o_fwd = oe.at(4).s == "T";
    const auto ec = args_of(g.ent(oe.at(3).r));             // (name, v_start, v_end, geometry, same_sense)
    const bool same = ec.at(4).s == "T";
    const V3 vs = g.vertex(ec.at(1).r), ve = g.vertex(ec.at(2).r);
    const bool closed = ec.at(1).r == ec.at(2).r;
    const bool along = (same == o_fwd);                      // traversal follows the curve parameter
    const V3 start = o_fwd ? vs : ve, end = o_fwd ? ve : vs;
    EdgePath ep;
    long cid = ec.at(3).r;
    const Entity* ce = &g.ent(cid);
    if (ce->type == "SURFACE_CURVE" || ce->type == "SEAM_CURVE") {
        cid = args_of(*ce).at(1).r;
        ce = &g.ent(cid);
    }
    if (ce->type == "LINE" || ce->type == "POLYLINE") {
        ep.kind = 0;
        ep.a = start;
        ep.b = end;
        ep.nsub = 2;
        return ep;
    }
    if (ce->type == "CIRCLE" || ce->type == "ELLIPSE") {
        auto a = args_of(*ce);
        ep.kind = 1;
        ep.f = g.frame(a.at(1).r);
        ep.ra = num(a.at(2));
        ep.rb = ce->type == "ELLIPSE" ? num(a.at(3)) : ep.ra;
        const double ts = angle_in(ep.f, start, ep.ra, ep.rb), te = angle_in(ep.f, end, ep.ra, ep.rb);
        double d;
        if (along) {
            d = std::fmod(te - ts + 4 * kPi, 2 * kPi);
            if (closed || d < 1e-12) d = 2 * kPi;
        } else {
            d = -std::fmod(ts - te + 4 * kPi, 2 * kPi);
            if (closed || d > -1e-12) d = -2 * kPi;
        }
        ep.t0 = ts;
        ep.t1 = ts + d;
        ep.nsub = std::max(2, static_cast<int>(std::ceil(std::fabs(d) / (kPi / 12))));
        return ep;
    }
    Nurbs nb;
    if (ce->type == "B_SPLINE_CURVE_WITH_KNOTS") {
        auto a = args_of(*ce);  // name, degree, pts, form, closed, self, mults, knots, spec
        std::vector<Val> bsc(a.begin() + 1, a.begin() + 6), kp(a.begin() + 6, a.end());
        nb = nurbs_from(bsc, kp, nullptr, g);
    } else if (ce->type.empty()) {
        auto parts = complex_parts(*ce);
        if (!parts.count("B_SPLINE_CURVE") || !parts.count("B_SPLINE_CURVE_WITH_KNOTS"))
            throw std::runtime_error("STEP: unsupported complex edge curve #" + std::to_string(cid));
        const std::vector<Val>* w = nullptr;
        if (parts.count("RATIONAL_B_SPLINE_CURVE")) w = &parts["RATIONAL_B_SPLINE_CURVE"].at(0).l;
        nb = nurbs_from(parts["B_SPLINE_CURVE"], parts["B_SPLINE_CURVE_WITH_KNOTS"], w, g);
    } else {
        throw std::runtime_error("STEP: unsupported edge curve type " + ce->type);
    }
    ep.kind = 2;
    double ts = nb.project(start), te = nb.project(end);
    if (closed) { ts = along ? nb.t0() : nb.t1(); te = along ? nb.t1() : nb.t0(); }
    ep.t0 = ts;
    ep.t1 = te;
    ep.nsub = std::max(8, 2 * nb.spans());
    ep.nb = std::move(nb);
    return ep;
}

// ------------------------------------------------------------------ surfaces
struct Surf {
    int kind = -1;  // 0 plane, 1 cylinder, 2 cone, 3 sphere, 4 torus, -1 other
    Geo::Frame f{};
    double R = 0, r = 0, ta = 0;
    bool periodic_v() const { return kind == 4; }
    bool periodic_u() const { return kind >= 1; }
    std::array<double, 2> uv(const V3& q) const {
        const V3 d = sub(q, f.o);
        const double x = dot(d, f.x), y = dot(d, f.y), z = dot(d, f.z);
        switch (kind) {
            case 0: return {x, y};
            case 1:
            case 2: return {std::atan2(y, x), z};
            case 3: return {std::atan2(y, x), std::atan2(z, std::hypot(x, y))};
            default: return {std::atan2(y, x), std::atan2(z, std::hypot(x, y) - R)};
        }
    }
    void eval(double u, double v, V3& S, V3& N) const {
        const double cu = std::cos(u), su = std::sin(u);
        const V3 e = add(scl(f.x, cu), scl(f.y, su)), ep = add(scl(f.x, -su), scl(f.y, cu));
        V3 Su, Sv;
        switch (kind) {
            case 0:
                S = add(f.o, add(scl(f.x, u), scl(f.y, v)));
                N = f.z;
                return;
            case 1:
                S = add(f.o, add(scl(e, R), scl(f.z, v)));
                Su = scl(ep, R);
                Sv = f.z;
                break;
            case 2: {
                const double rho = R + v * ta;
                S = add(f.o, add(scl(e, rho), scl(f.z, v)));
                Su = scl(ep, rho);
                Sv = add(scl(e, ta), f.z);
                break;
            }
            case 3: {
                const double cv = std::cos(v), sv = std::sin(v);
                S = add(f.o, scl(add(scl(e, cv), scl(f.z, sv)), R));
                Su = scl(ep, R * cv);
                Sv = scl(add(scl(e, -sv), scl(f.z, cv)), R);
                break;
            }
            default: {
                const double cv = std::cos(v), sv = std::sin(v), rho = R + r * cv;
                S = add(f.o, add(scl(e, rho), scl(f.z, r * sv)));
                Su = scl(ep, rho);
                Sv = scl(add(scl(e, -sv), scl(f.z, cv)), r);
                break;
            }
        }
        N = crs(Su, Sv);
    }
};

Surf surface(long id, Geo& g) {
    const Entity& s = g.ent(id);
    Surf out;
    if (s.type == "PLANE") { out.kind = 0; out.f = g.frame(args_of(s).at(1).r); }
    else if (s.type == "CYLINDRICAL_SURFACE") { auto a = args_of(s); out.kind = 1; out.f = g.frame(a.at(1).r); out.R = num(a.at(2)); }
    else if (s.type == "CONICAL_SURFACE") {
        auto a = args_of(s);
        out.kind = 2; out.f = g.frame(a.at(1).r); out.R = num(a.at(2)); out.ta = std::tan(num(a.at(3)));
    } else if (s.type == "SPHERICAL_SURFACE") { auto a = args_of(s); out.kind = 3; out.f = g.frame(a.at(1).r); out.R = num(a.at(2)); }
    else if (s.type == "TOROIDAL_SURFACE") {
        auto a = args_of(s);
        out.kind = 4; out.f = g.frame(a.at(1).r); out.R = num(a.at(2)); out.r = num(a.at(3));
    }
    return out;
}

double unwrap(double x, double ref) { return x + 2 * kPi * std::round((ref - x) / (2 * kPi)); }

// integrand: Phi_k(S - ref) * N_x for the 10 volume integrals
void integrand(const V3& S, const V3& N, const V3& ref, double* out, double scale) {
    const double x = S[0] - ref[0], y = S[1] - ref[1], z = S[2] - ref[2], nx = N[0] * scale;
    out[10] += y * N[1] * scale;  // the volume once more with fields (0, y, 0) and (0, 0, z):
    out[11] += z * N[2] * scale;  // a consistency check of the boundary representation
    out[12] += x * x * y * N[1] * scale;      // int x^2 dV with the field (0, x^2 y, 0)
    out[13] += y * y * y / 3.0 * N[1] * scale; // int y^2 dV with (0, y^3/3, 0)
    out[14] += z * z * z / 3.0 * N[2] * scale; // int z^2 dV with (0, 0, z^3/3)
    out[0] += x * nx;
    out[1] += 0.5 * x * x * nx;
    out[2] += x * y * nx;
    out[3] += x * z * nx;
    out[4] += x * x * x / 3.0 * nx;
    out[5] += x * y * y * nx;
    out[6] += x * z * z * nx;
    out[7] += 0.5 * x * x * y * nx;
    out[8] += x * y * z * nx;
    out[9] += 0.5 * x * x * z * nx;
}

struct Loop { std::vector<EdgePath> edges; };

// -oint G du with G(u, v) = int_{v0}^{v} w(u, s) ds
void face_green(const Surf& sf, const std::vector<Loop>& loops, double v0, const V3& ref, double* acc) {
    auto G = [&](double u, double v, double* out) {
        const double span = v - v0;
        if (span == 0) return;
        const bool ang = sf.kind == 3 || sf.kind == 4;
        const int nsub = ang ? std::max(1, static_cast<int>(std::ceil(std::fabs(span) / (kPi / 8)))) : 1;
        const double h = span / nsub;
        for (int k = 0; k < nsub; ++k) {
            const double a = v0 + k * h;
            for (int i = 0; i < 8; ++i) {
                const double s = a + 0.5 * h * (GX[i] + 1);
                V3 S, N;
                sf.eval(u, s, S, N);
                integrand(S, N, ref, out, 0.5 * h * GW[i]);
            }
        }
    };
    // torus: put all bounds on one branch of the periodic v around the circular mean of the
    // boundary (the face lies in a band narrower than a full turn of the tube)
    double vc = 0.0;
    if (sf.periodic_v()) {
        double cs = 0, sn = 0;
        for (const auto& lp : loops)
            for (const auto& ep : lp.edges)
                for (int i = 0; i < 4 * ep.nsub; ++i) {
                    const double v = sf.uv(ep.at(static_cast<double>(i) / (4 * ep.nsub)))[1];
                    cs += std::cos(v);
                    sn += std::sin(v);
                }
        vc = std::atan2(sn, cs);
    }
    for (const auto& lp : loops) {
        if (lp.edges.empty()) continue;
        auto prev = sf.uv(lp.edges.front().at(0.0));
        if (sf.periodic_v()) prev[1] = unwrap(prev[1], vc);
        for (const auto& ep : lp.edges) {
            for (int k = 0; k < ep.nsub; ++k) {
                const double sa = static_cast<double>(k) / ep.nsub, sb = static_cast<double>(k + 1) / ep.nsub;
                for (int i = 0; i < 8; ++i) {
                    const double s = sa + 0.5 * (sb - sa) * (GX[i] + 1);
                    auto q = sf.uv(ep.at(s));
                    if (sf.periodic_u()) q[0] = unwrap(q[0], prev[0]);
                    if (sf.periodic_v()) q[1] = unwrap(q[1], prev[1]);
                    const double hs = 1e-6 * (sb - sa);
                    auto qp = sf.uv(ep.at(std::min(1.0, s + hs))), qm = sf.uv(ep.at(std::max(0.0, s - hs)));
                    if (sf.periodic_u()) { qp[0] = unwrap(qp[0], q[0]); qm[0] = unwrap(qm[0], q[0]); }
                    const double dsp = std::min(1.0, s + hs) - std::max(0.0, s - hs);
                    const double dudt = (qp[0] - qm[0]) / dsp;
                    prev = q;
                    if (dudt == 0) continue;
                    double g[15] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
                    G(q[0], q[1], g);
                    const double wgt = -dudt * 0.5 * (sb - sa) * GW[i];
                    for (int c = 0; c < 15; ++c) acc[c] += g[c] * wgt;
                }
            }
        }
    }
}

/// total winding of the loops in the periodic coordinates (in turns)
std::array<double, 2> windings(const Surf& sf, const std::vector<Loop>& loops) {
    std::array<double, 2> w{0, 0};
    for (const auto& lp : loops) {
        if (lp.edges.empty()) continue;
        const auto first = sf.uv(lp.edges.front().at(0.0));
        auto prev = first;
        for (const auto& ep : lp.edges) {
            const int n = 8 * ep.nsub;
            for (int i = 1; i <= n; ++i) {
                auto q = sf.uv(ep.at(static_cast<double>(i) / n));
                q[0] = unwrap(q[0], prev[0]);
                q[1] = sf.periodic_v() ? unwrap(q[1], prev[1]) : q[1];
                prev = q;
            }
        }
        w[0] += (prev[0] - first[0]) / (2 * kPi);
        w[1] += (prev[1] - first[1]) / (2 * kPi);
    }
    return {std::round(w[0]), std::round(w[1])};
}

// fan of triangles over every bound (exact for planar bounds made of straight edges)
void face_fan(const std::vector<Loop>& loops, const V3& ref, double* acc) {
    static const double A[7][4] = {{1.0 / 3, 1.0 / 3, 1.0 / 3, 0.225},
                                   {0.059715871789770, 0.470142064105115, 0.470142064105115, 0.132394152788506},
                                   {0.470142064105115, 0.059715871789770, 0.470142064105115, 0.132394152788506},
                                   {0.470142064105115, 0.470142064105115, 0.059715871789770, 0.132394152788506},
                                   {0.797426985353087, 0.101286507323456, 0.101286507323456, 0.125939180544827},
                                   {0.101286507323456, 0.797426985353087, 0.101286507323456, 0.125939180544827},
                                   {0.101286507323456, 0.101286507323456, 0.797426985353087, 0.125939180544827}};
    for (const auto& lp : loops) {
        std::vector<V3> pts;
        for (const auto& ep : lp.edges) {
            const int n = 4 * ep.nsub;
            for (int i = 0; i < n; ++i) pts.push_back(ep.at(static_cast<double>(i) / n));
        }
        if (pts.size() < 3) continue;
        V3 c{0, 0, 0};
        for (const auto& p : pts) c = add(c, p);
        c = scl(c, 1.0 / static_cast<double>(pts.size()));
        for (std::size_t i = 0; i < pts.size(); ++i) {
            const V3& a = pts[i];
            const V3& b = pts[(i + 1) % pts.size()];
            const V3 N = crs(sub(a, c), sub(b, c));  // twice the area times the normal
            for (const auto& q : A) {
                const V3 S = add(add(scl(c, q[0]), scl(a, q[1])), scl(b, q[2]));
                integrand(S, N, ref, acc, 0.5 * q[3]);
            }
        }
    }
}

void integrate_face(long face_id, Geo& g, VolumeIntegrals& out, double sign) {
    const auto fa = args_of(g.ent(face_id));  // (name, bounds, surface, same_sense)
    std::vector<Loop> loops;
    for (const auto& b : fa.at(1).l) {
        const auto ba = args_of(g.ent(b.r));     // (name, loop, orientation)
        const bool fwd = ba.at(2).s == "T";
        const Entity& le = g.ent(ba.at(1).r);
        Loop lp;
        if (le.type == "EDGE_LOOP") {
            const auto la = args_of(le);
            for (const auto& oe : la.at(1).l) lp.edges.push_back(edge_path(oe.r, g));
        } else if (le.type == "POLY_LOOP") {
            auto pl = args_of(le).at(1).l;
            for (std::size_t i = 0; i < pl.size(); ++i) {
                EdgePath ep;
                ep.kind = 0;
                ep.a = g.point(pl[i].r);
                ep.b = g.point(pl[(i + 1) % pl.size()].r);
                lp.edges.push_back(ep);
            }
        } else {
            continue;  // VERTEX_LOOP: a degenerate bound (cone apex), no contribution
        }
        if (!fwd) {
            std::reverse(lp.edges.begin(), lp.edges.end());
            for (auto& ep : lp.edges) {
                if (ep.kind == 0) std::swap(ep.a, ep.b);
                else std::swap(ep.t0, ep.t1);
            }
        }
        loops.push_back(std::move(lp));
    }
    double acc[15] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    Surf sf = surface(fa.at(2).r, g);
    bool fan = sf.kind < 0;
    if (!fan) {
        double v0 = 0.0;
        if (sf.kind >= 1) {
            auto w = windings(sf, loops);
            const bool same = fa.at(3).s == "T";
            if (sf.kind == 1 && w[0] != 0) fan = true;
            if (sf.kind == 4 && (w[0] != 0 || w[1] != 0)) fan = true;
            // G must vanish at a degenerate point (cone apex, sphere pole) that the face reaches,
            // because the apex/pole is an implicit edge of the parameter domain
            if (sf.kind == 2 && sf.ta != 0) v0 = -sf.R / sf.ta;
            if (sf.kind == 3) {
                if (w[0] != 0) {
                    v0 = ((w[0] > 0) == same ? 1.0 : -1.0) * kPi / 2;  // a pole inside the face
                } else {
                    double vmax = -1e300, vmin = 1e300;
                    for (const auto& lp : loops)
                        for (const auto& ep : lp.edges)
                            for (int i = 0; i <= 8 * ep.nsub; ++i) {
                                const double v = sf.uv(ep.at(static_cast<double>(i) / (8 * ep.nsub)))[1];
                                vmax = std::max(vmax, v);
                                vmin = std::min(vmin, v);
                            }
                    const bool north = vmax > kPi / 2 - 1e-6, south = vmin < -kPi / 2 + 1e-6;
                    if (north && south) fan = true;  // a lune through both poles
                    else if (north) v0 = kPi / 2;
                    else if (south) v0 = -kPi / 2;
                }
            }
        }
        if (!fan) face_green(sf, loops, v0, out.ref, acc);
    }
    if (fan) {
        face_fan(loops, out.ref, acc);
        ++out.faces_triangulated;
    }
    for (int c = 0; c < 15; ++c) out.v[c] += sign * acc[c];
    ++out.faces;
}

void integrate_shell(long shell_id, Geo& g, VolumeIntegrals& out, double sign) {
    const Entity* s = &g.ent(shell_id);
    if (s->type == "ORIENTED_CLOSED_SHELL") {
        auto a = args_of(*s);  // (name, *, shell, orientation)
        if (a.at(3).s == "F") sign = -sign;
        s = &g.ent(a.at(2).r);
    }
    const auto sa = args_of(*s);
    for (const auto& f : sa.at(1).l) {
        try {
            integrate_face(f.r, g, out, sign);
        } catch (const std::exception& ex) {
            throw std::runtime_error(std::string(ex.what()) + " (face #" + std::to_string(f.r) + ")");
        }
    }
}

}  // namespace

VolumeIntegrals integrate_solids(const Entities& e, const std::vector<long>& reps) {
    Geo g{e, {}};
    VolumeIntegrals out;
    std::set<long> solids;
    for (long rep : reps) {
        const Entity* r = get(e, rep);
        if (!r) continue;
        auto a = args_of(*r);
        if (a.size() < 2 || a[1].t != 'L') continue;
        for (const auto& it : a[1].l) {
            if (it.t != 'R') continue;
            const Entity* x = get(e, it.r);
            if (x && (x->type == "MANIFOLD_SOLID_BREP" || x->type == "BREP_WITH_VOIDS")) solids.insert(it.r);
        }
    }
    bool ref_set = false;
    for (long sid : solids) {
        const Entity& s = g.ent(sid);
        auto a = args_of(s);
        if (!ref_set) {  // reference point: any vertex of the first solid (keeps the moments well conditioned)
            const auto shell = args_of(g.ent(a.at(1).r));
            for (const auto& f : shell.at(1).l) {
                auto fa = args_of(g.ent(f.r));
                auto ba = args_of(g.ent(fa.at(1).l.at(0).r));
                const Entity& le = g.ent(ba.at(1).r);
                if (le.type != "EDGE_LOOP") continue;
                const auto la = args_of(le);
                auto oe = args_of(g.ent(la.at(1).l.at(0).r));
                auto ec = args_of(g.ent(oe.at(3).r));
                out.ref = g.vertex(ec.at(1).r);
                ref_set = true;
                break;
            }
        }
        VolumeIntegrals one;
        one.ref = out.ref;
        integrate_shell(a.at(1).r, g, one, 1.0);
        if (one.v[0] < 0)  // tolerate an inside-out outer shell
            for (int c = 0; c < 15; ++c) one.v[c] = -one.v[c];
        if (s.type == "BREP_WITH_VOIDS") {
            // a void always removes material; exporters differ in how they combine the
            // orientation flag with the face orientation of the void shell, so the sign is
            // taken from the enclosed volume itself
            for (const auto& v : a.at(2).l) {
                VolumeIntegrals vd;
                vd.ref = out.ref;
                integrate_shell(v.r, g, vd, 1.0);
                const double sv = vd.v[0] > 0 ? -1.0 : 1.0;
                for (int c = 0; c < 15; ++c) one.v[c] += sv * vd.v[c];
                one.faces += vd.faces;
                one.faces_triangulated += vd.faces_triangulated;
            }
        }
        for (int c = 0; c < 15; ++c) out.v[c] += one.v[c];
        out.faces += one.faces;
        out.faces_triangulated += one.faces_triangulated;
        ++out.solids;
    }
    return out;
}

std::vector<SolidGeometry> solid_geometry(const Entities& e, const std::vector<long>& reps) {
    Geo g{e, {}};
    std::vector<SolidGeometry> out;
    std::set<long> solids;
    for (long rep : reps) {
        const Entity* r = get(e, rep);
        if (!r) continue;
        auto a = args_of(*r);
        if (a.size() < 2 || a[1].t != 'L') continue;
        for (const auto& it : a[1].l) {
            if (it.t != 'R') continue;
            const Entity* x = get(e, it.r);
            if (x && (x->type == "MANIFOLD_SOLID_BREP" || x->type == "BREP_WITH_VOIDS")) solids.insert(it.r);
        }
    }
    auto sample = [](const EdgePath& ep, Polyline3& pl, bool first) {
        const int n = std::max(1, ep.nsub);
        for (int i = first ? 0 : 1; i <= n; ++i) pl.push_back(ep.at(static_cast<double>(i) / n));
    };
    for (long sid : solids) {
        SolidGeometry sg;
        std::set<long> seen_edges;
        std::vector<long> shells;
        auto a = args_of(g.ent(sid));
        shells.push_back(a.at(1).r);
        if (g.ent(sid).type == "BREP_WITH_VOIDS")
            for (const auto& v : a.at(2).l) shells.push_back(v.r);
        for (long sh : shells) {
            const Entity* s = &g.ent(sh);
            if (s->type == "ORIENTED_CLOSED_SHELL") s = &g.ent(args_of(*s).at(2).r);
            const auto shell_args = args_of(*s);
            for (const auto& f : shell_args.at(1).l) {
                try {
                    const auto fa = args_of(g.ent(f.r));
                    std::vector<Polyline3> loops;
                    for (const auto& b : fa.at(1).l) {
                        const auto ba = args_of(g.ent(b.r));
                        const Entity& le = g.ent(ba.at(1).r);
                        Polyline3 lp;
                        if (le.type == "EDGE_LOOP") {
                            const auto loop_args = args_of(le);
                            for (const auto& oe : loop_args.at(1).l) {
                                const EdgePath ep = edge_path(oe.r, g);
                                sample(ep, lp, lp.empty());
                                const long ec = args_of(g.ent(oe.r)).at(3).r;
                                if (seen_edges.insert(ec).second) {
                                    Polyline3 el;
                                    sample(ep, el, true);
                                    sg.edges.push_back(std::move(el));
                                }
                            }
                        } else if (le.type == "POLY_LOOP") {
                            const auto poly_args = args_of(le);
                            for (const auto& p : poly_args.at(1).l) lp.push_back(g.point(p.r));
                            for (std::size_t i = 0; i < lp.size(); ++i) sg.edges.push_back({lp[i], lp[(i + 1) % lp.size()]});
                        }
                        if (lp.size() >= 3) loops.push_back(std::move(lp));
                    }
                    if (!loops.empty()) {
                        sg.faces.push_back(std::move(loops));
                        const Entity& surf = g.ent(fa.at(2).r);
                        sg.curved.push_back(surf.type == "PLANE" ? 0 : 1);
                    }
                } catch (const std::exception& ex) {
                    // a face that cannot be sampled is left out of the drawing
                    if (std::getenv("CYLLOAD_DEBUG")) std::fprintf(stderr, "face #%ld: %s\n", f.r, ex.what());
                }
            }
        }
        out.push_back(std::move(sg));
    }
    return out;
}

}  // namespace step_detail
}  // namespace cylload
