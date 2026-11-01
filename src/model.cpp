// SPDX-License-Identifier: MIT
#include "cylload/model.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace cylload {

namespace {
constexpr double kPi = 3.14159265358979323846;

std::string strip_comment(const std::string& s) {
    auto p = s.find('#');
    return p == std::string::npos ? s : s.substr(0, p);
}

std::vector<std::string> split_ws(const std::string& s) {
    std::istringstream is(s);
    std::vector<std::string> out;
    std::string t;
    while (is >> t) out.push_back(t);
    return out;
}

std::vector<std::string> split_csv(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::string cur;
    bool quoted = false;
    for (char ch : s) {
        if (ch == '"') { quoted = !quoted; continue; }
        if (ch == sep && !quoted) { out.push_back(cur); cur.clear(); continue; }
        if (ch != '\r') cur += ch;
    }
    out.push_back(cur);
    for (auto& f : out) {
        auto a = f.find_first_not_of(" \t");
        auto b = f.find_last_not_of(" \t");
        f = a == std::string::npos ? std::string() : f.substr(a, b - a + 1);
    }
    return out;
}

double num(const std::string& s, const std::string& where) {
    try {
        std::size_t pos = 0;
        double v = std::stod(s, &pos);
        if (pos != s.size()) throw std::invalid_argument(s);
        return v;
    } catch (...) {
        throw std::runtime_error(where + ": expected a number, got '" + s + "'");
    }
}

std::map<std::string, std::string> keyvals(const std::vector<std::string>& t, std::size_t from, const std::string& where) {
    std::map<std::string, std::string> kv;
    for (std::size_t i = from; i < t.size(); ++i) {
        auto p = t[i].find('=');
        if (p == std::string::npos) throw std::runtime_error(where + ": expected key=value, got '" + t[i] + "'");
        kv[t[i].substr(0, p)] = t[i].substr(p + 1);
    }
    return kv;
}

std::string dir_of(const std::string& path) {
    auto p = path.find_last_of("/\\");
    return p == std::string::npos ? std::string(".") : path.substr(0, p);
}

Vec2 rot(double th, const Vec2& p) {
    const double c = std::cos(th), s = std::sin(th);
    return {c * p[0] - s * p[1], s * p[0] + c * p[1]};
}
Vec2 drot(double th, const Vec2& p) {  // d/dtheta of rot(th, p)
    const double c = std::cos(th), s = std::sin(th);
    return {-s * p[0] - c * p[1], c * p[0] - s * p[1]};
}
}  // namespace

double Cylinder::area_piston() const { return kPi * bore * bore / 4.0; }
double Cylinder::area_annulus() const { return kPi * (bore * bore - rod * rod) / 4.0; }

int Model::body_index(const std::string& name) const {
    for (std::size_t i = 0; i < bodies.size(); ++i)
        if (bodies[i] == name) return static_cast<int>(i);
    throw std::runtime_error("unknown body '" + name + "'");
}

int Model::actuator_index(const std::string& name) const {
    for (std::size_t i = 0; i < actuators.size(); ++i)
        if (actuators[i].name == name) return static_cast<int>(i);
    throw std::runtime_error("unknown actuator '" + name + "'");
}

int Model::n_constraints() const {
    return 2 * static_cast<int>(revolutes.size() + prismatics.size()) + static_cast<int>(actuators.size());
}

std::vector<MassItem> Model::cylinder_masses() const {
    std::vector<MassItem> out;
    for (const auto& a : actuators) {
        if (!a.has_cylinder || a.cyl.mass <= 0) continue;
        out.push_back({a.name + "/pin_a", a.name, a.body_a, a.pa, a.cyl.mass / 2, 0.0});
        out.push_back({a.name + "/pin_b", a.name, a.body_b, a.pb, a.cyl.mass / 2, 0.0});
    }
    return out;
}

std::vector<MassItem> Model::all_masses() const {
    auto out = masses;
    auto c = cylinder_masses();
    out.insert(out.end(), c.begin(), c.end());
    return out;
}

double Model::total_mass() const {
    double s = 0;
    for (const auto& m : all_masses()) s += m.m;
    return s;
}

void Model::validate() const {
    if (n_moving() < 1) throw std::runtime_error("model has no moving body");
    if (n_constraints() != n_coords()) {
        std::ostringstream os;
        os << "model is not exactly driven: " << n_moving() << " moving bodies give " << n_coords()
           << " coordinates, joints remove " << n_coords() - dof() << ", so " << dof()
           << " actuator(s) are needed but " << actuators.size() << " are defined";
        throw std::runtime_error(os.str());
    }
    for (const auto& a : actuators) {
        if (a.l_ref <= 0) throw std::runtime_error("actuator " + a.name + " has zero length");
        if (a.has_cylinder) {
            if (a.cyl.bore <= a.cyl.rod || a.cyl.rod <= 0) throw std::runtime_error("cylinder " + a.name + ": need bore > rod > 0");
            if (a.cyl.lmax <= a.cyl.lmin) throw std::runtime_error("cylinder " + a.name + ": need max > min");
        }
    }
    for (const auto& p : prismatics)
        if (std::hypot(p.u[0], p.u[1]) < 1e-12) throw std::runtime_error("prismatic " + p.name + ": zero direction");
}

void Model::load_masses_csv(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open mass table '" + path + "'");
    std::string line;
    std::vector<std::string> head;
    char sep = ',';
    int ln = 0;
    while (std::getline(f, line)) {
        ++ln;
        if (ln == 1 && line.size() >= 3 && static_cast<unsigned char>(line[0]) == 0xEF) line = line.substr(3);  // BOM
        auto s = strip_comment(line);
        if (s.find_first_not_of(" \t\r") == std::string::npos) continue;
        if (head.empty()) {
            if (s.find(';') != std::string::npos && s.find(',') == std::string::npos) sep = ';';
            head = split_csv(s, sep);
            continue;
        }
        auto v = split_csv(s, sep);
        auto col = [&](const std::string& key, bool req) -> std::string {
            for (std::size_t i = 0; i < head.size(); ++i)
                if (head[i] == key) return i < v.size() ? v[i] : std::string();
            if (req) throw std::runtime_error(path + ": missing column '" + key + "'");
            return std::string();
        };
        const std::string where = path + ":" + std::to_string(ln);
        MassItem mi;
        mi.name = col("name", true);
        mi.body = body_index(col("body", true));
        mi.c = {num(col("x_mm", true), where), num(col("y_mm", true), where)};
        mi.m = num(col("mass_kg", true), where);
        auto j = col("J_kgm2", false);
        mi.j = j.empty() ? 0.0 : num(j, where);
        mi.group = col("group", false);
        if (mi.group.empty()) mi.group = bodies[mi.body];
        masses.push_back(mi);
    }
}

Model Model::from_file(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open model file '" + path + "'");
    std::stringstream ss;
    ss << f.rdbuf();
    return from_string(ss.str(), dir_of(path));
}

Model Model::from_string(const std::string& text, const std::string& base_dir) {
    Model m;
    std::istringstream is(text);
    std::string line;
    int ln = 0;
    auto point_of = [&](const std::vector<std::string>& t, std::size_t i, const std::string& w) {
        return Vec2{num(t.at(i), w), num(t.at(i + 1), w)};
    };
    while (std::getline(is, line)) {
        ++ln;
        auto t = split_ws(strip_comment(line));
        if (t.empty()) continue;
        const std::string w = "line " + std::to_string(ln);
        const std::string& k = t[0];
        try {
            if (k == "title") {
                m.title = line.substr(line.find("title") + 5);
                m.title = m.title.substr(m.title.find_first_not_of(" \t") == std::string::npos ? 0 : m.title.find_first_not_of(" \t"));
            } else if (k == "gravity") {
                m.gravity = point_of(t, 1, w);
            } else if (k == "body") {
                for (std::size_t i = 1; i < t.size(); ++i) {
                    if (std::find(m.bodies.begin(), m.bodies.end(), t[i]) != m.bodies.end())
                        throw std::runtime_error("body '" + t[i] + "' defined twice");
                    m.bodies.push_back(t[i]);
                }
            } else if (k == "revolute") {  // revolute NAME BODY_A BODY_B X Y
                if (t.size() != 6) throw std::runtime_error("usage: revolute NAME BODY_A BODY_B X Y");
                m.revolutes.push_back({t[1], m.body_index(t[2]), m.body_index(t[3]), point_of(t, 4, w)});
            } else if (k == "prismatic") {  // prismatic NAME BODY_A BODY_B X Y UX UY
                if (t.size() != 8) throw std::runtime_error("usage: prismatic NAME BODY_A BODY_B X Y UX UY");
                Prismatic p{t[1], m.body_index(t[2]), m.body_index(t[3]), point_of(t, 4, w), point_of(t, 6, w)};
                const double n = std::hypot(p.u[0], p.u[1]);
                if (n > 0) p.u = {p.u[0] / n, p.u[1] / n};
                m.prismatics.push_back(p);
            } else if (k == "actuator") {  // actuator NAME BODY_A XA YA BODY_B XB YB [key=value ...]
                if (t.size() < 8) throw std::runtime_error("usage: actuator NAME BODY_A XA YA BODY_B XB YB [bore= rod= min= max= mass=]");
                Actuator a;
                a.name = t[1];
                a.body_a = m.body_index(t[2]);
                a.pa = point_of(t, 3, w);
                a.body_b = m.body_index(t[5]);
                a.pb = point_of(t, 6, w);
                a.l_ref = std::hypot(a.pb[0] - a.pa[0], a.pb[1] - a.pa[1]);
                auto kv = keyvals(t, 8, w);
                if (!kv.empty()) {
                    a.has_cylinder = true;
                    for (const auto& [key, val] : kv) {
                        const double x = num(val, w);
                        if (key == "bore") a.cyl.bore = x;
                        else if (key == "rod") a.cyl.rod = x;
                        else if (key == "min") a.cyl.lmin = x;
                        else if (key == "max") a.cyl.lmax = x;
                        else if (key == "mass") a.cyl.mass = x;
                        else throw std::runtime_error("unknown actuator key '" + key + "'");
                    }
                }
                m.actuators.push_back(a);
            } else if (k == "mass") {  // mass NAME BODY X Y M [J] [group=..]
                if (t.size() < 6) throw std::runtime_error("usage: mass NAME BODY X Y MASS [J] [group=G]");
                MassItem mi;
                mi.name = t[1];
                mi.body = m.body_index(t[2]);
                mi.c = point_of(t, 3, w);
                mi.m = num(t[5], w);
                std::size_t i = 6;
                if (t.size() > 6 && t[6].find('=') == std::string::npos) { mi.j = num(t[6], w); i = 7; }
                auto kv = keyvals(t, i, w);
                mi.group = kv.count("group") ? kv["group"] : m.bodies[mi.body];
                m.masses.push_back(mi);
            } else if (k == "force") {  // force NAME BODY X Y FX FY (N, global directions)
                if (t.size() != 7) throw std::runtime_error("usage: force NAME BODY X Y FX FY");
                ExternalForce ef;
                ef.name = t[1];
                ef.body = m.body_index(t[2]);
                ef.p = point_of(t, 3, w);
                ef.f = {num(t[5], w), num(t[6], w)};
                m.forces.push_back(ef);
            } else if (k == "outline") {  // outline BODY X1 Y1 X2 Y2 ... (drawing only)
                if (t.size() < 4 || t.size() % 2 != 0) throw std::runtime_error("usage: outline BODY X1 Y1 [X2 Y2 ...]");
                const int b = m.body_index(t[1]);
                for (std::size_t i = 2; i + 1 < t.size(); i += 2) m.outlines[b].push_back(point_of(t, i, w));
            } else if (k == "label") {  // label BODY TEXT (drawing only)
                if (t.size() < 3) throw std::runtime_error("usage: label BODY TEXT");
                std::string s;
                for (std::size_t i = 2; i < t.size(); ++i) s += (i > 2 ? " " : "") + t[i];
                m.labels[m.body_index(t[1])] = s;
            } else if (k == "masses") {
                std::string p = t.at(1);
                if (!(p.size() > 0 && (p[0] == '/' || (p.size() > 1 && p[1] == ':')))) p = base_dir + "/" + p;
                m.load_masses_csv(p);
            } else if (k == "hydraulics") {
                auto kv = keyvals(t, 1, w);
                for (const auto& [key, val] : kv) {
                    if (key == "beta") m.hydraulics.beta = num(val, w);
                    else if (key == "dead_volume") m.hydraulics.dead_volume = num(val, w);
                    else throw std::runtime_error("unknown hydraulics key '" + key + "'");
                }
            } else {
                throw std::runtime_error("unknown keyword '" + k + "'");
            }
        } catch (const std::exception& e) {
            throw std::runtime_error("model " + w + ": " + e.what());
        }
    }
    m.validate();
    return m;
}

// ---------------------------------------------------------------- kinematics

Kinematics::Kinematics(const Model& m) : m_(m) { m_.validate(); }

Pose Kinematics::reference() const {
    Pose s;
    s.q.assign(m_.n_coords(), 0.0);
    for (const auto& a : m_.actuators) s.lengths.push_back(a.l_ref);
    return s;
}

Vec2 Kinematics::point(const Pose& s, int body, const Vec2& p) const {
    if (body == 0) return p;
    const double* q = &s.q[3 * (body - 1)];
    Vec2 r = rot(q[2], p);
    return {q[0] + r[0], q[1] + r[1]};
}

double Kinematics::angle(const Pose& s, int body) const { return body == 0 ? 0.0 : s.q[3 * (body - 1) + 2]; }

Vec2 Kinematics::actuator_dir(const Pose& s, int k) const {
    const auto& a = m_.actuators[k];
    Vec2 ra = point(s, a.body_a, a.pa), rb = point(s, a.body_b, a.pb);
    const double L = std::hypot(rb[0] - ra[0], rb[1] - ra[1]);
    return {(rb[0] - ra[0]) / L, (rb[1] - ra[1]) / L};
}

std::vector<double> Kinematics::residual(const Pose& s) const {
    std::vector<double> r;
    r.reserve(m_.n_constraints());
    for (const auto& j : m_.revolutes) {
        Vec2 a = point(s, j.body_a, j.p), b = point(s, j.body_b, j.p);
        r.push_back(a[0] - b[0]);
        r.push_back(a[1] - b[1]);
    }
    for (const auto& j : m_.prismatics) {
        const double ta = angle(s, j.body_a), tb = angle(s, j.body_b);
        Vec2 n = rot(ta, {-j.u[1], j.u[0]});
        Vec2 a = point(s, j.body_a, j.p), b = point(s, j.body_b, j.p);
        r.push_back(ta - tb);
        r.push_back(n[0] * (b[0] - a[0]) + n[1] * (b[1] - a[1]));
    }
    for (std::size_t k = 0; k < m_.actuators.size(); ++k) {
        const auto& a = m_.actuators[k];
        Vec2 ra = point(s, a.body_a, a.pa), rb = point(s, a.body_b, a.pb);
        r.push_back(std::hypot(rb[0] - ra[0], rb[1] - ra[1]) - s.lengths[k]);
    }
    return r;
}

Matrix Kinematics::jacobian(const Pose& s) const {
    const int n = m_.n_coords();
    Matrix J(m_.n_constraints(), n);
    // d(point)/dq of body b: columns 3(b-1)+{0,1,2}
    auto add_point = [&](std::size_t row, int body, const Vec2& p, const Vec2& w, double sign) {
        if (body == 0) return;  // ground does not move
        const int c = 3 * (body - 1);
        Vec2 d = drot(angle(s, body), p);
        J(row, c + 0) += sign * w[0];
        J(row, c + 1) += sign * w[1];
        J(row, c + 2) += sign * (w[0] * d[0] + w[1] * d[1]);
    };
    std::size_t row = 0;
    for (const auto& j : m_.revolutes) {
        add_point(row, j.body_a, j.p, {1, 0}, +1);
        add_point(row, j.body_b, j.p, {1, 0}, -1);
        ++row;
        add_point(row, j.body_a, j.p, {0, 1}, +1);
        add_point(row, j.body_b, j.p, {0, 1}, -1);
        ++row;
    }
    for (const auto& j : m_.prismatics) {
        if (j.body_a) J(row, 3 * (j.body_a - 1) + 2) += 1;
        if (j.body_b) J(row, 3 * (j.body_b - 1) + 2) -= 1;
        ++row;
        const double ta = angle(s, j.body_a);
        Vec2 n0{-j.u[1], j.u[0]};
        Vec2 nv = rot(ta, n0), dn = drot(ta, n0);
        Vec2 a = point(s, j.body_a, j.p), b = point(s, j.body_b, j.p);
        add_point(row, j.body_b, j.p, nv, +1);
        add_point(row, j.body_a, j.p, nv, -1);
        if (j.body_a) J(row, 3 * (j.body_a - 1) + 2) += dn[0] * (b[0] - a[0]) + dn[1] * (b[1] - a[1]);
        ++row;
    }
    for (std::size_t k = 0; k < m_.actuators.size(); ++k) {
        const auto& a = m_.actuators[k];
        Vec2 u = actuator_dir(s, static_cast<int>(k));
        add_point(row, a.body_b, a.pb, u, +1);
        add_point(row, a.body_a, a.pa, u, -1);
        ++row;
    }
    return J;
}

bool Kinematics::newton(Pose& s, int max_iter, double tol) const {
    for (int it = 0; it < max_iter; ++it) {
        auto r = residual(s);
        if (norm_inf(r) < tol) return true;
        LU lu(jacobian(s));
        if (lu.singular()) return false;
        auto dx = lu.solve(r);
        double step = norm_inf(dx);
        for (std::size_t i = 0; i < s.q.size(); ++i) s.q[i] -= dx[i];
        if (!std::isfinite(step) || step > 1e6) return false;
    }
    return norm_inf(residual(s)) < 1e-8;
}

Pose Kinematics::solve(const Pose& start, const std::vector<double>& lengths, double max_step_mm) const {
    if (lengths.size() != m_.actuators.size()) throw std::invalid_argument("solve: wrong number of actuator lengths");
    Pose s = start;
    if (s.q.size() != static_cast<std::size_t>(m_.n_coords())) s = reference();
    // continuation in steps of at most max_step_mm; a step that fails is halved (up to 1/1024)
    double t = 0.0;
    const auto l0 = s.lengths;
    double dmax = 0.0;
    for (std::size_t k = 0; k < lengths.size(); ++k) dmax = std::max(dmax, std::fabs(lengths[k] - l0[k]));
    if (dmax == 0.0 && newton(s)) return s;
    double dt = dmax > 0 ? std::min(1.0, max_step_mm / dmax) : 1.0;
    while (t < 1.0) {
        const double t1 = std::min(1.0, t + dt);
        Pose trial = s;
        for (std::size_t k = 0; k < lengths.size(); ++k) trial.lengths[k] = l0[k] + (lengths[k] - l0[k]) * t1;
        if (newton(trial)) {
            s = trial;
            t = t1;
        } else {
            dt *= 0.5;
            if (dt * dmax < max_step_mm / 1024.0 || dt < 1e-9) {
                std::ostringstream os;
                os << "position solver failed (unreachable or singular position) at lengths:";
                for (std::size_t k = 0; k < lengths.size(); ++k) os << ' ' << m_.actuators[k].name << '=' << trial.lengths[k];
                throw std::runtime_error(os.str());
            }
        }
    }
    return s;
}

Matrix Kinematics::sensitivity(const Pose& s) const {
    LU lu(jacobian(s));
    if (lu.singular()) throw std::runtime_error("sensitivity: singular Jacobian (dead point of the mechanism)");
    const int n = m_.n_coords(), na = static_cast<int>(m_.actuators.size());
    const int row0 = m_.n_constraints() - na;
    Matrix dq(n, na);
    for (int k = 0; k < na; ++k) {
        std::vector<double> e(m_.n_constraints(), 0.0);
        e[row0 + k] = 1.0;  // dPhi/dL_k = -1  ->  Phi_q dq = +e_k
        auto x = lu.solve(e);
        for (int i = 0; i < n; ++i) dq(i, k) = x[i];
    }
    return dq;
}

Vec2 Kinematics::point_rate(const Pose& s, const Matrix& dq, int k, int body, const Vec2& p) const {
    if (body == 0) return {0.0, 0.0};
    const int c = 3 * (body - 1);
    Vec2 d = drot(angle(s, body), p);
    const double th = dq(c + 2, k);
    return {dq(c, k) + d[0] * th, dq(c + 1, k) + d[1] * th};
}

std::map<std::string, double> parse_assignments(const std::vector<std::string>& args) {
    std::map<std::string, double> out;
    for (const auto& a : args) {
        auto p = a.find('=');
        if (p == std::string::npos) throw std::runtime_error("expected NAME=VALUE, got '" + a + "'");
        out[a.substr(0, p)] = num(a.substr(p + 1), a);
    }
    return out;
}

}  // namespace cylload
