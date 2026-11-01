// SPDX-License-Identifier: MIT
// Minimal STEP Part 21 reader for assembly trees and bounding boxes. It follows
// NEXT_ASSEMBLY_USAGE_OCCURRENCE links, composes the placements given by
// ITEM_DEFINED_TRANSFORMATION, and collects CARTESIAN_POINTs reachable from each
// part's shape representation (vertices and curve points of the B-rep).
#include "cylload/step.hpp"
#include "step_entities.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace cylload {

namespace {

}  // namespace

namespace step_detail {
Mat4 identity() { return {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}; }

Mat4 mul(const Mat4& a, const Mat4& b) {
    Mat4 c{};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            double s = 0;
            for (int k = 0; k < 4; ++k) s += a[4 * i + k] * b[4 * k + j];
            c[4 * i + j] = s;
        }
    return c;
}

Mat4 rigid_inverse(const Mat4& m) {
    Mat4 r = identity();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) r[4 * i + j] = m[4 * j + i];
    for (int i = 0; i < 3; ++i) r[4 * i + 3] = -(r[4 * i] * m[3] + r[4 * i + 1] * m[7] + r[4 * i + 2] * m[11]);
    return r;
}

std::vector<long> refs(const std::string& s) {
    std::vector<long> out;
    bool str = false;
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '\'') { str = !str; continue; }
        if (str || s[i] != '#') continue;
        std::size_t j = i + 1;
        long v = 0;
        bool any = false;
        while (j < s.size() && std::isdigit(static_cast<unsigned char>(s[j]))) { v = v * 10 + (s[j] - '0'); ++j; any = true; }
        if (any) out.push_back(v);
        i = j - 1;
    }
    return out;
}

std::vector<double> numbers(const std::string& s) {
    std::vector<double> out;
    bool str = false;
    for (std::size_t i = 0; i < s.size(); ++i) {
        const char c = s[i];
        if (c == '\'') { str = !str; continue; }
        if (str) continue;
        if (c == '#') {  // skip reference
            while (i + 1 < s.size() && std::isdigit(static_cast<unsigned char>(s[i + 1]))) ++i;
            continue;
        }
        const bool start = std::isdigit(static_cast<unsigned char>(c)) ||
                           ((c == '-' || c == '+' || c == '.') && i + 1 < s.size() &&
                            (std::isdigit(static_cast<unsigned char>(s[i + 1])) || s[i + 1] == '.'));
        if (!start) continue;
        std::string t;
        std::size_t j = i;
        while (j < s.size() && (std::isdigit(static_cast<unsigned char>(s[j])) || s[j] == '.' || s[j] == '-' || s[j] == '+' ||
                                s[j] == 'E' || s[j] == 'e' || s[j] == 'D' || s[j] == 'd')) {
            t += (s[j] == 'D' || s[j] == 'd') ? 'E' : s[j];
            ++j;
        }
        char* end = nullptr;
        const double v = std::strtod(t.c_str(), &end);
        if (end != t.c_str()) out.push_back(v);
        i = j - 1;
    }
    return out;
}

std::vector<std::string> strings(const std::string& s) {
    std::vector<std::string> out;
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (s[i] != '\'') continue;
        std::string t;
        std::size_t j = i + 1;
        for (; j < s.size(); ++j) {
            if (s[j] == '\'') {
                if (j + 1 < s.size() && s[j + 1] == '\'') { t += '\''; ++j; continue; }
                break;
            }
            t += s[j];
        }
        out.push_back(t);
        i = j;
    }
    return out;
}

void utf8(std::string& o, unsigned cp) {
    if (cp < 0x80) o += static_cast<char>(cp);
    else if (cp < 0x800) { o += static_cast<char>(0xC0 | (cp >> 6)); o += static_cast<char>(0x80 | (cp & 0x3F)); }
    else { o += static_cast<char>(0xE0 | (cp >> 12)); o += static_cast<char>(0x80 | ((cp >> 6) & 0x3F)); o += static_cast<char>(0x80 | (cp & 0x3F)); }
}

Entities parse(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open STEP file '" + path + "'");
    std::stringstream ss;
    ss << f.rdbuf();
    const std::string txt = ss.str();
    auto p = txt.find("DATA;");
    if (p == std::string::npos) throw std::runtime_error("not a STEP Part 21 file (no DATA section): " + path);
    p += 5;
    Entities ent;
    const std::size_t n = txt.size();
    while (p < n) {
        // find start of next statement
        while (p < n && txt[p] != '#') {
            if (txt.compare(p, 7, "ENDSEC;") == 0) return ent;
            ++p;
        }
        if (p >= n) break;
        // statement end: ';' outside strings
        std::size_t q = p;
        bool str = false;
        for (; q < n; ++q) {
            if (txt[q] == '\'') str = !str;
            else if (!str && txt[q] == ';') break;
        }
        const std::string st = txt.substr(p, q - p);
        p = q + 1;
        const auto eq = st.find('=');
        if (eq == std::string::npos) continue;
        const long id = std::atol(st.c_str() + 1);
        std::size_t a = st.find_first_not_of(" \t\r\n", eq + 1);
        if (a == std::string::npos) continue;
        Entity e;
        const auto open = st.find('(', a);
        const auto close = st.rfind(')');
        if (open == std::string::npos || close == std::string::npos || close < open) continue;
        if (st[a] != '(') {
            e.type = st.substr(a, open - a);
            e.type.erase(e.type.find_last_not_of(" \t\r\n") + 1);
        }
        e.args = st.substr(open + 1, close - open - 1);
        ent.emplace(id, std::move(e));
    }
    return ent;
}

const Entity* get(const Entities& e, long id) {
    auto it = e.find(id);
    return it == e.end() ? nullptr : &it->second;
}

std::array<double, 3> vec3(const Entities& e, long id, std::array<double, 3> dflt) {
    const Entity* x = get(e, id);
    if (!x) return dflt;
    auto v = numbers(x->args);
    std::array<double, 3> r{0, 0, 0};
    for (std::size_t i = 0; i < 3 && i < v.size(); ++i) r[i] = v[i];
    return r;
}

Mat4 placement(const Entities& e, long id) {
    const Entity* ax = get(e, id);
    if (!ax) return identity();
    auto r = refs(ax->args);
    std::array<double, 3> P = r.size() > 0 ? vec3(e, r[0], {0, 0, 0}) : std::array<double, 3>{0, 0, 0};
    std::array<double, 3> Z = r.size() > 1 ? vec3(e, r[1], {0, 0, 1}) : std::array<double, 3>{0, 0, 1};
    std::array<double, 3> X = r.size() > 2 ? vec3(e, r[2], {1, 0, 0}) : std::array<double, 3>{1, 0, 0};
    auto nrm = [](std::array<double, 3>& v, std::array<double, 3> d) {
        const double n = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        if (n > 0) for (auto& c : v) c /= n; else v = d;
    };
    nrm(Z, {0, 0, 1});
    const double dz = X[0] * Z[0] + X[1] * Z[1] + X[2] * Z[2];
    for (int i = 0; i < 3; ++i) X[i] -= dz * Z[i];
    nrm(X, {1, 0, 0});
    std::array<double, 3> Y{Z[1] * X[2] - Z[2] * X[1], Z[2] * X[0] - Z[0] * X[2], Z[0] * X[1] - Z[1] * X[0]};
    return {X[0], Y[0], Z[0], P[0], X[1], Y[1], Z[1], P[1], X[2], Y[2], Z[2], P[2], 0, 0, 0, 1};
}

}  // namespace step_detail

using namespace step_detail;

std::string decode_step_string(const std::string& s) {
    std::string o;
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (s.compare(i, 4, "\\X2\\") == 0) {
            auto end = s.find("\\X0\\", i + 4);
            if (end == std::string::npos) break;
            for (std::size_t k = i + 4; k + 4 <= end; k += 4) utf8(o, static_cast<unsigned>(std::stoul(s.substr(k, 4), nullptr, 16)));
            i = end + 3;
        } else if (s.compare(i, 3, "\\X\\") == 0 && i + 5 <= s.size()) {
            utf8(o, static_cast<unsigned>(std::stoul(s.substr(i + 3, 2), nullptr, 16)));
            i += 4;
        } else if (s.compare(i, 2, "''") == 0) {
            o += '\'';
            ++i;
        } else {
            o += s[i];
        }
    }
    return o;
}

StepModel read_step(const std::string& path, const std::string& root_name, bool mass_properties, bool geometry) {
    const Entities ent = parse(path);
    StepModel out;
    out.n_entities = ent.size();
    auto is = [&](long id, const char* t) { const Entity* x = get(ent, id); return x && x->type == t; };

    std::unordered_map<long, std::string> product_name;
    for (const auto& [id, e] : ent)
        if (e.type == "PRODUCT") {
            auto s = strings(e.args);
            product_name[id] = s.empty() ? std::string() : decode_step_string(s[0]);
        }
    std::unordered_map<long, std::string> pd_name;
    for (const auto& [id, e] : ent) {
        if (e.type != "PRODUCT_DEFINITION") continue;
        for (long r : refs(e.args)) {
            const Entity* x = get(ent, r);
            if (!x || x->type.rfind("PRODUCT_DEFINITION_FORMATION", 0) != 0) continue;
            auto rr = refs(x->args);
            if (!rr.empty() && product_name.count(rr[0])) pd_name[id] = product_name[rr[0]];
        }
    }
    std::map<long, std::pair<long, long>> nauo;  // id -> (parent pd, child pd)
    for (const auto& [id, e] : ent)
        if (e.type == "NEXT_ASSEMBLY_USAGE_OCCURRENCE") {
            auto r = refs(e.args);
            if (r.size() >= 2) nauo[id] = {r[0], r[1]};
        }
    std::unordered_map<long, long> pds_nauo, pds_pd;
    for (const auto& [id, e] : ent)
        if (e.type == "PRODUCT_DEFINITION_SHAPE") {
            auto r = refs(e.args);
            if (r.empty()) continue;
            if (nauo.count(r[0])) pds_nauo[id] = r[0];
            else if (is(r[0], "PRODUCT_DEFINITION")) pds_pd[id] = r[0];
        }
    // placement of each NAUO: parent <- child
    std::unordered_map<long, Mat4> tr;
    for (const auto& [id, e] : ent) {
        if (e.type != "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION") continue;
        auto r = refs(e.args);
        if (r.size() < 2 || !pds_nauo.count(r[1])) continue;
        const Entity* rel = get(ent, r[0]);
        if (!rel) continue;
        auto k = rel->args.find("REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION");
        if (k == std::string::npos) continue;
        auto o = rel->args.find('(', k), c = rel->args.find(')', k);
        if (o == std::string::npos || c == std::string::npos) continue;
        auto idt = refs(rel->args.substr(o, c - o));
        if (idt.empty()) continue;
        const Entity* t = get(ent, idt[0]);
        if (!t) continue;
        auto ra = refs(t->args);
        if (ra.size() < 2) continue;
        tr[pds_nauo[r[1]]] = mul(placement(ent, ra[0]), rigid_inverse(placement(ent, ra[1])));
    }
    // shape representations per product definition
    std::unordered_map<long, std::vector<long>> srr;
    for (const auto& [id, e] : ent)
        if (e.type == "SHAPE_REPRESENTATION_RELATIONSHIP") {
            auto r = refs(e.args);
            if (r.size() >= 2) { srr[r[0]].push_back(r[1]); srr[r[1]].push_back(r[0]); }
        }
    std::unordered_map<long, std::vector<long>> reps;
    for (const auto& [id, e] : ent) {
        if (e.type != "SHAPE_DEFINITION_REPRESENTATION") continue;
        auto r = refs(e.args);
        if (r.size() < 2 || !pds_pd.count(r[0])) continue;
        const long pd = pds_pd[r[0]];
        reps[pd].push_back(r[1]);
        for (long x : srr[r[1]])
            if (is(x, "ADVANCED_BREP_SHAPE_REPRESENTATION")) reps[pd].push_back(x);
    }
    // tree
    std::map<long, std::vector<std::pair<long, long>>> children;  // parent pd -> (nauo, child pd)
    std::set<long> is_child;
    for (const auto& [id, pc] : nauo) {
        children[pc.first].push_back({id, pc.second});
        is_child.insert(pc.second);
    }
    long root = 0;
    for (const auto& [pd, ch] : children) {
        const std::string nm = pd_name.count(pd) ? pd_name[pd] : std::string();
        if (!root_name.empty() ? nm == root_name : !is_child.count(pd)) { root = pd; break; }
    }
    if (!root) {
        std::string avail;
        for (const auto& [pd, ch] : children) if (!is_child.count(pd)) avail += " '" + pd_name[pd] + "'";
        throw std::runtime_error("STEP root assembly '" + root_name + "' not found; top assemblies:" + avail);
    }
    out.root = pd_name[root];

    std::unordered_map<long, std::vector<std::array<double, 3>>> pts_cache;
    auto points = [&](long pd) -> const std::vector<std::array<double, 3>>& {
        auto it = pts_cache.find(pd);
        if (it != pts_cache.end()) return it->second;
        std::vector<std::array<double, 3>> res;
        std::vector<long> stack(reps[pd].begin(), reps[pd].end());
        std::set<long> seen;
        while (!stack.empty()) {
            const long x = stack.back();
            stack.pop_back();
            if (!seen.insert(x).second) continue;
            const Entity* e = get(ent, x);
            if (!e || e->type.empty()) continue;  // complex entities are not followed
            if (e->type == "CARTESIAN_POINT") {
                auto v = numbers(e->args);
                if (v.size() == 3) res.push_back({v[0], v[1], v[2]});
                continue;
            }
            if (e->type == "DIRECTION" || e->type == "AXIS2_PLACEMENT_3D") continue;
            for (long r : refs(e->args)) stack.push_back(r);
        }
        return pts_cache.emplace(pd, std::move(res)).first->second;
    };

    std::unordered_map<long, VolumeIntegrals> vol_cache;
    auto volumes = [&](long pd) -> const VolumeIntegrals& {
        auto it = vol_cache.find(pd);
        if (it != vol_cache.end()) return it->second;
        VolumeIntegrals vi;
        try {
            vi = integrate_solids(ent, reps[pd]);
        } catch (const std::exception& ex) {
            throw std::runtime_error("mass properties of '" + pd_name[pd] + "': " + ex.what());
        }
        return vol_cache.emplace(pd, vi).first->second;
    };

    std::unordered_map<long, std::vector<SolidGeometry>> geo_cache;
    struct Item { long pd; Mat4 M; std::vector<std::string> path; };
    std::vector<Item> stack{{root, identity(), {out.root}}};
    while (!stack.empty()) {
        Item it = stack.back();
        stack.pop_back();
        for (const auto& [n, c] : children[it.pd]) {
            Mat4 M = mul(it.M, tr.count(n) ? tr[n] : identity());
            auto path2 = it.path;
            path2.push_back(pd_name.count(c) ? pd_name[c] : std::string("?"));
            ++out.n_tree_nodes;
            const auto& P = points(c);
            if (!P.empty()) {
                StepOccurrence o;
                o.path = path2;
                o.product_definition = c;
                o.n_points = P.size();
                o.lo = {1e300, 1e300, 1e300};
                o.hi = {-1e300, -1e300, -1e300};
                for (const auto& p : P) {
                    for (int i = 0; i < 3; ++i) {
                        const double v = M[4 * i] * p[0] + M[4 * i + 1] * p[1] + M[4 * i + 2] * p[2] + M[4 * i + 3];
                        o.lo[i] = std::min(o.lo[i], v);
                        o.hi[i] = std::max(o.hi[i], v);
                    }
                }
                if (mass_properties) {
                    const auto& vi = volumes(c);
                    o.faces = vi.faces;
                    o.faces_triangulated = vi.faces_triangulated;
                    if (vi.v[0] > 0) {
                        const double V = vi.v[0];
                        const std::array<double, 3> d{vi.v[1] / V, vi.v[2] / V, vi.v[3] / V};  // centroid - ref (local)
                        const std::array<double, 3> cl{vi.ref[0] + d[0], vi.ref[1] + d[1], vi.ref[2] + d[2]};
                        // local covariance about the centroid
                        const double Sxx = vi.v[4], Syy = vi.v[5], Szz = vi.v[6], Sxy = vi.v[7], Syz = vi.v[8], Szx = vi.v[9];
                        const double C[9] = {Sxx - V * d[0] * d[0], Sxy - V * d[0] * d[1], Szx - V * d[0] * d[2],
                                             Sxy - V * d[0] * d[1], Syy - V * d[1] * d[1], Syz - V * d[1] * d[2],
                                             Szx - V * d[0] * d[2], Syz - V * d[1] * d[2], Szz - V * d[2] * d[2]};
                        o.volume = V;
                        o.volume_check = std::max(std::fabs(vi.v[10] - V), std::fabs(vi.v[11] - V)) / V;
                        const double trace = vi.v[4] + vi.v[5] + vi.v[6];
                        o.moment_check = trace > 0 ? std::max({std::fabs(vi.v[12] - vi.v[4]), std::fabs(vi.v[13] - vi.v[5]),
                                                               std::fabs(vi.v[14] - vi.v[6])}) / trace : 0.0;
                        for (int i = 0; i < 3; ++i)
                            o.centroid[i] = M[4 * i] * cl[0] + M[4 * i + 1] * cl[1] + M[4 * i + 2] * cl[2] + M[4 * i + 3];
                        for (int i = 0; i < 3; ++i)
                            for (int j = 0; j < 3; ++j) {
                                double sum = 0;
                                for (int k = 0; k < 3; ++k)
                                    for (int l = 0; l < 3; ++l) sum += M[4 * i + k] * C[3 * k + l] * M[4 * j + l];
                                o.cov[3 * i + j] = sum;
                            }
                    }
                }
                if (geometry) {
                    auto it2 = geo_cache.find(c);
                    if (it2 == geo_cache.end()) it2 = geo_cache.emplace(c, solid_geometry(ent, reps[c])).first;
                    auto tf = [&](const std::array<double, 3>& p) {
                        std::array<double, 3> q{};
                        for (int i = 0; i < 3; ++i) q[i] = M[4 * i] * p[0] + M[4 * i + 1] * p[1] + M[4 * i + 2] * p[2] + M[4 * i + 3];
                        return q;
                    };
                    for (const auto& sgm : it2->second) {
                        StepSolid so;
                        for (const auto& f : sgm.faces) {
                            std::vector<std::vector<std::array<double, 3>>> loops;
                            for (const auto& lp : f) {
                                std::vector<std::array<double, 3>> l2;
                                for (const auto& p : lp) l2.push_back(tf(p));
                                loops.push_back(std::move(l2));
                            }
                            so.faces.push_back(std::move(loops));
                        }
                        so.curved = sgm.curved;
                        for (const auto& ed : sgm.edges) {
                            std::vector<std::array<double, 3>> e2;
                            for (const auto& p : ed) e2.push_back(tf(p));
                            so.edges.push_back(std::move(e2));
                        }
                        o.solids.push_back(std::move(so));
                    }
                }
                out.occurrences.push_back(o);
            }
            if (children.count(c)) stack.push_back({c, M, path2});
        }
    }
    return out;
}

void MassProps::add(const MassProps& o) {
    if (o.volume <= 0) return;
    const double V = volume + o.volume;
    std::array<double, 3> c;
    for (int i = 0; i < 3; ++i) c[i] = (volume * centroid[i] + o.volume * o.centroid[i]) / V;
    std::array<double, 9> C{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            const double a = centroid[i] - c[i], b = centroid[j] - c[j];
            const double oa = o.centroid[i] - c[i], ob = o.centroid[j] - c[j];
            C[3 * i + j] = cov[3 * i + j] + volume * a * b + o.cov[3 * i + j] + o.volume * oa * ob;
        }
    volume = V;
    centroid = c;
    cov = C;
}

std::vector<StepGroup> group_by_level(const StepModel& m, std::size_t level) {
    std::vector<StepGroup> out;
    std::map<std::string, std::size_t> idx;
    for (const auto& o : m.occurrences) {
        const std::string& k = o.path[std::min(level, o.path.size() - 1)];
        auto it = idx.find(k);
        if (it == idx.end()) {
            idx[k] = out.size();
            StepGroup g;
            g.key = k;
            g.lo = o.lo;
            g.hi = o.hi;
            out.push_back(g);
            it = idx.find(k);
        }
        auto& g = out[it->second];
        ++g.n_occurrences;
        g.n_points += o.n_points;
        g.faces += o.faces;
        g.faces_triangulated += o.faces_triangulated;
        g.volume_check = std::max(g.volume_check, o.volume_check);
        g.moment_check = std::max(g.moment_check, o.moment_check);
        MassProps mp;
        mp.volume = o.volume;
        mp.centroid = o.centroid;
        mp.cov = o.cov;
        g.mass.add(mp);
        for (int i = 0; i < 3; ++i) {
            g.lo[i] = std::min(g.lo[i], o.lo[i]);
            g.hi[i] = std::max(g.hi[i], o.hi[i]);
        }
    }
    return out;
}

StepGroup union_by_prefix(const std::vector<StepGroup>& groups, const std::string& prefix) {
    StepGroup u;
    u.key = prefix;
    bool any = false;
    for (const auto& g : groups) {
        if (g.key.rfind(prefix, 0) != 0) continue;
        if (!any) { u.lo = g.lo; u.hi = g.hi; any = true; }
        for (int i = 0; i < 3; ++i) {
            u.lo[i] = std::min(u.lo[i], g.lo[i]);
            u.hi[i] = std::max(u.hi[i], g.hi[i]);
        }
        u.n_occurrences += g.n_occurrences;
        u.n_points += g.n_points;
        u.faces += g.faces;
        u.faces_triangulated += g.faces_triangulated;
        u.volume_check = std::max(u.volume_check, g.volume_check);
        u.moment_check = std::max(u.moment_check, g.moment_check);
        u.mass.add(g.mass);
    }
    if (!any) throw std::runtime_error("no STEP group starts with '" + prefix + "'");
    return u;
}

}  // namespace cylload
