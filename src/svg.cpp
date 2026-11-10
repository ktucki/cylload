// SPDX-License-Identifier: MIT
#include "cylload/svg.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace cylload {

namespace {

constexpr double kPt = 0.352778;  // mm per point

std::string fmt(double v, int prec = 3) {
    std::ostringstream os;
    os.setf(std::ios::fixed);
    os.precision(prec);
    os << v;
    return os.str();
}

std::string esc_plain(const std::string& s);

/// escape for SVG text and turn _{..} / _x into subscripts and ^{..} / ^x into superscripts
std::string esc(const std::string& s) {
    std::string o;
    for (std::size_t i = 0; i < s.size(); ++i) {
        const char c = s[i];
        if ((c == '_' || c == '^') && i + 1 < s.size()) {
            std::string part;
            if (s[i + 1] == '{') {
                const auto e = s.find('}', i + 2);
                if (e == std::string::npos) { o += esc_plain(std::string(1, c)); continue; }
                part = s.substr(i + 2, e - i - 2);
                i = e;
            } else {
                part = s.substr(i + 1, 1);
                ++i;
            }
            o += std::string("<tspan font-size=\"70%\" baseline-shift=\"") + (c == '_' ? "sub" : "super") + "\">" + esc_plain(part) + "</tspan>";
            continue;
        }
        o += esc_plain(std::string(1, c));
    }
    return o;
}

std::string esc_plain(const std::string& s) {
    std::string o;
    for (char c : s) {
        if (c == '&') o += "&amp;";
        else if (c == '<') o += "&lt;";
        else if (c == '>') o += "&gt;";
        else if (c == '"') o += "&quot;";
        else o += c;
    }
    return o;
}

std::string dash(LineStyle s, double w) {
    const double k = std::max(w, 0.15) / 0.2;
    switch (s) {
        case LineStyle::Dashed: return " stroke-dasharray=\"" + fmt(1.6 * k) + "," + fmt(0.8 * k) + "\"";
        case LineStyle::Dotted: return " stroke-dasharray=\"" + fmt(0.2 * k) + "," + fmt(0.6 * k) + "\" stroke-linecap=\"round\"";
        case LineStyle::DashDot:
            return " stroke-dasharray=\"" + fmt(1.6 * k) + "," + fmt(0.6 * k) + "," + fmt(0.2 * k) + "," + fmt(0.6 * k) + "\"";
        default: return "";
    }
}

std::string marker_svg(Marker m, double x, double y, double r = 0.7, const std::string& col = "#000") {
    const std::string open = " fill=\"#fff\" stroke=\"" + col + "\" stroke-width=\"0.2\"";
    const std::string full = " fill=\"" + col + "\" stroke=\"" + col + "\" stroke-width=\"0.18\"";
    auto poly = [&](const std::vector<std::pair<double, double>>& p, bool filled) {
        std::string s = "<polygon points=\"";
        for (const auto& q : p) s += fmt(x + q.first) + "," + fmt(y + q.second) + " ";
        return s + "\"" + (filled ? full : open) + "/>";
    };
    switch (m) {
        case Marker::Circle: return "<circle cx=\"" + fmt(x) + "\" cy=\"" + fmt(y) + "\" r=\"" + fmt(r * 0.85) + "\"" + full + "/>";
        case Marker::CircleOpen: return "<circle cx=\"" + fmt(x) + "\" cy=\"" + fmt(y) + "\" r=\"" + fmt(r * 0.85) + "\"" + open + "/>";
        case Marker::Square:
        case Marker::SquareOpen:
            return poly({{-r * 0.8, -r * 0.8}, {r * 0.8, -r * 0.8}, {r * 0.8, r * 0.8}, {-r * 0.8, r * 0.8}}, m == Marker::Square);
        case Marker::Triangle:
        case Marker::TriangleOpen: return poly({{0, -r}, {r * 0.9, r * 0.7}, {-r * 0.9, r * 0.7}}, m == Marker::Triangle);
        case Marker::Diamond:
        case Marker::DiamondOpen: return poly({{0, -r}, {r, 0}, {0, r}, {-r, 0}}, m == Marker::Diamond);
        case Marker::Cross:
            return "<path d=\"M" + fmt(x - r) + " " + fmt(y - r) + "L" + fmt(x + r) + " " + fmt(y + r) + "M" + fmt(x - r) + " " + fmt(y + r) +
                   "L" + fmt(x + r) + " " + fmt(y - r) + "\" stroke=\"" + col + "\" stroke-width=\"0.22\"/>";
        default: return "";
    }
}

int decimals_for(double step) {
    int d = 0;
    while (d < 6 && std::fabs(step * std::pow(10.0, d) - std::round(step * std::pow(10.0, d))) > 1e-9 * std::pow(10.0, d)) ++d;
    return d;
}

std::string tick_label(double v, int dec) {
    if (std::fabs(v) < 1e-12) v = 0.0;
    std::string s = fmt(v, dec);
    if (s.size() > 1 && s[0] == '-') s = "\xE2\x88\x92" + s.substr(1);  // typographic minus
    return s;
}

struct Ticks {
    double lo = 0, hi = 1;
    std::vector<double> at;
    std::vector<std::string> label;
    bool log = false;
};

Ticks make_ticks(const Axis& a, double dlo, double dhi) {
    Ticks t;
    t.log = a.log;
    if (a.log) {
        double lo = a.fixed ? a.min : dlo, hi = a.fixed ? a.max : dhi;
        if (lo <= 0) lo = hi > 0 ? hi / 1000 : 1;
        const int k0 = static_cast<int>(std::floor(std::log10(lo) + 1e-9)), k1 = static_cast<int>(std::ceil(std::log10(hi) - 1e-9));
        t.lo = std::pow(10.0, k0);
        t.hi = std::pow(10.0, std::max(k1, k0 + 1));
        for (int k = k0; k <= std::max(k1, k0 + 1); ++k) {
            t.at.push_back(std::pow(10.0, k));
            t.label.push_back(k >= -3 && k <= 5 ? tick_label(std::pow(10.0, k), k < 0 ? -k : 0) : "1e" + std::to_string(k));
        }
        return t;
    }
    double amin, amax, step;
    if (a.fixed) {
        amin = a.min;
        amax = a.max;
        step = a.step > 0 ? a.step : 0;
        if (step <= 0) { double x0, x1; nice_axis(amin, amax, 6, x0, x1, step); }
    } else {
        nice_axis(dlo, dhi, 6, amin, amax, step);
    }
    t.lo = amin;
    t.hi = amax;
    const int dec = decimals_for(step);
    const int n = static_cast<int>(std::round((amax - amin) / step));
    for (int i = 0; i <= n; ++i) {
        const double v = amin + i * step;
        t.at.push_back(v);
        t.label.push_back(tick_label(v, dec));
    }
    return t;
}

LineStyle parse_line(const std::string& s) {
    if (s == "solid") return LineStyle::Solid;
    if (s == "dashed") return LineStyle::Dashed;
    if (s == "dotted") return LineStyle::Dotted;
    if (s == "dashdot") return LineStyle::DashDot;
    if (s == "none") return LineStyle::None;
    throw std::runtime_error("unknown line style '" + s + "'");
}

Marker parse_marker(const std::string& s) {
    static const std::map<std::string, Marker> m{{"none", Marker::None},         {"circle", Marker::Circle},
                                                 {"circle-open", Marker::CircleOpen}, {"square", Marker::Square},
                                                 {"square-open", Marker::SquareOpen}, {"triangle", Marker::Triangle},
                                                 {"triangle-open", Marker::TriangleOpen}, {"diamond", Marker::Diamond},
                                                 {"diamond-open", Marker::DiamondOpen}, {"cross", Marker::Cross}};
    auto it = m.find(s);
    if (it == m.end()) throw std::runtime_error("unknown marker '" + s + "'");
    return it->second;
}

std::vector<std::string> tokens(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool q = false, any = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"') { q = !q; any = true; continue; }
        // '#' starts a comment only as a separate word, so colours such as color=#0072B2 survive
        if (!q && c == '#' && !any && (i + 1 == line.size() || line[i + 1] == ' ' || line[i + 1] == '\t')) break;
        if (!q && (c == ' ' || c == '\t' || c == '\r')) {
            if (any) out.push_back(cur);
            cur.clear();
            any = false;
            continue;
        }
        cur += c;
        any = true;
    }
    if (any) out.push_back(cur);
    return out;
}

struct Csv {
    std::vector<std::string> head;
    std::vector<std::vector<double>> rows;
    int col(const std::string& n) const {
        for (std::size_t i = 0; i < head.size(); ++i)
            if (head[i] == n) return static_cast<int>(i);
        throw std::runtime_error("column '" + n + "' not found");
    }
};

Csv read_numeric_csv(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    Csv c;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> v;
        std::string cur;
        for (char ch : line) {
            if (ch == ',') { v.push_back(cur); cur.clear(); }
            else if (ch != '\r' && ch != '"') cur += ch;
        }
        v.push_back(cur);
        if (c.head.empty()) { c.head = v; continue; }
        std::vector<double> r;
        for (const auto& x : v) {
            try { r.push_back(std::stod(x)); } catch (...) { r.push_back(NAN); }
        }
        c.rows.push_back(r);
    }
    return c;
}

}  // namespace

std::vector<std::string> okabe_ito() { return {"#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#000000"}; }

void nice_axis(double lo, double hi, int target, double& amin, double& amax, double& step) {
    if (!(hi > lo)) { const double d = std::fabs(lo) > 0 ? std::fabs(lo) * 0.1 : 1.0; lo -= d; hi += d; }
    const double raw = (hi - lo) / std::max(1, target);
    const double mag = std::pow(10.0, std::floor(std::log10(raw)));
    double best = mag;
    for (double m : {1.0, 2.0, 2.5, 5.0, 10.0})
        if (m * mag >= raw * 0.999) { best = m * mag; break; }
    step = best;
    amin = std::floor(lo / step + 1e-9) * step;
    amax = std::ceil(hi / step - 1e-9) * step;
}

std::string render_svg(const Figure& f) {
    const double fs = f.font_pt * kPt;
    const int cols = std::max(1, f.columns);
    const int rows = static_cast<int>((f.panels.size() + static_cast<std::size_t>(cols) - 1) / static_cast<std::size_t>(cols));
    const double cw = f.width_mm / cols, ch = f.height_mm / std::max(1, rows);
    std::ostringstream o;
    o << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << fmt(f.width_mm) << "mm\" height=\"" << fmt(f.height_mm)
      << "mm\" viewBox=\"0 0 " << fmt(f.width_mm) << " " << fmt(f.height_mm) << "\" font-family=\"" << esc(f.font) << "\" font-size=\""
      << fmt(fs) << "\">\n<rect width=\"100%\" height=\"100%\" fill=\"#fff\"/>\n";
    for (std::size_t pi = 0; pi < f.panels.size(); ++pi) {
        const Panel& p = f.panels[pi];
        std::map<const Series*, std::string> scol;
        std::size_t pal = 0;
        for (const auto& se : p.series)
            scol[&se] = !se.color.empty() ? se.color : (f.palette.empty() ? std::string("#000") : f.palette[pal++ % f.palette.size()]);
        auto colour_of = [&](const Series& se) { return scol[&se]; };
        const double ox = cw * static_cast<double>(pi % static_cast<std::size_t>(cols));
        const double oy = ch * static_cast<double>(pi / static_cast<std::size_t>(cols));
        // data range
        double xlo = 1e300, xhi = -1e300, ylo = 1e300, yhi = -1e300;
        for (const auto& s : p.series)
            for (std::size_t i = 0; i < s.x.size(); ++i) {
                if (!std::isfinite(s.x[i]) || !std::isfinite(s.y[i])) continue;
                if (p.x.log && s.x[i] <= 0) continue;
                if (p.y.log && s.y[i] <= 0) continue;
                xlo = std::min(xlo, s.x[i]); xhi = std::max(xhi, s.x[i]);
                ylo = std::min(ylo, s.y[i]); yhi = std::max(yhi, s.y[i]);
            }
        for (const auto& b : p.bands)
            for (std::size_t i = 0; i < b.x.size(); ++i) {
                if (!std::isfinite(b.lo[i]) || !std::isfinite(b.hi[i])) continue;
                xlo = std::min(xlo, b.x[i]); xhi = std::max(xhi, b.x[i]);
                ylo = std::min(ylo, b.lo[i]); yhi = std::max(yhi, b.hi[i]);
            }
        for (const auto& r : p.refs) {
            if (r.horizontal) { ylo = std::min(ylo, r.value); yhi = std::max(yhi, r.value); }
            else { xlo = std::min(xlo, r.value); xhi = std::max(xhi, r.value); }
        }
        if (xlo > xhi) { xlo = 0; xhi = 1; }
        if (ylo > yhi) { ylo = 0; yhi = 1; }
        const Ticks tx = make_ticks(p.x, xlo, xhi), ty = make_ticks(p.y, ylo, yhi);
        // legend layout (below)
        std::vector<const Series*> leg;
        for (const auto& s : p.series)
            if (!s.label.empty()) leg.push_back(&s);
        std::vector<const Band*> legb;
        for (const auto& b : p.bands)
            if (!b.label.empty()) legb.push_back(&b);
        std::vector<const RefLine*> legr;
        for (const auto& r : p.refs)
            if (!r.label.empty() && p.legend != "none") legr.push_back(&r);
        const double entry_w = 8.0;
        auto text_w = [&](const std::string& s) { return 0.52 * fs * static_cast<double>(s.size()); };
        int leg_lines = 0;
        if (p.legend == "below" && (!leg.empty() || !legr.empty() || !legb.empty())) {
            double x = 0;
            leg_lines = 1;
            const double avail = cw - 8;
            int count = 0;
            auto place = [&](const std::string& lab) {
                const double w = entry_w + text_w(lab) + 3;
                const bool full = p.legend_columns > 0 && count > 0 && count % p.legend_columns == 0;
                if (x > 0 && (full || x + w > avail)) { ++leg_lines; x = 0; }
                x += w;
                ++count;
            };
            for (auto* b : legb) place(b->label);
            for (auto* s : leg) place(s->label);
            for (auto* r : legr) place(r->label);
        }
        if (p.legend == "below") leg_lines = std::max(leg_lines, p.legend_rows);  // same plot height as the panel beside it
        // margins
        double ytw = 0;
        for (const auto& l : ty.label) ytw = std::max(ytw, text_w(l));
        const double ml = ytw + 2.0 + (p.y.label.empty() ? 0 : fs * 1.6), mr = 3.0;
        const double mt = 2.0 + (p.title.empty() ? 0 : fs * 1.5);
        const double mb = fs * 1.6 + (p.x.label.empty() ? 0 : fs * 1.4) + (p.letter.empty() ? 0 : fs * 1.6) + leg_lines * fs * 1.5 + 1.0;
        const double X0 = ox + ml, X1 = ox + cw - mr, Y0 = oy + mt, Y1 = oy + ch - mb;
        auto mapx = [&](double v) {
            return tx.log ? X0 + (std::log10(v) - std::log10(tx.lo)) / (std::log10(tx.hi) - std::log10(tx.lo)) * (X1 - X0)
                          : X0 + (v - tx.lo) / (tx.hi - tx.lo) * (X1 - X0);
        };
        auto mapy = [&](double v) {
            return ty.log ? Y1 - (std::log10(v) - std::log10(ty.lo)) / (std::log10(ty.hi) - std::log10(ty.lo)) * (Y1 - Y0)
                          : Y1 - (v - ty.lo) / (ty.hi - ty.lo) * (Y1 - Y0);
        };
        const std::string clip = "clip" + std::to_string(pi);
        o << "<defs><clipPath id=\"" << clip << "\"><rect x=\"" << fmt(X0) << "\" y=\"" << fmt(Y0) << "\" width=\"" << fmt(X1 - X0)
          << "\" height=\"" << fmt(Y1 - Y0) << "\"/></clipPath></defs>\n";
        if (!p.title.empty())
            o << "<text x=\"" << fmt((X0 + X1) / 2) << "\" y=\"" << fmt(oy + mt - 1.2) << "\" text-anchor=\"middle\">" << esc(p.title) << "</text>\n";
        // frame and ticks
        o << "<rect x=\"" << fmt(X0) << "\" y=\"" << fmt(Y0) << "\" width=\"" << fmt(X1 - X0) << "\" height=\"" << fmt(Y1 - Y0)
          << "\" fill=\"none\" stroke=\"#000\" stroke-width=\"0.2\"/>\n";
        for (std::size_t i = 0; i < tx.at.size(); ++i) {
            const double x = mapx(tx.at[i]);
            o << "<line x1=\"" << fmt(x) << "\" y1=\"" << fmt(Y1) << "\" x2=\"" << fmt(x) << "\" y2=\"" << fmt(Y1 - 1.0)
              << "\" stroke=\"#000\" stroke-width=\"0.2\"/>";
            o << "<text x=\"" << fmt(x) << "\" y=\"" << fmt(Y1 + fs * 1.15) << "\" text-anchor=\"middle\">" << esc(tx.label[i]) << "</text>\n";
        }
        for (std::size_t i = 0; i < ty.at.size(); ++i) {
            const double y = mapy(ty.at[i]);
            o << "<line x1=\"" << fmt(X0) << "\" y1=\"" << fmt(y) << "\" x2=\"" << fmt(X0 + 1.0) << "\" y2=\"" << fmt(y)
              << "\" stroke=\"#000\" stroke-width=\"0.2\"/>";
            o << "<text x=\"" << fmt(X0 - 1.0) << "\" y=\"" << fmt(y + fs * 0.35) << "\" text-anchor=\"end\">" << esc(ty.label[i]) << "</text>\n";
        }
        if (tx.log || ty.log) {  // minor ticks on log axes
            auto minor = [&](const Ticks& t, bool xaxis) {
                for (std::size_t i = 0; i + 1 < t.at.size(); ++i)
                    for (int k = 2; k <= 9; ++k) {
                        const double v = t.at[i] * k;
                        if (xaxis) { const double x = mapx(v); o << "<line x1=\"" << fmt(x) << "\" y1=\"" << fmt(Y1) << "\" x2=\"" << fmt(x) << "\" y2=\"" << fmt(Y1 - 0.5) << "\" stroke=\"#000\" stroke-width=\"0.15\"/>"; }
                        else { const double y = mapy(v); o << "<line x1=\"" << fmt(X0) << "\" y1=\"" << fmt(y) << "\" x2=\"" << fmt(X0 + 0.5) << "\" y2=\"" << fmt(y) << "\" stroke=\"#000\" stroke-width=\"0.15\"/>"; }
                    }
            };
            if (tx.log) minor(tx, true);
            if (ty.log) minor(ty, false);
            o << "\n";
        }
        if (!p.x.label.empty())
            o << "<text x=\"" << fmt((X0 + X1) / 2) << "\" y=\"" << fmt(Y1 + fs * 2.6) << "\" text-anchor=\"middle\">" << esc(p.x.label) << "</text>\n";
        if (!p.y.label.empty()) {
            const double xl = ox + fs * 0.9, yl = (Y0 + Y1) / 2;
            o << "<text transform=\"translate(" << fmt(xl) << "," << fmt(yl) << ") rotate(-90)\" text-anchor=\"middle\">" << esc(p.y.label) << "</text>\n";
        }
        if (f.grid) {
            for (double v : tx.at) {
                const double x = mapx(v);
                if (x > X0 + 0.05 && x < X1 - 0.05)
                    o << "<line x1=\"" << fmt(x) << "\" y1=\"" << fmt(Y0) << "\" x2=\"" << fmt(x) << "\" y2=\"" << fmt(Y1) << "\" stroke=\"#d9d9d9\" stroke-width=\"0.12\"/>";
            }
            for (double v : ty.at) {
                const double y = mapy(v);
                if (y > Y0 + 0.05 && y < Y1 - 0.05)
                    o << "<line x1=\"" << fmt(X0) << "\" y1=\"" << fmt(y) << "\" x2=\"" << fmt(X1) << "\" y2=\"" << fmt(y) << "\" stroke=\"#d9d9d9\" stroke-width=\"0.12\"/>";
            }
            o << "\n";
        }
        if (!p.bands.empty()) {
            o << "<g clip-path=\"url(#" << clip << ")\">\n";
            for (const auto& b : p.bands) {
                std::string d;
                for (std::size_t i = 0; i < b.x.size(); ++i) d += (i ? "L" : "M") + fmt(mapx(b.x[i])) + " " + fmt(mapy(b.hi[i]));
                for (std::size_t i = b.x.size(); i-- > 0;) d += "L" + fmt(mapx(b.x[i])) + " " + fmt(mapy(b.lo[i]));
                o << "<path d=\"" << d << "Z\" fill=\"" << b.color << "\" fill-opacity=\"" << fmt(b.opacity) << "\" stroke=\"none\"/>\n";
            }
            o << "</g>\n";
        }
        // reference lines
        for (const auto& r : p.refs) {
            if (r.horizontal) {
                const double y = mapy(r.value);
                o << "<line x1=\"" << fmt(X0) << "\" y1=\"" << fmt(y) << "\" x2=\"" << fmt(X1) << "\" y2=\"" << fmt(y)
                  << "\" stroke=\"" << r.color << "\" stroke-width=\"0.16\"" << dash(r.line, 0.16) << "/>\n";
            } else {
                const double x = mapx(r.value);
                o << "<line x1=\"" << fmt(x) << "\" y1=\"" << fmt(Y0) << "\" x2=\"" << fmt(x) << "\" y2=\"" << fmt(Y1)
                  << "\" stroke=\"" << r.color << "\" stroke-width=\"0.16\"" << dash(r.line, 0.16) << "/>\n";
            }
        }
        // series
        o << "<g clip-path=\"url(#" << clip << ")\">\n";
        for (const auto& s : p.series) {
            if (s.line != LineStyle::None) {
                std::string d;
                bool pen = false;
                for (std::size_t i = 0; i < s.x.size(); ++i) {
                    const bool ok = std::isfinite(s.x[i]) && std::isfinite(s.y[i]) && (!tx.log || s.x[i] > 0) && (!ty.log || s.y[i] > 0);
                    if (!ok) { pen = false; continue; }
                    d += (pen ? "L" : "M") + fmt(mapx(s.x[i])) + " " + fmt(mapy(s.y[i]));
                    pen = true;
                }
                o << "<path d=\"" << d << "\" fill=\"none\" stroke=\"" << colour_of(s) << "\" stroke-width=\"" << fmt(s.width) << "\" stroke-linejoin=\"round\""
                  << dash(s.line, s.width) << "/>\n";
            }
        }
        o << "</g>\n";
        for (const auto& s : p.series) {
            if (s.marker == Marker::None) continue;
            for (std::size_t i = 0; i < s.x.size(); i += static_cast<std::size_t>(std::max(1, s.marker_every))) {
                if (!std::isfinite(s.x[i]) || !std::isfinite(s.y[i])) continue;
                const double x = mapx(s.x[i]), y = mapy(s.y[i]);
                if (x < X0 - 0.01 || x > X1 + 0.01 || y < Y0 - 0.01 || y > Y1 + 0.01) continue;
                o << marker_svg(s.marker, x, y, 0.7, colour_of(s));
            }
            o << "\n";
        }
        // legend
        auto legend_entry = [&](double x, double y, LineStyle ls, Marker mk, double w, const std::string& lab, const std::string& col) {
            if (ls != LineStyle::None)
                o << "<line x1=\"" << fmt(x) << "\" y1=\"" << fmt(y) << "\" x2=\"" << fmt(x + 6) << "\" y2=\"" << fmt(y)
                  << "\" stroke=\"" << col << "\" stroke-width=\"" << fmt(w) << "\"" << dash(ls, w) << "/>";
            if (mk != Marker::None) o << marker_svg(mk, x + 3, y, 0.7, col);
            o << "<text x=\"" << fmt(x + entry_w - 1) << "\" y=\"" << fmt(y + fs * 0.35) << "\">" << esc(lab) << "</text>\n";
        };
        if (p.legend == "below" && leg_lines > 0) {
            double x = 0;
            int line = 0;
            const double ybase = Y1 + fs * 1.6 + (p.x.label.empty() ? 0 : fs * 1.4) + fs * 0.9;
            const double avail = cw - 8;
            int count = 0;
            auto next_entry = [&](double ew) {  // new line when the row is full or has legend_columns entries
                const bool full = p.legend_columns > 0 && count > 0 && count % p.legend_columns == 0;
                if (x > 0 && (full || x + ew > avail)) { ++line; x = 0; }
                ++count;
            };
            auto place = [&](LineStyle ls, Marker mk, double w, const std::string& lab, const std::string& col) {
                const double ew = entry_w + text_w(lab) + 3;
                next_entry(ew);
                legend_entry(X0 + x, ybase + line * fs * 1.5, ls, mk, w, lab, col);
                x += ew;
            };
            for (auto* b : legb) {  // band: a filled swatch
                const double ew = entry_w + text_w(b->label) + 3;
                next_entry(ew);
                const double yy = ybase + line * fs * 1.5;
                o << "<rect x=\"" << fmt(X0 + x) << "\" y=\"" << fmt(yy - 1.0) << "\" width=\"6\" height=\"2\" fill=\"" << b->color
                  << "\" fill-opacity=\"" << fmt(b->opacity) << "\"/><text x=\"" << fmt(X0 + x + entry_w - 1) << "\" y=\"" << fmt(yy + fs * 0.35)
                  << "\">" << esc(b->label) << "</text>\n";
                x += ew;
            }
            for (auto* s : leg) place(s->line, s->marker, s->width, s->label, colour_of(*s));
            for (auto* r : legr) place(r->line, Marker::None, 0.2, r->label, r->color);
        } else if (p.legend == "inside" && (!leg.empty() || !legr.empty())) {
            double wmax = 0;
            for (auto* s : leg) wmax = std::max(wmax, text_w(s->label));
            for (auto* r : legr) wmax = std::max(wmax, text_w(r->label));
            const double bw = entry_w + wmax + 1.5, bh = (static_cast<double>(leg.size() + legr.size())) * fs * 1.4 + 1.0;
            const double bx = X1 - bw - 1.0, by = Y0 + 1.0;
            o << "<rect x=\"" << fmt(bx) << "\" y=\"" << fmt(by) << "\" width=\"" << fmt(bw) << "\" height=\"" << fmt(bh) << "\" fill=\"#fff\"/>\n";
            double y = by + fs * 0.9;
            for (auto* s : leg) { legend_entry(bx + 1, y, s->line, s->marker, s->width, s->label, colour_of(*s)); y += fs * 1.4; }
            for (auto* r : legr) { legend_entry(bx + 1, y, r->line, Marker::None, 0.2, r->label, r->color); y += fs * 1.4; }
        }
        if (!p.letter.empty())
            o << "<text x=\"" << fmt((X0 + X1) / 2) << "\" y=\"" << fmt(oy + ch - 1.0) << "\" text-anchor=\"middle\">" << esc(p.letter) << "</text>\n";
    }
    o << "</svg>\n";
    return o.str();
}

Figure figure_from_spec(const std::string& text, const std::string& base_dir) {
    Figure f;
    std::istringstream is(text);
    std::string line;
    int ln = 0;
    std::map<std::string, Csv> cache;
    auto path_of = [&](const std::string& p) {
        if (!p.empty() && (p[0] == '/' || (p.size() > 1 && p[1] == ':'))) return p;
        return base_dir + "/" + p;
    };
    auto panel = [&]() -> Panel& {
        if (f.panels.empty()) f.panels.emplace_back();
        return f.panels.back();
    };
    auto axis = [&](Axis& a, const std::vector<std::string>& t) {
        a.fixed = true;
        a.min = std::stod(t.at(1));
        a.max = std::stod(t.at(2));
        a.step = t.size() > 3 ? std::stod(t[3]) : 0;
    };
    while (std::getline(is, line)) {
        ++ln;
        auto t = tokens(line);
        if (t.empty()) continue;
        try {
            const auto& k = t[0];
            std::map<std::string, std::string> kv;
            std::vector<std::string> pos;
            for (std::size_t i = 1; i < t.size(); ++i) {
                auto e = t[i].find('=');
                if (e != std::string::npos && e > 0 && k != "title" && k != "xlabel" && k != "ylabel") kv[t[i].substr(0, e)] = t[i].substr(e + 1);
                else pos.push_back(t[i]);
            }
            if (k == "size") { f.width_mm = std::stod(t.at(1)); f.height_mm = std::stod(t.at(2)); }
            else if (k == "columns") f.columns = std::stoi(t.at(1));
            else if (k == "font_pt") f.font_pt = std::stod(t.at(1));
            else if (k == "font") f.font = t.at(1);
            else if (k == "panel") { f.panels.emplace_back(); if (!pos.empty()) f.panels.back().letter = pos[0]; }
            else if (k == "title") panel().title = t.at(1);
            else if (k == "xlabel") panel().x.label = t.at(1);
            else if (k == "ylabel") panel().y.label = t.at(1);
            else if (k == "xlog") panel().x.log = true;
            else if (k == "ylog") panel().y.log = true;
            else if (k == "xrange") axis(panel().x, t);
            else if (k == "yrange") axis(panel().y, t);
            else if (k == "legend") panel().legend = t.at(1);
            else if (k == "legend_columns") panel().legend_columns = std::stoi(t.at(1));
            else if (k == "legend_rows") panel().legend_rows = std::stoi(t.at(1));
            else if (k == "palette") {
                if (t.at(1) == "okabe-ito") f.palette = okabe_ito();
                else if (t.at(1) == "black") f.palette.clear();
                else f.palette.assign(t.begin() + 1, t.end());
            } else if (k == "grid") f.grid = t.size() < 2 || t[1] != "off";
            else if (k == "band") {
                // band FILE XCOL LOCOL HICOL [label=] [color=] [opacity=] [xscale=] [xoffset=] [yscale=]
                const std::string file = path_of(pos.at(0));
                if (!cache.count(file)) cache[file] = read_numeric_csv(file);
                const Csv& c = cache[file];
                const int cx = c.col(pos.at(1)), cl = c.col(pos.at(2)), chh = c.col(pos.at(3));
                Band b;
                b.label = kv.count("label") ? kv["label"] : "";
                if (kv.count("color")) b.color = kv["color"];
                if (kv.count("opacity")) b.opacity = std::stod(kv["opacity"]);
                const double xs = kv.count("xscale") ? std::stod(kv["xscale"]) : 1.0, xo = kv.count("xoffset") ? std::stod(kv["xoffset"]) : 0.0;
                const double ys = kv.count("yscale") ? std::stod(kv["yscale"]) : 1.0;
                for (const auto& r : c.rows) {
                    b.x.push_back(r[static_cast<std::size_t>(cx)] * xs + xo);
                    b.lo.push_back(r[static_cast<std::size_t>(cl)] * ys);
                    b.hi.push_back(r[static_cast<std::size_t>(chh)] * ys);
                }
                panel().bands.push_back(b);
            }
            else if (k == "series") {
                // series FILE XCOL YCOL [label= line= marker= every= width= xscale= yscale= xoffset= where=COL:VALUE]
                const std::string file = path_of(pos.at(0));
                if (!cache.count(file)) cache[file] = read_numeric_csv(file);
                const Csv& c = cache[file];
                const int cx = c.col(pos.at(1)), cy = c.col(pos.at(2));
                Series s;
                s.label = kv.count("label") ? kv["label"] : "";
                s.line = kv.count("line") ? parse_line(kv["line"]) : LineStyle::Solid;
                s.marker = kv.count("marker") ? parse_marker(kv["marker"]) : Marker::None;
                s.marker_every = kv.count("every") ? std::stoi(kv["every"]) : 1;
                s.width = kv.count("width") ? std::stod(kv["width"]) : 0.2;
                if (kv.count("color")) s.color = kv["color"];
                const double xs = kv.count("xscale") ? std::stod(kv["xscale"]) : 1.0, ys = kv.count("yscale") ? std::stod(kv["yscale"]) : 1.0;
                const double xo = kv.count("xoffset") ? std::stod(kv["xoffset"]) : 0.0, yo = kv.count("yoffset") ? std::stod(kv["yoffset"]) : 0.0;
                int wcol = -1;
                double wval = 0;
                if (kv.count("where")) {
                    auto c2 = kv["where"].find(':');
                    if (c2 == std::string::npos) throw std::runtime_error("where=COLUMN:VALUE");
                    wcol = c.col(kv["where"].substr(0, c2));
                    wval = std::stod(kv["where"].substr(c2 + 1));
                }
                for (const auto& r : c.rows) {
                    if (wcol >= 0 && std::fabs(r[static_cast<std::size_t>(wcol)] - wval) > 1e-6 * std::max(1.0, std::fabs(wval))) continue;
                    s.x.push_back(r[static_cast<std::size_t>(cx)] * xs + xo);
                    s.y.push_back(r[static_cast<std::size_t>(cy)] * ys + yo);
                }
                panel().series.push_back(s);
            } else if (k == "point") {  // point X Y [label=] [marker=]  (single published value etc.)
                Series s;
                s.x = {std::stod(pos.at(0))};
                s.y = {std::stod(pos.at(1))};
                s.label = kv.count("label") ? kv["label"] : "";
                s.line = LineStyle::None;
                s.marker = kv.count("marker") ? parse_marker(kv["marker"]) : Marker::CircleOpen;
                s.color = kv.count("color") ? kv["color"] : "#000";
                panel().series.push_back(s);
            } else if (k == "hline" || k == "vline") {
                RefLine r;
                r.horizontal = k == "hline";
                r.value = std::stod(pos.at(0));
                r.label = kv.count("label") ? kv["label"] : "";
                r.line = kv.count("line") ? parse_line(kv["line"]) : LineStyle::Dashed;
                if (kv.count("color")) r.color = kv["color"];
                panel().refs.push_back(r);
            } else {
                throw std::runtime_error("unknown keyword '" + k + "'");
            }
        } catch (const std::exception& e) {
            throw std::runtime_error("plot spec line " + std::to_string(ln) + ": " + e.what());
        }
    }
    return f;
}

// ------------------------------------------------------------------ mechanism drawing
namespace {
std::vector<Vec2> hull(std::vector<Vec2> p) {
    std::sort(p.begin(), p.end());
    p.erase(std::unique(p.begin(), p.end(), [](const Vec2& a, const Vec2& b) { return std::hypot(a[0] - b[0], a[1] - b[1]) < 1e-6; }), p.end());
    if (p.size() < 3) return p;
    auto cross = [](const Vec2& o, const Vec2& a, const Vec2& b) { return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]); };
    std::vector<Vec2> h(2 * p.size());
    std::size_t k = 0;
    for (std::size_t i = 0; i < p.size(); ++i) {
        while (k >= 2 && cross(h[k - 2], h[k - 1], p[i]) <= 0) --k;
        h[k++] = p[i];
    }
    for (std::size_t i = p.size() - 1, t = k + 1; i > 0; --i) {
        while (k >= t && cross(h[k - 2], h[k - 1], p[i - 1]) <= 0) --k;
        h[k++] = p[i - 1];
    }
    h.resize(k - 1);
    return h;
}
}  // namespace

std::string draw_mechanism(const Kinematics& kin, const std::vector<Pose>& poses, double W, double H, bool labels, double font_pt, bool colour) {
    DrawOptions o;
    o.width_mm = W; o.height_mm = H; o.labels = labels; o.font_pt = font_pt; o.colour = colour;
    return draw_mechanism(kin, poses, o);
}

std::string draw_mechanism(const Kinematics& kin, const std::vector<Pose>& poses, const DrawOptions& opt) {
    double W = opt.width_mm, H = opt.height_mm;
    const bool labels = opt.labels, colour = opt.colour;
    const double font_pt = opt.font_pt, scale_mm_per_m = opt.scale_mm_per_m, y_bottom = opt.y_bottom;
    // colour scheme: links in light blue-grey with a dark outline, cylinders in vermilion (they are the subject)
    const Model& m = kin.model();
    const double fs = font_pt * kPt;
    const double rr = 1.0;  // rounding radius of the links on paper [mm]
    std::vector<std::vector<Vec2>> bp(m.bodies.size());
    for (const auto& j : m.revolutes) { bp[static_cast<std::size_t>(j.body_a)].push_back(j.p); bp[static_cast<std::size_t>(j.body_b)].push_back(j.p); }
    for (const auto& j : m.prismatics) { bp[static_cast<std::size_t>(j.body_a)].push_back(j.p); bp[static_cast<std::size_t>(j.body_b)].push_back(j.p); }
    for (const auto& a : m.actuators) {
        if (!a.has_cylinder) continue;
        bp[static_cast<std::size_t>(a.body_a)].push_back(a.pa);
        bp[static_cast<std::size_t>(a.body_b)].push_back(a.pb);
    }
    for (const auto& [b, pts] : m.outlines)
        for (const auto& q : pts) bp[static_cast<std::size_t>(b)].push_back(q);
    for (std::size_t b = 1; b < bp.size(); ++b)
        if (bp[b].size() < 2)
            for (const auto& mi : m.masses)
                if (mi.body == static_cast<int>(b)) bp[b].push_back(mi.c);
    double xlo = 1e300, xhi = -1e300, ylo = 1e300, yhi = -1e300;
    for (const auto& s : poses)
        for (std::size_t b = 0; b < m.bodies.size(); ++b)
            for (const auto& p : bp[b]) {
                const Vec2 q = kin.point(s, static_cast<int>(b), p);
                xlo = std::min(xlo, q[0]); xhi = std::max(xhi, q[0]);
                ylo = std::min(ylo, q[1]); yhi = std::max(yhi, q[1]);
            }
    if (opt.view.size() == 4) { xlo = opt.view[0]; ylo = opt.view[1]; xhi = opt.view[2]; yhi = opt.view[3]; }
    const double ml = 5.0, mb = 9.0;  // margins: sides and bottom (ground and scale bar)
    const double mt = opt.title.empty() ? ml : ml + 1.6 * font_pt * kPt;  // top margin
    double sc;
    if (scale_mm_per_m > 0) {
        sc = scale_mm_per_m / 1000.0;
        if (std::isfinite(y_bottom)) ylo = opt.view.size() == 4 ? y_bottom : std::min(ylo, y_bottom);
        W = sc * (xhi - xlo) + 2 * ml;
        H = std::max(H, sc * (yhi - ylo) + mt + mb);  // extra height goes above the mechanism
    } else {
        sc = std::min((W - 2 * ml) / std::max(1.0, xhi - xlo), (H - mt - mb) / std::max(1.0, yhi - ylo));
    }
    const double cx = W / 2 - sc * (xlo + xhi) / 2;
    const double cy = scale_mm_per_m > 0 ? H - mb + sc * ylo : mt + (H - mt - mb) / 2 + sc * (ylo + yhi) / 2;
    auto X = [&](const Vec2& q) { return cx + sc * q[0]; };
    auto Y = [&](const Vec2& q) { return cy - sc * q[1]; };
    struct P2 { double x, y; double w = 3.0; };
    std::vector<P2> obst;  // paper points that labels should avoid
    auto seg_obst = [&](P2 a, P2 b, double w = 3.0) {  // weight: 3 for marks that must stay visible, less for outlines
        const double L = std::hypot(b.x - a.x, b.y - a.y);
        const int n = std::max(1, static_cast<int>(L / 0.6));
        for (int i = 0; i <= n; ++i) obst.push_back({a.x + (b.x - a.x) * i / n, a.y + (b.y - a.y) * i / n, w});
    };
    std::ostringstream o, lab;
    o << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << fmt(W) << "mm\" height=\"" << fmt(H) << "mm\" viewBox=\"0 0 " << fmt(W) << " "
      << fmt(H) << "\" font-family=\"Arial, Helvetica, sans-serif\" font-size=\"" << fmt(fs) << "\">\n<rect width=\"100%\" height=\"100%\" fill=\"#fff\"/>\n";
    if (opt.view.size() == 4)
        o << "<clipPath id=\"view\"><rect x=\"" << fmt(X({xlo, 0})) << "\" y=\"" << fmt(Y({0, yhi})) << "\" width=\"" << fmt(sc * (xhi - xlo))
          << "\" height=\"" << fmt(sc * (yhi - ylo)) << "\"/></clipPath>\n<g clip-path=\"url(#view)\">\n";
    std::vector<std::vector<P2>> barrel_pts(m.actuators.size());
    std::map<int, std::vector<std::vector<P2>>> body_polys;  // drawn faces of each link (first pose), for item anchors
    struct Lbl { std::string text; P2 at; std::vector<P2> cand; std::vector<double> pen; int body = -1; };
    std::vector<Lbl> lbls;
    std::vector<std::pair<int, std::vector<P2>>> hulls;  // links of the first pose, paper coordinates
    auto inside = [](const std::vector<P2>& h, P2 q) {  // inside a convex polygon of either orientation
        if (h.size() < 3) return false;
        int sg = 0;
        for (std::size_t i = 0; i < h.size(); ++i) {
            const P2 a = h[i], c = h[(i + 1) % h.size()];
            const double cr = (c.x - a.x) * (q.y - a.y) - (c.y - a.y) * (q.x - a.x);
            const int si = cr > 0 ? 1 : (cr < 0 ? -1 : 0);
            if (si == 0) continue;
            if (sg == 0) sg = si; else if (si != sg) return false;
        }
        return true;
    };
    for (std::size_t pi = 0; pi < poses.size(); ++pi) {
        const Pose& s = poses[pi];
        const bool main_pose = pi == 0;
        const std::string ds = main_pose ? "" : " stroke-dasharray=\"1.2,0.8\"";
        const double lw = main_pose ? 0.35 : 0.18;
        o << "<g>\n";
        // ground: one base line under all ground joints (not when the solid model or a sketch draws the frame itself)
        bool frame_drawn = false;
        for (const auto& pt : opt.parts) frame_drawn = frame_drawn || (pt.body == 0 && pt.actuator < 0);
        if (main_pose && !frame_drawn) {
            std::vector<P2> g;
            for (const auto& j : m.revolutes)
                if (j.body_a == 0 || j.body_b == 0) { const Vec2 q = kin.point(s, j.body_b, j.p); g.push_back({X(q), Y(q)}); }
            if (!g.empty()) {
                double gx0 = 1e9, gx1 = -1e9, gy = -1e9;
                for (auto& q : g) { gx0 = std::min(gx0, q.x); gx1 = std::max(gx1, q.x); gy = std::max(gy, q.y); }
                gy += 2.2; gx0 -= 2.6; gx1 += 2.6;
                for (auto& q : g) {
                    o << "<polygon points=\"" << fmt(q.x) << "," << fmt(q.y) << " " << fmt(q.x - 1.4) << "," << fmt(gy) << " " << fmt(q.x + 1.4) << ","
                      << fmt(gy) << "\" fill=\"#fff\" stroke=\"#000\" stroke-width=\"0.2\" stroke-linejoin=\"round\"/>\n";
                    seg_obst({q.x - 1.4, gy}, {q.x + 1.4, gy});
                }
                o << "<line x1=\"" << fmt(gx0) << "\" y1=\"" << fmt(gy) << "\" x2=\"" << fmt(gx1) << "\" y2=\"" << fmt(gy) << "\" stroke=\"#000\" stroke-width=\"0.3\"/>\n";
                for (double hx = gx0 + 0.9; hx <= gx1 + 1e-9; hx += 0.9)
                    o << "<line x1=\"" << fmt(hx) << "\" y1=\"" << fmt(gy) << "\" x2=\"" << fmt(hx - 0.8) << "\" y2=\"" << fmt(gy + 0.9)
                      << "\" stroke=\"#000\" stroke-width=\"0.15\"/>";
                o << "\n";
                seg_obst({gx0, gy}, {gx1, gy});
            }
        }
        // solid-model parts (first pose): rigid motion of each part from the reference configuration
        const bool have_parts = !opt.parts.empty();
        std::vector<bool> act_drawn(m.actuators.size(), false);
        if (have_parts && main_pose) {
            const Pose ref = kin.reference();
            struct Tf { double c, s, tx, ty; };
            auto body_tf = [&](int b) {
                const Vec2 q0 = kin.point(s, b, {0, 0}), q1 = kin.point(s, b, {1000, 0});
                const double ang = std::atan2(q1[1] - q0[1], q1[0] - q0[0]);
                return Tf{std::cos(ang), std::sin(ang), q0[0], q0[1]};
            };
            auto act_tf = [&](int k, bool rod) {
                const auto& a = m.actuators[static_cast<std::size_t>(k)];
                const Vec2 u0 = kin.actuator_dir(ref, k), u1 = kin.actuator_dir(s, k);
                const double ang = std::atan2(u1[1], u1[0]) - std::atan2(u0[1], u0[0]);
                const Vec2 p0 = rod ? kin.point(ref, a.body_b, a.pb) : kin.point(ref, a.body_a, a.pa);
                const Vec2 p1 = rod ? kin.point(s, a.body_b, a.pb) : kin.point(s, a.body_a, a.pa);
                const double c = std::cos(ang), sn = std::sin(ang);
                return Tf{c, sn, p1[0] - (c * p0[0] - sn * p0[1]), p1[1] - (sn * p0[0] + c * p0[1])};
            };
            std::vector<std::size_t> idx(opt.parts.size());
            for (std::size_t i = 0; i < idx.size(); ++i) idx[i] = i;
            std::stable_sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return opt.parts[a].depth < opt.parts[b].depth; });
            static const char* tint[] = {"#CFE3F3", "#CDEBDD", "#F2D7E6", "#D6ECF8", "#F8E3BF", "#E6E6E6"};
            for (std::size_t ii : idx) {
                const DrawPart& pt = opt.parts[ii];
                const Tf t = pt.actuator >= 0 ? act_tf(pt.actuator, pt.rod) : body_tf(pt.body);
                if (pt.actuator >= 0) act_drawn[static_cast<std::size_t>(pt.actuator)] = true;
                auto tp = [&](const Vec2& q) { return P2{X({t.c * q[0] - t.s * q[1] + t.tx, 0}), Y({0, t.s * q[0] + t.c * q[1] + t.ty})}; };
                std::string fill;
                if (pt.actuator >= 0) fill = pt.rod ? "#F4F4F4" : (colour ? "#E8875A" : "#fff");
                else if (pt.body == 0) fill = "#DADADA";
                else if (m.outlines.count(pt.body)) fill = "#E9E9E9";
                else fill = colour ? tint[(pt.body - 1) % 5] : "#fff";
                o << "<g fill=\"" << fill << "\" stroke=\"none\">";
                for (std::size_t fi = 0; fi < pt.faces.size(); ++fi) {
                    const auto& f = pt.faces[fi];
                    o << "<path fill-rule=\"evenodd\" d=\"";
                    if (pt.actuator < 0) {  // outer loop (or hull) of the face, for item anchors
                        std::vector<Vec2> src = f.front();
                        if (fi < pt.curved.size() && pt.curved[fi]) {
                            std::vector<Vec2> all;
                            for (const auto& lp : f) all.insert(all.end(), lp.begin(), lp.end());
                            src = hull(all);
                        }
                        std::vector<P2> poly;
                        for (const auto& q : src) poly.push_back(tp(q));
                        if (poly.size() >= 3) body_polys[pt.body].push_back(std::move(poly));
                    }
                    if (fi < pt.curved.size() && pt.curved[fi]) {
                        std::vector<Vec2> all;
                        for (const auto& lp : f) all.insert(all.end(), lp.begin(), lp.end());
                        const auto hh = hull(all);
                        for (std::size_t i = 0; i < hh.size(); ++i) {
                            const P2 q = tp(hh[i]);
                            o << (i ? "L" : "M") << fmt(q.x) << " " << fmt(q.y);
                        }
                        o << "Z";
                    } else {
                        for (const auto& lp : f) {
                            for (std::size_t i = 0; i < lp.size(); ++i) {
                                const P2 q = tp(lp[i]);
                                o << (i ? "L" : "M") << fmt(q.x) << " " << fmt(q.y);
                            }
                            o << "Z";
                        }
                    }
                    o << "\"/>";
                }
                // line work without faces (a 2D sketch): coloured by link, the frame in grey
                static const char* sk_col[] = {"#0072B2", "#009E73", "#CC79A7", "#56B4E9", "#E69F00", "#000000"};
                const bool sketch = pt.faces.empty() && pt.actuator < 0;
                const std::string edge_col = !sketch ? "#2b2b2b" : !colour ? "#000" : pt.body == 0 ? "#8C8C8C" : sk_col[(pt.body - 1) % 6];
                o << "</g>\n<path fill=\"none\" stroke=\"" << edge_col << "\" stroke-width=\"" << (sketch ? (pt.body == 0 ? "0.08" : "0.13") : "0.05")
                  << "\" stroke-linejoin=\"round\" stroke-linecap=\"round\" d=\"";
                for (const auto& e : pt.edges) {
                    for (std::size_t i = 0; i < e.size(); ++i) {
                        const P2 q = tp(e[i]);
                        o << (i ? "L" : "M") << fmt(q.x) << " " << fmt(q.y);
                        if (i % 4 == 0) obst.push_back({q.x, q.y, 0.02});
                    }
                }
                o << "\"/>\n";
                if (pt.actuator >= 0 && !pt.rod) {  // barrel extent for the cylinder label
                    for (const auto& e : pt.edges)
                        for (const auto& q : e) barrel_pts[static_cast<std::size_t>(pt.actuator)].push_back(tp(q));
                }
            }
        }
        // links: convex hull of the body points offset by the rounding radius; larger links first
        std::vector<std::pair<double, std::size_t>> order;
        for (std::size_t b = 1; b < m.bodies.size(); ++b) {
            std::vector<Vec2> pts;
            for (const auto& p : bp[b]) pts.push_back(kin.point(s, static_cast<int>(b), p));
            const auto h = hull(pts);
            double area = 0, len = 0;
            for (std::size_t i = 0; i < h.size(); ++i) {
                const Vec2& a = h[i]; const Vec2& c = h[(i + 1) % h.size()];
                area += a[0] * c[1] - a[1] * c[0];
                len = std::max(len, std::hypot(c[0] - a[0], c[1] - a[1]));
            }
            order.push_back({std::fabs(area) / 2 + len * 50.0, b});
        }
        std::sort(order.begin(), order.end(), [](const auto& a, const auto& c) { return a.first > c.first; });
        for (const auto& [area_, b] : order) {
            std::vector<Vec2> pts;
            for (const auto& p : bp[b]) pts.push_back(kin.point(s, static_cast<int>(b), p));
            auto h = hull(pts);
            if (h.empty()) continue;
            std::vector<P2> hp;
            for (const auto& q : h) hp.push_back({X(q), Y(q)});
            // the rounded outline (hull offset by rr) is drawn as a wide round-joined stroke under a narrower one
            std::ostringstream d;
            for (std::size_t i = 0; i < hp.size(); ++i) d << (i ? "L" : "M") << fmt(hp[i].x) << " " << fmt(hp[i].y);
            if (hp.size() == 1) d << "h0.001";
            if (hp.size() >= 3) d << "Z";
            for (std::size_t i = 0; i < hp.size(); ++i) {
                const P2 a = hp[i], c = hp[(i + 1) % hp.size()];
                const double L = std::max(1e-9, std::hypot(c.x - a.x, c.y - a.y));
                const P2 nn{-(c.y - a.y) / L, (c.x - a.x) / L};
                for (double sg : {1.0, -1.0}) seg_obst({a.x + sg * rr * nn.x, a.y + sg * rr * nn.y}, {c.x + sg * rr * nn.x, c.y + sg * rr * nn.y}, 0.02);
            }
            if (main_pose) hulls.push_back({static_cast<int>(b), hp});
            if (labels && opt.items && main_pose && m.labels.count(static_cast<int>(b))) {
                Lbl l{m.labels.at(static_cast<int>(b)), {0, 0}, {}, {}, static_cast<int>(b)};
                lbls.push_back(l);  // candidates are generated once all links are known
            }
            static const char* link_col[] = {"#0072B2", "#009E73", "#CC79A7", "#56B4E9", "#E69F00", "#000000"};
            const bool big = m.outlines.count(static_cast<int>(b)) > 0;  // structure drawn from outlines (mast, hammer)
            const std::string ec = !colour ? "#000" : big ? "#7F7F7F" : link_col[(b - 1) % 6];
            const std::string fc = !colour ? "#fff" : big ? "#F2F2F2" : "#FFFFFF";
            if (main_pose && have_parts) {
                // drawn from the solid model
            } else if (main_pose) {
                o << "<path d=\"" << d.str() << "\" fill=\"" << ec << "\" stroke=\"" << ec << "\" stroke-width=\"" << fmt(2 * rr + lw)
                  << "\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>\n";
                o << "<path d=\"" << d.str() << "\" fill=\"" << fc << "\" stroke=\"" << fc << "\" stroke-width=\"" << fmt(2 * rr - lw)
                  << "\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>\n";
            } else {
                o << "<path d=\"" << d.str() << "\" fill=\"none\" stroke=\"" << ec << "\" stroke-width=\"" << fmt(lw) << "\" stroke-linejoin=\"round\""
                  << ds << "/>\n";
            }
        }
        // cylinders: barrel from pin a, rod to pin b
        for (std::size_t k = 0; k < m.actuators.size(); ++k) {
            const auto& a = m.actuators[k];
            if (!a.has_cylinder) continue;
            const Vec2 pa = kin.point(s, a.body_a, a.pa), pb = kin.point(s, a.body_b, a.pb);
            const Vec2 u = kin.actuator_dir(s, static_cast<int>(k));
            const double lb = std::max(a.cyl.lmin - 0.08 * a.cyl.stroke(), 0.4 * a.cyl.lmin) * 0.92;
            P2 A{X(pa), Y(pa)}, B{X(pb), Y(pb)};
            const P2 uu{u[0], -u[1]}, nn{-uu.y, uu.x};
            P2 E{A.x + uu.x * lb * sc, A.y + uu.y * lb * sc};
            double hw = std::max(0.75, a.cyl.bore * 1.3 * sc) / 2, rw = std::max(0.3, a.cyl.rod * sc) / 2;
            const bool from_model = act_drawn[k];
            if (from_model && !barrel_pts[k].empty()) {  // the barrel as drawn from the solid model
                const auto& bp2 = barrel_pts[k];
                P2 c{0, 0};
                for (const auto& q : bp2) { c.x += q.x / static_cast<double>(bp2.size()); c.y += q.y / static_cast<double>(bp2.size()); }
                double t0 = 1e9, t1 = -1e9, wmax = 0;
                for (const auto& q : bp2) {
                    const double t = (q.x - c.x) * uu.x + (q.y - c.y) * uu.y, w = std::fabs((q.x - c.x) * nn.x + (q.y - c.y) * nn.y);
                    t0 = std::min(t0, t); t1 = std::max(t1, t); wmax = std::max(wmax, w);
                }
                A = {c.x + uu.x * t0, c.y + uu.y * t0};
                E = {c.x + uu.x * t1, c.y + uu.y * t1};
                hw = wmax;
            }
            auto poly = [&](P2 c0, P2 c1, double w, const char* fill, double sw) {
                o << "<polygon points=\"";
                for (auto [c, sg] : {std::pair<P2, double>{c0, 1}, {c1, 1}, {c1, -1}, {c0, -1}})
                    o << fmt(c.x + nn.x * w * sg) << "," << fmt(c.y + nn.y * w * sg) << " ";
                o << "\" fill=\"" << fill << "\" stroke=\"#000\" stroke-width=\"" << fmt(sw) << "\" stroke-linejoin=\"round\"" << ds << "/>\n";
            };
            if (!from_model) {
                poly(E, B, rw, "#fff", lw * 0.8);
                poly(A, E, hw, colour ? "#D55E00" : "#fff", lw);
            }
            seg_obst(A, B);
            for (double sg : {1.0, -1.0}) seg_obst({A.x + nn.x * hw * sg, A.y + nn.y * hw * sg}, {E.x + nn.x * hw * sg, E.y + nn.y * hw * sg});
            if (labels && main_pose) {
                Lbl l{a.name, {(A.x + E.x) / 2, (A.y + E.y) / 2}, {}, {}};
                for (double off : {0.4, 1.4, 2.6, 3.8})
                    for (double t : {0.5, 0.3, 0.7, 0.15, 0.85})
                        for (double sg : {1.0, -1.0}) {
                            const P2 c{A.x + (E.x - A.x) * t, A.y + (E.y - A.y) * t};
                            l.cand.push_back({c.x + nn.x * sg * (hw + off), c.y + nn.y * sg * (hw + off)});
                            l.pen.push_back(0.3 * off + 0.5 * std::fabs(t - 0.5));
                        }
                lbls.push_back(l);
            }
        }
        // joints and pins
        for (const auto& j : m.revolutes) {
            const Vec2 q = kin.point(s, j.body_b, j.p);
            o << "<circle cx=\"" << fmt(X(q)) << "\" cy=\"" << fmt(Y(q)) << "\" r=\"0.65\" fill=\"#fff\" stroke=\"#000\" stroke-width=\"0.22\"/>\n";
            obst.push_back({X(q), Y(q)});
            if (labels && main_pose) {
                Lbl l{j.name, {X(q), Y(q)}, {}, {}};
                for (int k = 0; k < 12; ++k) {
                    const double ang = k * M_PI / 6;
                    l.cand.push_back({X(q) + 1.3 * std::cos(ang), Y(q) - 1.3 * std::sin(ang)});
                }
                lbls.push_back(l);
            }
        }
        for (std::size_t k = 0; k < m.actuators.size(); ++k) {
            const auto& a = m.actuators[k];
            if (!a.has_cylinder || act_drawn[k]) continue;  // a cylinder drawn from the solid model shows its own eyes
            for (const Vec2 q : {kin.point(s, a.body_a, a.pa), kin.point(s, a.body_b, a.pb)}) {
                o << "<circle cx=\"" << fmt(X(q)) << "\" cy=\"" << fmt(Y(q)) << "\" r=\"0.5\" fill=\"#fff\" stroke=\"#000\" stroke-width=\"0.18\"/>";
                obst.push_back({X(q), Y(q)});
            }
            o << "\n";
        }
        // centres of gravity of the bodies (first pose only)
        if (main_pose && opt.cg)
            for (std::size_t b = 1; b < m.bodies.size(); ++b) {
                double mm = 0, gx = 0, gy = 0;
                for (const auto& mi : m.masses)
                    if (mi.body == static_cast<int>(b)) {
                        const Vec2 q = kin.point(s, mi.body, mi.c);
                        mm += mi.m; gx += mi.m * q[0]; gy += mi.m * q[1];
                    }
                if (mm <= 0) continue;
                const double x = cx + sc * gx / mm, y = cy - sc * gy / mm, r = 0.75;
                o << "<circle cx=\"" << fmt(x) << "\" cy=\"" << fmt(y) << "\" r=\"" << fmt(r) << "\" fill=\"#fff\" stroke=\"#000\" stroke-width=\"0.18\"/>"
                  << "<path d=\"M" << fmt(x) << " " << fmt(y) << "L" << fmt(x + r) << " " << fmt(y) << "A" << fmt(r) << " " << fmt(r) << " 0 0 0 "
                  << fmt(x) << " " << fmt(y - r) << "ZM" << fmt(x) << " " << fmt(y) << "L" << fmt(x - r) << " " << fmt(y) << "A" << fmt(r) << " "
                  << fmt(r) << " 0 0 0 " << fmt(x) << " " << fmt(y + r) << "Z\" fill=\"#000\"/>\n";
                for (int k = 0; k < 8; ++k) obst.push_back({x + r * std::cos(k * M_PI / 4), y + r * std::sin(k * M_PI / 4)});
                obst.push_back({x, y});
            }
        o << "</g>\n";
    }
    // body labels: candidates on a grid inside the link (and the view), preferably not inside another link
    {
        double vx0 = 0, vx1 = W, vy0 = 0, vy1 = H;
        if (opt.view.size() == 4) { vx0 = X({xlo, 0}); vx1 = X({xhi, 0}); vy0 = Y({0, yhi}); vy1 = Y({0, ylo}); }
        for (auto& l : lbls) {
            if (l.body < 0) continue;
            const std::vector<P2>* own = nullptr;
            for (const auto& hb : hulls) if (hb.first == l.body) own = &hb.second;
            if (!own || own->empty()) continue;
            double x0 = 1e9, x1 = -1e9, y0 = 1e9, y1 = -1e9;
            for (const auto& q : *own) { x0 = std::min(x0, q.x); x1 = std::max(x1, q.x); y0 = std::min(y0, q.y); y1 = std::max(y1, q.y); }
            x0 = std::max(x0, vx0 + 2); x1 = std::min(x1, vx1 - 2); y0 = std::max(y0, vy0 + 2); y1 = std::min(y1, vy1 - 2);
            P2 g{(x0 + x1) / 2, (y0 + y1) / 2};
            l.at = g;
            const double step = body_polys.empty() ? 0.8 : 0.4;
            for (double gx = x0; gx <= x1 + 1e-9; gx += step)
                for (double gy = y0; gy <= y1 + 1e-9; gy += step) {
                    const P2 q{gx, gy};
                    if (own->size() >= 3 && !inside(*own, q)) continue;
                    if (body_polys.count(l.body)) {  // on the drawn part itself (point in polygon, even-odd)
                        bool on = false;
                        for (const auto& poly : body_polys.at(l.body)) {
                            bool in = false;
                            for (std::size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++)
                                if ((poly[i].y > q.y) != (poly[j].y > q.y) &&
                                    q.x < (poly[j].x - poly[i].x) * (q.y - poly[i].y) / (poly[j].y - poly[i].y) + poly[i].x)
                                    in = !in;
                            if (in) { on = true; break; }
                        }
                        if (!on) continue;
                    }
                    double pen = 0.01 * std::hypot(gx - g.x, gy - g.y);
                    const double bw = fs * 0.62 * static_cast<double>(l.text.size()) / 2 + 0.3, bh = fs * 0.45;
                    for (const auto& hb : hulls) {
                        if (hb.first == l.body) continue;
                        for (const P2 cq : {q, P2{q.x - bw, q.y - bh}, P2{q.x + bw, q.y - bh}, P2{q.x - bw, q.y + bh}, P2{q.x + bw, q.y + bh}})
                            if (inside(hb.second, cq)) pen += 1.0;
                    }
                    l.cand.push_back(q);
                    l.pen.push_back(pen);
                }
            if (l.cand.empty() && own->size() < 3) { l.cand.push_back(g); l.pen.push_back(0); }
        }
    }
    if (opt.view.size() == 4) {
        o << "</g>\n";
        // break marks where a link drawn from outlines leaves the view
        const double bx0 = X({xlo, 0}), bx1 = X({xhi, 0}), by0 = Y({0, yhi}), by1 = Y({0, ylo});
        std::vector<std::pair<double, double>> cuts[4];
        for (const auto& hb : hulls) {
            if (!m.outlines.count(hb.first)) continue;
            const auto& h = hb.second;
            for (int e = 0; e < 4; ++e) {
                const bool vert = e < 2;
                const double c = e == 0 ? bx0 : e == 1 ? bx1 : e == 2 ? by0 : by1;
                std::vector<double> hits;
                for (std::size_t i = 0; i < h.size(); ++i) {
                    const P2 a = h[i], d = h[(i + 1) % h.size()];
                    const double ua = vert ? a.x : a.y, ud = vert ? d.x : d.y;
                    if ((ua - c) * (ud - c) > 0 || ua == ud) continue;
                    const double t = (c - ua) / (ud - ua);
                    hits.push_back(vert ? a.y + t * (d.y - a.y) : a.x + t * (d.x - a.x));
                }
                if (hits.size() < 2) continue;
                const double lo = *std::min_element(hits.begin(), hits.end()) - rr - 0.6, hi = *std::max_element(hits.begin(), hits.end()) + rr + 0.6;
                const double lim0 = vert ? by0 : bx0, lim1 = vert ? by1 : bx1;
                const double a0 = std::max(lo, lim0), a1 = std::min(hi, lim1);
                if (a1 > a0) cuts[e].push_back({a0, a1});
            }
        }
        for (int e = 0; e < 4; ++e) {
            auto& cv = cuts[e];
            std::sort(cv.begin(), cv.end());
            std::vector<std::pair<double, double>> merged;
            for (const auto& iv : cv) {
                if (!merged.empty() && iv.first <= merged.back().second) merged.back().second = std::max(merged.back().second, iv.second);
                else merged.push_back(iv);
            }
            const bool vert = e < 2;
            const double c = e == 0 ? bx0 : e == 1 ? bx1 : e == 2 ? by0 : by1;
            for (const auto& [a0, a1] : merged) {
                // zigzag across the cut
                std::ostringstream z;
                const int n = std::max(2, static_cast<int>((a1 - a0) / 1.2));
                for (int k = 0; k <= n; ++k) {
                    const double u = a0 + (a1 - a0) * k / n, off = (k % 2 ? 0.5 : -0.5);
                    const double px = vert ? c + off : u, py = vert ? u : c + off;
                    z << (k ? "L" : "M") << fmt(px) << " " << fmt(py);
                }
                o << "<path d=\"" << z.str() << "\" fill=\"none\" stroke=\"#fff\" stroke-width=\"1.0\"/><path d=\"" << z.str()
                  << "\" fill=\"none\" stroke=\"#000\" stroke-width=\"0.25\"/>\n";
            }
        }
        // labels must stay inside the view
        const double vx0 = X({xlo, 0}), vx1 = X({xhi, 0}), vy0 = Y({0, yhi}), vy1 = Y({0, ylo});
        for (auto& l : lbls)
            if (l.at.x < vx0 || l.at.x > vx1 || l.at.y < vy0 || l.at.y > vy1) l.cand.clear();
    }
    // labels: candidate positions scored against drawn points and earlier labels
    struct Box { double x0, y0, x1, y1; };
    std::vector<Box> placed;
    std::stable_partition(lbls.begin(), lbls.end(), [](const Lbl& l) { return l.body < 0; });  // link labels last
    double lv0 = 0.5, lv1 = W - 0.5, lw0 = 0.5, lw1 = H - 0.5;  // label area
    if (opt.view.size() == 4) { lv0 = X({xlo, 0}); lv1 = X({xhi, 0}); lw0 = Y({0, yhi}); lw1 = Y({0, ylo}); }
    std::ostringstream leaders;
    for (const auto& l : lbls) {
        if (l.cand.empty()) continue;
        const double w = fs * 0.62 * static_cast<double>(l.text.size()) + 0.3, h = fs * 0.78;
        if (l.body >= 0) continue;  // link labels: see below
        double best = 1e300; Box bb{};
        for (std::size_t ci = 0; ci < l.cand.size(); ++ci) {
            const P2 c = l.cand[ci];
            // place the box so that its nearest side touches the candidate point, away from the anchor
            const double dx = c.x - l.at.x, dy = c.y - l.at.y;
            double bx = c.x - w / 2, by = c.y - h / 2;
            if (std::fabs(dx) > 1e-9) bx = dx > 0 ? c.x - w / 2 + w / 2 * std::min(1.0, std::fabs(dx) / (std::hypot(dx, dy) + 1e-12)) : c.x - w / 2 - w / 2 * std::min(1.0, std::fabs(dx) / (std::hypot(dx, dy) + 1e-12));
            if (std::fabs(dy) > 1e-9) by = dy > 0 ? c.y - h / 2 + h / 2 * std::min(1.0, std::fabs(dy) / (std::hypot(dx, dy) + 1e-12)) : c.y - h / 2 - h / 2 * std::min(1.0, std::fabs(dy) / (std::hypot(dx, dy) + 1e-12));
            Box b{bx, by, bx + w, by + h};
            double cost = ci < l.pen.size() ? l.pen[ci] : 0.01 * static_cast<double>(ci);
            for (const auto& q : obst)
                if (q.x > b.x0 - 0.25 && q.x < b.x1 + 0.25 && q.y > b.y0 - 0.25 && q.y < b.y1 + 0.25) cost += q.w;
            for (const auto& p : placed) {
                const double ox = std::min(b.x1, p.x1) - std::max(b.x0, p.x0), oy = std::min(b.y1, p.y1) - std::max(b.y0, p.y0);
                if (ox > -0.3 && oy > -0.3) cost += 20;
            }
            if (b.x0 < 0.5 || b.y0 < 0.5 || b.x1 > W - 0.5 || b.y1 > H - 0.5) cost += 50;
            if (cost < best) { best = cost; bb = b; }
        }
        placed.push_back(bb);
        lab << "<text x=\"" << fmt((bb.x0 + bb.x1) / 2) << "\" y=\"" << fmt(bb.y1 - fs * 0.08) << "\" text-anchor=\"middle\" stroke=\"#fff\" stroke-width=\"0.5\" "
            << "stroke-linejoin=\"round\" paint-order=\"stroke\">" << esc(l.text) << "</text>\n";
    }
    // link labels as item numbers: a column beside the linkage with leaders, ordered by height so that leaders do not cross
    {
        struct Item { std::string text; P2 anc; double w; };
        std::vector<Item> items;
        double mx0 = 1e9, mx1 = -1e9;  // extent of everything drawn, clipped to the view
        for (const auto& hb : hulls)
            for (const auto& q : hb.second) { mx0 = std::min(mx0, std::max(q.x, lv0)); mx1 = std::max(mx1, std::min(q.x, lv1)); }
        for (const auto& q : obst) { mx0 = std::min(mx0, std::max(q.x, lv0)); mx1 = std::max(mx1, std::min(q.x, lv1)); }
        const double mid = (mx0 + mx1) / 2;
        for (const auto& l : lbls) {
            if (l.body < 0 || l.cand.empty()) continue;
            P2 anc = l.cand[0];
            double ab = 1e300;
            for (std::size_t ci = 0; ci < l.cand.size(); ++ci) {
                const P2 c = l.cand[ci];
                double cost = l.pen[ci];
                for (const auto& q : obst) {
                    const double d = std::hypot(q.x - c.x, q.y - c.y);
                    if (q.w >= 1 && d < 1.5) cost += 1;
                    else if (q.w < 1 && d < 1.5) cost += 0.3;  // well inside the link, away from its outline
                }
                for (const auto& pb : placed)  // clear of the joint and cylinder labels
                    if (c.x > pb.x0 - 0.8 && c.x < pb.x1 + 0.8 && c.y > pb.y0 - 0.8 && c.y < pb.y1 + 0.8) cost += 5;
                // the leader runs about horizontally to the nearer side: keep it off labels and pins
                const double xe = c.x < mid ? mx0 : mx1;
                for (double x = std::min(c.x, xe); x <= std::max(c.x, xe); x += 0.5) {
                    for (const auto& pb : placed)
                        if (x > pb.x0 - 0.3 && x < pb.x1 + 0.3 && c.y > pb.y0 - 0.3 && c.y < pb.y1 + 0.3) cost += 2;
                    for (const auto& q : obst)
                        if (q.w >= 1 && std::fabs(q.x - x) < 0.4 && std::fabs(q.y - c.y) < 0.4) cost += 0.3;
                }
                if (cost < ab) { ab = cost; anc = c; }
            }
            items.push_back({l.text, anc, fs * 0.62 * static_cast<double>(l.text.size()) + 0.3});
        }
        if (!items.empty()) {
            // two columns, left and right of everything drawn (clipped to the view); each item goes to its nearer side
            const double gap = 3.0, hh = fs * 0.78, dy_min = hh + 1.0;
            double lmin = 1e9, rmin = 1e9;  // free space beside the drawing
            double maxw = 0;
            for (const auto& it : items) maxw = std::max(maxw, it.w);
            lmin = mx0 - lv0; rmin = lv1 - mx1;
            std::sort(items.begin(), items.end(), [](const Item& a, const Item& c) { return a.anc.y < c.anc.y; });
            for (int side : {-1, 1}) {
            std::vector<Item> col;
            for (const auto& it : items) {
                bool left = it.anc.x < mid;
                if (left && lmin < maxw + gap) left = false;
                if (!left && rmin < maxw + gap && lmin >= maxw + gap) left = true;
                if ((side < 0) == left) col.push_back(it);
            }
            if (col.empty()) continue;
            const double colx = side > 0 ? mx1 + gap : mx0 - gap;
            std::vector<double> ys;
            for (const auto& it : col) ys.push_back(it.anc.y);
            for (std::size_t i = 1; i < ys.size(); ++i) ys[i] = std::max(ys[i], ys[i - 1] + dy_min);
            for (std::size_t i = ys.size(); i-- > 1;) ys[i - 1] = std::min(ys[i - 1], ys[i] - dy_min);  // keep the column compact
            for (std::size_t i = 0; i < col.size(); ++i) {
                const auto& it = col[i];
                const double y = std::clamp(ys[i], lw0 + hh, lw1 - hh);
                const double xe = colx - side * 0.6;              // end of the leader at the text
                const double xk = colx - side * 2.2;              // knee: the last part of the leader is horizontal
                leaders << "<polyline points=\"" << fmt(it.anc.x) << "," << fmt(it.anc.y) << " " << fmt(xk) << "," << fmt(y) << " " << fmt(xe) << ","
                        << fmt(y) << "\" fill=\"none\" stroke=\"#000\" stroke-width=\"0.15\"/><circle cx=\"" << fmt(it.anc.x) << "\" cy=\""
                        << fmt(it.anc.y) << "\" r=\"0.4\" fill=\"#000\"/>\n";
                lab << "<text x=\"" << fmt(colx) << "\" y=\"" << fmt(y + fs * 0.35) << "\" text-anchor=\"" << (side > 0 ? "start" : "end")
                    << "\" stroke=\"#fff\" stroke-width=\"0.5\" stroke-linejoin=\"round\" paint-order=\"stroke\">" << esc(it.text) << "</text>\n";
            }
            }
        }
    }
    o << leaders.str() << lab.str();
    if (!opt.title.empty())
        o << "<text x=\"" << fmt(W / 2) << "\" y=\"" << fmt(ml * 0.5 + fs * 1.1) << "\" text-anchor=\"middle\" font-size=\"" << fmt(fs * 1.15)
          << "\">" << esc(opt.title) << "</text>\n";
    // scale bar: 1 m
    {
        const double L = 1000 * sc, x0 = W - ml - L, y0 = H - 2.2;
        o << "<path d=\"M" << fmt(x0) << " " << fmt(y0 - 0.8) << "V" << fmt(y0) << "H" << fmt(x0 + L) << "V" << fmt(y0 - 0.8)
          << "\" fill=\"none\" stroke=\"#000\" stroke-width=\"0.25\"/><text x=\"" << fmt(x0 + L / 2) << "\" y=\"" << fmt(y0 - 0.6)
          << "\" text-anchor=\"middle\">1 m</text>\n";
    }
    o << "</svg>\n";
    return o.str();
}

}  // namespace cylload
