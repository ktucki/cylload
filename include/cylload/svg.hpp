// SPDX-License-Identifier: MIT
// Publication figures as SVG without dependencies: black by default (optional Okabe-Ito
// palette and grid), thin lines, axes that end
// on a labelled tick, panel letters below the panels, legend below or inside.
#pragma once
#include <limits>
#include <string>
#include <vector>

#include "cylload/model.hpp"

namespace cylload {

enum class LineStyle { Solid, Dashed, Dotted, DashDot, None };
enum class Marker { None, Circle, CircleOpen, Square, SquareOpen, Triangle, TriangleOpen, Diamond, DiamondOpen, Cross };

struct Series {
    std::vector<double> x, y;
    std::string label;
    LineStyle line = LineStyle::Solid;
    Marker marker = Marker::None;
    int marker_every = 1;
    double width = 0.2;  ///< line width [mm]
    std::string color;   ///< empty: black, or the figure palette when it is set
};

/// Shaded band between two curves (e.g. the 5th and 95th percentiles).
struct Band {
    std::vector<double> x, lo, hi;
    std::string label;
    std::string color = "#0072B2";
    double opacity = 0.18;
};

struct RefLine {
    bool horizontal = true;
    double value = 0;
    std::string label;
    LineStyle line = LineStyle::Dashed;
    std::string color = "#000";
};

struct Axis {
    std::string label;
    bool log = false;
    bool fixed = false;          ///< use min/max/step as given
    double min = 0, max = 0, step = 0;
};

struct Panel {
    Axis x, y;
    std::vector<Series> series;
    std::vector<RefLine> refs;
    std::vector<Band> bands;
    std::string letter;          ///< "(a)", printed below the panel
    std::string title;           ///< optional, above the plot area
    std::string legend = "below";  ///< below | inside | none
    int legend_columns = 0;        ///< entries per legend row below the plot (0: as many as fit)
    int legend_rows = 0;           ///< rows reserved for the legend below the plot (aligns panels of different legends)
};

struct Figure {
    double width_mm = 160, height_mm = 70;
    int columns = 1;
    double font_pt = 8;
    std::string font = "Arial, Helvetica, sans-serif";
    std::vector<Panel> panels;
    std::vector<std::string> palette;  ///< colours given to series without an explicit colour
    bool grid = false;                 ///< light grid at the major ticks
};
/// Okabe-Ito palette, distinguishable with colour-vision deficiencies.
std::vector<std::string> okabe_ito();

std::string render_svg(const Figure& f);
/// Parse a plot specification (see docs/input_format.md); CSV paths are relative to base_dir.
Figure figure_from_spec(const std::string& text, const std::string& base_dir = ".");

/// "Nice" axis: limits are multiples of the step, so the axis ends on a labelled tick.
void nice_axis(double lo, double hi, int target_ticks, double& amin, double& amax, double& step);

/// A part of the solid model to draw with the mechanism, projected on the machine plane in the
/// reference configuration. It moves with its link, or with the barrel or rod of a cylinder.
struct DrawPart {
    int body = 0;            ///< link (0 = ground)
    int actuator = -1;       ///< >= 0: moves with this cylinder instead of a link
    bool rod = false;        ///< with `actuator`: the rod (moves with pin b), else the barrel (pin a)
    double depth = 0;        ///< distance towards the viewer; nearer parts are drawn later
    std::vector<std::vector<std::vector<Vec2>>> faces;  ///< boundary loops of each face
    std::vector<char> curved;  ///< curved faces are filled as the hull of their loops (a tube seen side-on has two edge loops)
    std::vector<std::vector<Vec2>> edges;
};

/// Options of a mechanism drawing.
struct DrawOptions {
    double width_mm = 120, height_mm = 100;  ///< size when scale_mm_per_m is 0 (height is a minimum otherwise)
    double scale_mm_per_m = 0;               ///< fixed scale; the width follows from the mechanism
    double y_bottom = std::numeric_limits<double>::quiet_NaN();  ///< model ordinate at the lower margin
    bool labels = true;                      ///< joint, cylinder and body labels
    bool colour = false;                     ///< tinted links and cylinder barrels
    bool cg = false;                         ///< centre-of-gravity symbol of each link
    bool items = true;                       ///< link labels (`label` in the model) as item numbers with leaders
    double font_pt = 7;
    std::string title;                       ///< optional title above the drawing
    std::vector<double> view;                ///< xmin, ymin, xmax, ymax (model mm): drawn region, the rest is cut off
    std::vector<DrawPart> parts;             ///< solid-model parts drawn instead of the link outlines
};
/// Drawing of a mechanism in one or more poses (first pose solid, the others dashed). A link is the rounded
/// convex hull of its joints, cylinder pins and `outline` points (its centres of mass when it has fewer than two
/// such points); larger links are drawn first. Labels are placed where they cover the fewest marks, and a 1 m
/// scale bar is drawn. With a fixed scale and y_bottom several drawings share one scale and one ground level.
std::string draw_mechanism(const Kinematics& kin, const std::vector<Pose>& poses, const DrawOptions& opt);
/// Shorthand with a fitted scale.
std::string draw_mechanism(const Kinematics& kin, const std::vector<Pose>& poses, double width_mm, double height_mm,
                           bool labels = true, double font_pt = 7, bool colour = false);

}  // namespace cylload
