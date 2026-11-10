// SPDX-License-Identifier: MIT
#pragma once
#include <array>
#include <string>
#include <vector>

namespace cylload {

/// One occurrence of a part or sub-assembly in a STEP (AP203/AP214) assembly tree,
/// with the axis-aligned bounding box of its B-rep points in the root coordinate system.
/// Drawing geometry of one solid in root coordinates (with read_step(..., geometry = true)).
struct StepSolid {
    std::vector<std::vector<std::vector<std::array<double, 3>>>> faces;  ///< sampled boundary loops of each face
    std::vector<char> curved;                                            ///< the face lies on a curved surface
    std::vector<std::vector<std::array<double, 3>>> edges;               ///< sampled edges
};

struct StepOccurrence {
    std::vector<std::string> path;  ///< product names from the root to this occurrence
    long product_definition = 0;
    std::size_t n_points = 0;
    std::array<double, 3> lo{}, hi{};
    std::array<double, 3> centre() const { return {(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2}; }
    // mass properties of the solids (only with read_step(..., mass_properties = true))
    double volume = 0;                 ///< [mm^3]
    std::array<double, 3> centroid{};  ///< volume centroid, root coordinates [mm]
    std::array<double, 9> cov{};       ///< int (r-c)(r-c)^T dV, root coordinates, row-major [mm^5]
    std::size_t faces = 0, faces_triangulated = 0, faces_parametric = 0;  ///< all, triangle fans, NURBS in (u, v)
    double area = 0, area_triangulated = 0;  ///< surface area of all faces and of the triangulated faces [mm^2]
    double volume_check = 0;  ///< max relative difference of the volume from the x, y and z fields
    double moment_check = 0;  ///< same for the second moments (two fields each), relative to the trace
    std::vector<StepSolid> solids;  ///< drawing geometry (only with read_step(..., geometry = true))
};

/// Mass properties of a set of occurrences: volume, centroid and second moments.
struct MassProps {
    double volume = 0;
    std::array<double, 3> centroid{};
    std::array<double, 9> cov{};
    /// squared radius of gyration about the axis through the centroid normal to the
    /// plane spanned by root axes (i, j); J = m * k2 [m^2]
    double k2(int i = 0, int j = 1) const { return volume > 0 ? (cov[4 * i] + cov[4 * j]) / volume * 1e-6 : 0.0; }
    void add(const MassProps& o);
};

struct StepGroup {
    std::string key;
    std::size_t n_occurrences = 0, n_points = 0;
    std::array<double, 3> lo{}, hi{};
    std::array<double, 3> centre() const { return {(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2}; }
    MassProps mass;
    std::size_t faces = 0, faces_triangulated = 0, faces_parametric = 0;
    double area = 0, area_triangulated = 0;
    double volume_check = 0, moment_check = 0;
};

struct StepModel {
    std::string root;
    std::size_t n_entities = 0;
    std::vector<StepOccurrence> occurrences;  ///< only occurrences with geometry
    std::size_t n_tree_nodes = 0;             ///< all occurrences in the tree
    double length_unit_mm = 1.0;              ///< file length unit; all results are converted to mm
    std::string length_unit = "mm";           ///< as read from the file ("mm (assumed)" if none is declared)
};

/// Read a STEP file; `root` selects the top assembly by product name (empty = the
/// assembly that is not used inside any other). With `mass_properties`, the volume,
/// centroid and second moments of every solid are integrated from its B-rep.
/// With `geometry`, the sampled faces and edges of every solid are kept for drawing.
StepModel read_step(const std::string& path, const std::string& root = "", bool mass_properties = false, bool geometry = false);
/// Union of bounding boxes grouped by the path element at `level` (0 = root;
/// shallower occurrences use their last element).
std::vector<StepGroup> group_by_level(const StepModel& m, std::size_t level);
/// Union of all groups whose key starts with `prefix` (throws if none).
StepGroup union_by_prefix(const std::vector<StepGroup>& groups, const std::string& prefix);
/// B-spline surface faces are integrated in their parameter domain (default, `true`) or, as in
/// versions up to 1.1, replaced by fans of triangles spanning their boundary (`false`).
void set_bspline_integration(bool parametric);
/// Decode STEP string escapes (\X2\...\X0\, \X\hh, '') to UTF-8.
std::string decode_step_string(const std::string& s);

}  // namespace cylload
