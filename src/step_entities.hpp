// SPDX-License-Identifier: MIT
// Internal: STEP Part 21 entity table and small parsing helpers shared by step.cpp and brep.cpp.
#pragma once
#include <array>
#include <string>
#include <unordered_map>
#include <vector>

namespace cylload {
namespace step_detail {

struct Entity {
    std::string type;  ///< empty for complex entities "( A() B() )"
    std::string args;  ///< text between the outer parentheses
};
using Entities = std::unordered_map<long, Entity>;
using Mat4 = std::array<double, 16>;  // row-major rigid transform

Entities parse(const std::string& path);
const Entity* get(const Entities& e, long id);
std::vector<long> refs(const std::string& s);
std::vector<double> numbers(const std::string& s);
std::vector<std::string> strings(const std::string& s);
Mat4 identity();
Mat4 mul(const Mat4& a, const Mat4& b);
Mat4 rigid_inverse(const Mat4& m);
Mat4 placement(const Entities& e, long id);

/// Volume integrals of one solid in its own coordinates: V, first moments about `ref`
/// and second moments about `ref` (xx, yy, zz, xy, yz, zx), all in mm^k.
struct VolumeIntegrals {
    std::array<double, 3> ref{};
    /// [10], [11]: volume from the y and z fields; [12..14]: xx, yy, zz from other fields (checks)
    double v[15] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::size_t faces = 0, faces_triangulated = 0, solids = 0;
};
/// Integrate all solids (MANIFOLD_SOLID_BREP, BREP_WITH_VOIDS) found in the given
/// shape representations; throws on unsupported curve types.
VolumeIntegrals integrate_solids(const Entities& e, const std::vector<long>& reps);

/// Drawing geometry of one solid in its own coordinates: the boundary loops of every face
/// and every edge, as polylines sampled along the curves.
using Polyline3 = std::vector<std::array<double, 3>>;
struct SolidGeometry {
    std::vector<std::vector<Polyline3>> faces;  ///< loops of each face
    std::vector<char> curved;                   ///< face lies on a curved surface
    std::vector<Polyline3> edges;               ///< each edge once
};
std::vector<SolidGeometry> solid_geometry(const Entities& e, const std::vector<long>& reps);

}  // namespace step_detail
}  // namespace cylload
