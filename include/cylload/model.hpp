// SPDX-License-Identifier: MIT
#pragma once
#include <array>
#include <map>
#include <string>
#include <vector>

#include "cylload/linalg.hpp"

namespace cylload {

using Vec2 = std::array<double, 2>;

/// Hydraulic cylinder data attached to an actuator (all lengths in mm).
struct Cylinder {
    double bore = 0.0;     ///< piston diameter D [mm]
    double rod = 0.0;      ///< rod diameter d [mm]
    double lmin = 0.0;     ///< pin-to-pin length, fully retracted [mm]
    double lmax = 0.0;     ///< pin-to-pin length, fully extended [mm]
    double mass = 0.0;     ///< cylinder mass [kg]; half is lumped at each pin
    double p_back = 0.0;   ///< back-pressure in the unloaded chamber [bar] (from `hydraulics back_pressure=`)
    double area_piston() const;  ///< [mm^2]
    double area_annulus() const; ///< [mm^2]
    double stroke() const { return lmax - lmin; }
};

/// Two-point force element whose pin-to-pin length is a driving coordinate.
/// Without cylinder data it is a generic actuator (winch, screw jack, ...).
struct Actuator {
    std::string name;
    int body_a = 0, body_b = 0;  ///< body indices (0 = ground)
    Vec2 pa{}, pb{};             ///< pin positions in the reference configuration [mm]
    double l_ref = 0.0;          ///< pin distance in the reference configuration [mm]
    bool has_cylinder = false;
    Cylinder cyl;
};

/// Revolute joint between two bodies at a point of the reference configuration.
struct Revolute {
    std::string name;
    int body_a = 0, body_b = 0;
    Vec2 p{};
};

/// Prismatic joint: body_b slides relative to body_a along direction u (reference
/// configuration) through point p; relative rotation is locked.
struct Prismatic {
    std::string name;
    int body_a = 0, body_b = 0;
    Vec2 p{}, u{};
};

/// Lumped mass with its centre of gravity in the reference configuration.
struct MassItem {
    std::string name;
    std::string group;   ///< free label used for grouping and Monte Carlo
    int body = 0;
    Vec2 c{};            ///< centre of gravity [mm]
    double m = 0.0;      ///< [kg]
    double j = 0.0;      ///< own moment of inertia about the CoG, axis normal to the plane [kg m^2]
};

/// Constant external force [N] applied at a point of a body (reference coordinates) with a fixed
/// direction in the global frame, e.g. a hitch load or a digging force; it acts like a weight.
struct ExternalForce {
    std::string name;
    int body = 0;
    Vec2 p{};  ///< point of application [mm]
    Vec2 f{};  ///< force [N]
};

struct Hydraulics {
    double beta = 1.0e9;        ///< effective bulk modulus of oil and hoses [Pa]
    double dead_volume = 0.3;   ///< dead volume on the piston side [L] (`dead_volume=` sets both sides)
    double dead_volume_rod = 0.3;  ///< dead volume on the rod side [L]
    double back_pressure = 0.0;    ///< pressure in the unloaded chamber [bar], e.g. from a counterbalance valve
    bool pin_masses_in_reduced_mass = true;  ///< `reduced_mass_pins=no` leaves the lumped cylinder masses out of m_r
};

/// Planar multibody model. Every body pose is (x, y, theta) relative to the
/// reference configuration, so local coordinates equal reference global coordinates.
class Model {
public:
    std::string title;
    Vec2 gravity{0.0, -9.81};                 ///< [m/s^2]
    std::vector<std::string> bodies{"ground"}; ///< index 0 is ground
    std::vector<Revolute> revolutes;
    std::vector<Prismatic> prismatics;
    std::vector<Actuator> actuators;
    std::vector<MassItem> masses;
    std::vector<ExternalForce> forces;          ///< constant external forces
    Hydraulics hydraulics;
    std::map<int, std::vector<Vec2>> outlines;  ///< drawing only: extra outline points of a body
    std::map<int, std::string> labels;          ///< drawing only: label of a body

    static Model from_file(const std::string& path);
    static Model from_string(const std::string& text, const std::string& base_dir = ".");
    void load_masses_csv(const std::string& path);

    int body_index(const std::string& name) const;   ///< throws if unknown
    int actuator_index(const std::string& name) const;
    int n_moving() const { return static_cast<int>(bodies.size()) - 1; }
    int n_coords() const { return 3 * n_moving(); }
    int n_constraints() const;
    /// number of degrees of freedom left after joints, must equal the number of actuators
    int dof() const { return n_coords() - 2 * static_cast<int>(revolutes.size() + prismatics.size()); }
    /// mass items including the lumped cylinder halves
    std::vector<MassItem> all_masses() const;
    double total_mass() const;
    void validate() const;

private:
    std::vector<MassItem> cylinder_masses() const;
};

/// Configuration q = (x, y, theta) for each moving body, lengths in mm, angles in rad.
struct Pose {
    std::vector<double> q;
    std::vector<double> lengths;  ///< actuator lengths [mm]
};

/// Kinematic solver for a Model.
class Kinematics {
public:
    explicit Kinematics(const Model& m);
    const Model& model() const { return m_; }

    Pose reference() const;
    /// Solve positions for the given actuator lengths, starting from `start`.
    /// Large changes are followed by continuation (steps <= max_step_mm, halved on failure)
    /// so the assembly branch of the start pose is kept.
    Pose solve(const Pose& start, const std::vector<double>& lengths, double max_step_mm = 50.0) const;

    Vec2 point(const Pose& s, int body, const Vec2& p_ref) const;
    double angle(const Pose& s, int body) const;
    /// constraint residuals and Jacobian d(Phi)/dq (rows: joints then actuators)
    std::vector<double> residual(const Pose& s) const;
    Matrix jacobian(const Pose& s) const;
    /// dq/dL for each actuator (column k = sensitivity to actuator k), other actuators locked
    Matrix sensitivity(const Pose& s) const;
    /// velocity of a material point per unit actuator length change, dr/dL_k (dimensionless)
    Vec2 point_rate(const Pose& s, const Matrix& dq, int k, int body, const Vec2& p_ref) const;
    /// actuator unit vector from pin a to pin b in the current pose
    Vec2 actuator_dir(const Pose& s, int k) const;

private:
    const Model& m_;
    bool newton(Pose& s, int max_iter = 50, double tol = 1e-11) const;
};

std::map<std::string, double> parse_assignments(const std::vector<std::string>& args);

}  // namespace cylload
