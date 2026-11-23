// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2026 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Salman Husain
// =============================================================================

#ifndef CH_HYDRAULIC_CYLINDER_MODEL_H
#define CH_HYDRAULIC_CYLINDER_MODEL_H

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "ChActuatorModel.h"

namespace chrono {
namespace actuator_models {

/// Double-acting hydraulic cylinder driven by an ideal critical-center 4/3 proportional valve (linear actuator).
///
/// The command is the normalized spool position u in [-1, 1]; u > 0 connects supply to the cap chamber (A) and the rod
/// chamber (B) to tank, extending the cylinder. Each metering edge passes Q = u K s(dp), where
///   s(dp) = dp / (dp^2 + p_t^2)^(1/4),   K = nominal_flow / sqrt(nominal_pressure_drop),
/// is a turbulent square-root law regularized to a laminar slope below the transition pressure p_t.
/// Chamber pressures follow
///   p_A' = beta / V_A (Q_A - A_A v - Q_leak),   p_B' = beta / V_B (Q_B + A_B v + Q_leak),   Q_leak = C (p_A - p_B),
/// with V_A = V_A0 + A_A x and V_B = V_B0 + A_B (stroke - x), x = initial_position + position (clamped to the stroke for
/// volumes; end stops are not modeled). Ideal relief and anti-cavitation valves bound each pressure to
/// [tank_pressure, relief_pressure]: outward derivatives are blocked at a bound and effort uses bounded pressures.
/// Output effort is p_A A_A - p_B A_B minus friction. Disabling centers the spool (u = 0).
///
/// States: [p_A, p_B] [Pa]. All evaluations are pure.
class ChHydraulicCylinderModel : public ChActuatorModel {
  public:
    enum StateIndex { P_A = 0, P_B = 1 };

    /// SI parameters (gauge pressures). Defaults: 50/28 mm bore/rod, 0.5 m stroke, 200 bar supply, 30 L/min valve.
    struct Parameters {
        double cap_area = 1.9635e-3;             ///< Piston area on the cap side A_A [m^2], positive.
        double annulus_area = 1.3477e-3;         ///< Annulus area on the rod side A_B [m^2], positive.
        double stroke = 0.5;                     ///< Stroke [m], positive.
        double initial_position = 0.25;          ///< Piston extension at the reference configuration, in [0, stroke].
        double dead_volume_cap = 5e-5;           ///< Cap-side dead volume incl. hose [m^3], positive.
        double dead_volume_rod = 5e-5;           ///< Rod-side dead volume incl. hose [m^3], positive.
        double bulk_modulus = 1e9;               ///< Effective bulk modulus beta [Pa], positive.
        double supply_pressure = 20e6;           ///< Supply pressure [Pa], greater than tank_pressure.
        double tank_pressure = 0;                ///< Tank pressure [Pa], also the anti-cavitation lower bound.
        double relief_pressure = 25e6;           ///< Relief upper bound [Pa], not less than supply_pressure.
        double initial_pressure = 0;             ///< Initial pressure of both chambers [Pa], within the bounds.
        double nominal_flow = 5e-4;              ///< Flow per metering edge at full spool and nominal drop [m^3/s].
        double nominal_pressure_drop = 3.5e6;    ///< Pressure drop for nominal_flow [Pa], positive.
        double transition_pressure = 1e5;        ///< Laminar-turbulent transition pressure p_t [Pa], positive.
        double leakage_coefficient = 0;          ///< Internal piston leakage C [m^3/(s*Pa)], nonnegative.
        double coulomb_friction = 0;             ///< Seal Coulomb friction magnitude [N], nonnegative.
        double viscous_friction = 0;             ///< Viscous coefficient [N*s/m], nonnegative.
        double friction_smoothing_speed = 1e-3;  ///< Positive Coulomb smoothing speed [m/s].
        double moving_mass = 0;                  ///< Piston and rod mass not modeled as a body [kg], nonnegative.
    };

    ChHydraulicCylinderModel() : ChHydraulicCylinderModel(Parameters{}) {}

    /// Copy and validate parameters; throw std::invalid_argument for invalid values.
    explicit ChHydraulicCylinderModel(const Parameters& parameters) : m_p(parameters) { ValidateParameters(); }

    /// Quasi-static (incompressible, leak-free) effort: chamber pressures that pass exactly the displaced flows,
    /// bounded to [tank, relief]. With a centered spool the chambers are locked (moving loads meet relief and tank
    /// pressures); at rest with a centered spool the indeterminate pressures are taken as the supply-tank mean.
    double ComputeEffort(const ActuatorCommand& command, const ActuatorState& state) const override {
        double pa, pb;
        QuasiStaticPressures(command, state.velocity, pa, pb);
        return EffortFromPressures(pa, pb, state.velocity);
    }

    std::size_t GetNumStates() const override { return 2; }

    void SetInitialConditions(std::vector<double>& y0) const override { y0.assign(2, m_p.initial_pressure); }

    void CalculateRHS(double, const std::vector<double>& y, const ActuatorCommand& command, const ActuatorState& state, std::vector<double>& rhs) const override {
        CheckStates(y);
        double pa = BoundPressure(y[P_A]), pb = BoundPressure(y[P_B]);
        Flows f = ValveFlows(Spool(command), pa, pb);
        double leak = m_p.leakage_coefficient * (pa - pb);
        double va, vb;
        Volumes(state.position, va, vb);
        rhs.assign(2, 0);
        rhs[P_A] = Blocked(y[P_A], m_p.bulk_modulus / va * (f.a - m_p.cap_area * state.velocity - leak));
        rhs[P_B] = Blocked(y[P_B], m_p.bulk_modulus / vb * (f.b + m_p.annulus_area * state.velocity + leak));
    }

    /// Analytic Jacobian. Finite differences with an absolute perturbation are unusable on pressure states (~1e7 Pa).
    bool CalculateJacobian(double, const std::vector<double>& y, const ActuatorCommand& command, const ActuatorState& state, std::vector<std::vector<double>>& jacobian) const override {
        CheckStates(y);
        double pa = BoundPressure(y[P_A]), pb = BoundPressure(y[P_B]);
        double u = Spool(command);
        Flows f = ValveFlows(u, pa, pb);
        double leak = m_p.leakage_coefficient * (pa - pb);
        double va, vb;
        Volumes(state.position, va, vb);
        double ka = m_p.bulk_modulus / va, kb = m_p.bulk_modulus / vb, cl = m_p.leakage_coefficient;
        double dfa, dfb;  // d(f.a)/dpa, d(f.b)/dpb
        if (u >= 0) {
            dfa = -u * OrificeFlowDerivative(m_p.supply_pressure - pa);
            dfb = -u * OrificeFlowDerivative(pb - m_p.tank_pressure);
        } else {
            dfa = u * OrificeFlowDerivative(pa - m_p.tank_pressure);
            dfb = u * OrificeFlowDerivative(m_p.supply_pressure - pb);
        }
        // Columns use the interior slope at the bounded pressure, also on or past a bound, so the restoring flow stays
        // implicit (a zero column would turn the next step explicit and let large steps overshoot between the bounds).
        // Rows blocked at a bound are frozen.
        double raw_a = f.a - m_p.cap_area * state.velocity - leak, raw_b = f.b + m_p.annulus_area * state.velocity + leak;
        double ra = Blocked(y[P_A], raw_a) == raw_a ? 1 : 0, rb = Blocked(y[P_B], raw_b) == raw_b ? 1 : 0;
        jacobian = {{ra * ka * (dfa - cl), ra * ka * cl}, {rb * kb * cl, rb * kb * (dfb - cl)}};
        return true;
    }

    /// Oil compressibility makes the pressure dynamics stiff.
    bool IsStiff() const override { return true; }

    double GetReflectedInertia() const override { return m_p.moving_mass; }

    double EffortFromStates(const std::vector<double>& y, const ActuatorCommand&, const ActuatorState& state) const override {
        CheckStates(y);
        return EffortFromPressures(BoundPressure(y[P_A]), BoundPressure(y[P_B]), state.velocity);
    }

    /// Quasi-static telemetry; electrical_power reports hydraulic input power p_s Q_s - p_T Q_T.
    ActuatorTelemetry ComputeTelemetry(const ActuatorCommand& command, const ActuatorState& state) const override {
        double pa, pb;
        QuasiStaticPressures(command, state.velocity, pa, pb);
        return Telemetry(pa, pb, command, state);
    }

    ActuatorTelemetry ComputeTelemetryFromStates(const std::vector<double>& y, const ActuatorCommand& command, const ActuatorState& state) const override {
        CheckStates(y);
        return Telemetry(BoundPressure(y[P_A]), BoundPressure(y[P_B]), command, state);
    }

    /// Stall force at full extension command: p_s A_A - p_T A_B [N].
    double GetMaxExtendForce() const { return m_p.supply_pressure * m_p.cap_area - m_p.tank_pressure * m_p.annulus_area; }

    /// Stall force at full retraction command: p_T A_A - p_s A_B [N] (negative).
    double GetMaxRetractForce() const { return m_p.tank_pressure * m_p.cap_area - m_p.supply_pressure * m_p.annulus_area; }

    /// Locked-chamber hydraulic stiffness beta (A_A^2 / V_A + A_B^2 / V_B) at the given position [N/m].
    double GetHydraulicStiffness(double position) const {
        double va, vb;
        Volumes(position, va, vb);
        return m_p.bulk_modulus * (m_p.cap_area * m_p.cap_area / va + m_p.annulus_area * m_p.annulus_area / vb);
    }

    /// Regularized metering-edge flow at full spool opening for pressure drop dp [m^3/s].
    double OrificeFlow(double dp) const { return FlowGain() * dp / std::pow(dp * dp + m_p.transition_pressure * m_p.transition_pressure, 0.25); }

    const Parameters& GetParameters() const { return m_p; }

  private:
    struct Flows {
        double a = 0, b = 0;            // into chambers A and B
        double supply = 0, tank = 0;    // from supply, to tank
    };

    void ValidateParameters() const {
        auto positive = [](double x) { return std::isfinite(x) && x > 0; };
        auto nonnegative = [](double x) { return std::isfinite(x) && x >= 0; };
        if (!positive(m_p.cap_area) || !positive(m_p.annulus_area) || !positive(m_p.stroke) || !std::isfinite(m_p.initial_position) ||
            m_p.initial_position < 0 || m_p.initial_position > m_p.stroke || !positive(m_p.dead_volume_cap) || !positive(m_p.dead_volume_rod) ||
            !positive(m_p.bulk_modulus) || !positive(m_p.nominal_flow) || !positive(m_p.nominal_pressure_drop) || !positive(m_p.transition_pressure) ||
            !nonnegative(m_p.leakage_coefficient) || !nonnegative(m_p.coulomb_friction) || !nonnegative(m_p.viscous_friction) ||
            !positive(m_p.friction_smoothing_speed) || !nonnegative(m_p.moving_mass))
            throw std::invalid_argument("Invalid hydraulic cylinder model parameters");
        if (!std::isfinite(m_p.tank_pressure) || !std::isfinite(m_p.relief_pressure) || !(m_p.tank_pressure < m_p.supply_pressure) ||
            !(m_p.supply_pressure <= m_p.relief_pressure) || !(m_p.initial_pressure >= m_p.tank_pressure && m_p.initial_pressure <= m_p.relief_pressure))
            throw std::invalid_argument("Hydraulic pressures must satisfy tank < supply <= relief and tank <= initial <= relief");
    }

    void CheckStates(const std::vector<double>& y) const {
        if (y.size() != 2)
            throw std::invalid_argument("Hydraulic cylinder model requires 2 states");
    }

    // d OrificeFlow / d dp = K (dp^2/2 + p_t^2) / (dp^2 + p_t^2)^(5/4)
    double OrificeFlowDerivative(double dp) const {
        double s = dp * dp + m_p.transition_pressure * m_p.transition_pressure;
        return FlowGain() * (0.5 * dp * dp + m_p.transition_pressure * m_p.transition_pressure) / std::pow(s, 1.25);
    }

    double FlowGain() const { return m_p.nominal_flow / std::sqrt(m_p.nominal_pressure_drop); }

    double Spool(const ActuatorCommand& command) const { return command.enabled ? std::clamp(command.value, -1.0, 1.0) : 0; }

    double BoundPressure(double p) const { return std::clamp(p, m_p.tank_pressure, m_p.relief_pressure); }

    double Blocked(double p, double rate) const {
        if ((p >= m_p.relief_pressure && rate > 0) || (p <= m_p.tank_pressure && rate < 0))
            return 0;
        return rate;
    }

    void Volumes(double position, double& va, double& vb) const {
        double x = std::clamp(m_p.initial_position + position, 0.0, m_p.stroke);
        va = m_p.dead_volume_cap + m_p.cap_area * x;
        vb = m_p.dead_volume_rod + m_p.annulus_area * (m_p.stroke - x);
    }

    Flows ValveFlows(double u, double pa, double pb) const {
        Flows f;
        if (u >= 0) {
            f.supply = u * OrificeFlow(m_p.supply_pressure - pa);
            f.tank = u * OrificeFlow(pb - m_p.tank_pressure);
            f.a = f.supply;
            f.b = -f.tank;
        } else {
            f.supply = -u * OrificeFlow(m_p.supply_pressure - pb);
            f.tank = -u * OrificeFlow(pa - m_p.tank_pressure);
            f.a = -f.tank;
            f.b = f.supply;
        }
        return f;
    }

    // Inverse of OrificeFlow scaled by the spool opening: pressure drop passing flow q through opening u.
    double PressureDrop(double q, double u) const {
        double w = q / (u * FlowGain());
        if (w == 0)
            return 0;
        if (!std::isfinite(w))
            return w;
        double w4 = w * w * w * w;
        double pt = m_p.transition_pressure;
        double dp2 = 0.5 * w4 * (1 + std::sqrt(1 + 4 * pt * pt / w4));
        return std::copysign(std::sqrt(dp2), w);
    }

    void QuasiStaticPressures(const ActuatorCommand& command, double v, double& pa, double& pb) const {
        double u = Spool(command);
        if (u == 0 && v == 0) {
            pa = pb = 0.5 * (m_p.supply_pressure + m_p.tank_pressure);
            return;
        }
        double qa = m_p.cap_area * v;       // required inflow to A
        double qb = -m_p.annulus_area * v;  // required inflow to B
        if (u > 0) {
            pa = m_p.supply_pressure - PressureDrop(qa, u);
            pb = m_p.tank_pressure + PressureDrop(-qb, u);
        } else if (u < 0) {
            pa = m_p.tank_pressure + PressureDrop(-qa, -u);
            pb = m_p.supply_pressure - PressureDrop(qb, -u);
        } else {
            pa = v > 0 ? m_p.tank_pressure : m_p.relief_pressure;
            pb = v > 0 ? m_p.relief_pressure : m_p.tank_pressure;
        }
        pa = BoundPressure(pa);
        pb = BoundPressure(pb);
    }

    double EffortFromPressures(double pa, double pb, double v) const {
        return pa * m_p.cap_area - pb * m_p.annulus_area - m_p.coulomb_friction * std::tanh(v / m_p.friction_smoothing_speed) - m_p.viscous_friction * v;
    }

    ActuatorTelemetry Telemetry(double pa, double pb, const ActuatorCommand& command, const ActuatorState& state) const {
        ActuatorTelemetry t;
        Flows f = ValveFlows(Spool(command), pa, pb);
        t.applied_effort = EffortFromPressures(pa, pb, state.velocity);
        t.mechanical_power = t.applied_effort * state.velocity;
        t.electrical_power = m_p.supply_pressure * f.supply - m_p.tank_pressure * f.tank;
        if (t.mechanical_power > 0 && t.electrical_power > 0)
            t.efficiency = std::clamp(t.mechanical_power / t.electrical_power, 0.0, 1.0);
        else if (t.mechanical_power < 0 && t.electrical_power < 0)
            t.efficiency = std::clamp(t.electrical_power / t.mechanical_power, 0.0, 1.0);
        t.command_saturated = command.enabled && std::abs(command.value) > 1;
        t.effort_saturated = pa >= m_p.relief_pressure || pb >= m_p.relief_pressure;
        return t;
    }

    Parameters m_p;
};

}  // namespace actuator_models
}  // namespace chrono

#endif
