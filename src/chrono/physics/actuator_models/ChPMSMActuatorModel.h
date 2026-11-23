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

#ifndef CH_PMSM_ACTUATOR_MODEL_H
#define CH_PMSM_ACTUATOR_MODEL_H

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "ChActuatorModel.h"
#include "ChActuatorSmoothing.h"

namespace chrono {
namespace actuator_models {

/// Permanent-magnet synchronous machine with field-oriented current control (torque mode) and a passive gearbox.
///
/// Rotor-frame (dq, amplitude-invariant Park) model:
///   Ld di_d/dt = v_d - R i_d + w_e Lq i_q
///   Lq di_q/dt = v_q - R i_q - w_e (Ld i_d + psi),     w_e = p * N * output_speed
///   T_motor    = 1.5 p (psi i_q + (Ld - Lq) i_d i_q)
/// Two PI current loops track i_d* = 0 and i_q* = command * max_current, with gains Kp = w_c L and Ki = w_c R
/// (pole-zero cancellation, first-order closed loop of bandwidth w_c) and speed-voltage decoupling feedforward.
/// The requested voltage vector is scaled onto the linear space-vector limit bus_voltage / sqrt(3); back-calculation
/// (gain 1/Kp) prevents integrator windup. Field weakening and MTPA are not modeled (i_d* = 0).
/// Disabling the command selects zero phase voltage (active short circuit), and integrators relax to zero at w_c.
///
/// States: [i_d, i_q, z_d, z_q], with z the integrals of the current errors [A*s]. All evaluations are pure.
class ChPMSMActuatorModel : public ChActuatorModel {
  public:
    enum StateIndex { I_D = 0, I_Q = 1, Z_D = 2, Z_Q = 3 };

    /// SI parameters. Electrical quantities are per phase and motor side; effort and speed are output side.
    struct Parameters {
        int pole_pairs = 4;                      ///< Number of pole pairs p, positive.
        double resistance = 0.5;                 ///< Phase resistance [ohm], positive.
        double inductance_d = 1e-3;              ///< d-axis inductance [H], positive.
        double inductance_q = 1e-3;              ///< q-axis inductance [H], positive.
        double flux_linkage = 0.02;              ///< Permanent-magnet flux linkage amplitude [Wb], positive.
        double bus_voltage = 48;                 ///< DC-link voltage [V], positive; phase amplitude limit is V/sqrt(3).
        double max_current = 10;                 ///< Phase current amplitude limit [A], positive (sets i_q* range).
        double current_bandwidth = 2000;         ///< Current-loop closed-loop bandwidth w_c [rad/s], positive.
        double gear_ratio = 1;                   ///< Motor speed/output speed, positive.
        double motoring_efficiency = 1;          ///< Gearbox motoring efficiency in (0, 1].
        double generating_efficiency = 1;        ///< Gearbox generating efficiency in (0, 1].
        double loss_smoothing_speed = 1e-2;      ///< Positive transmission crossover speed [rad/s or m/s].
        double coulomb_friction = 0;             ///< Output Coulomb friction magnitude [N*m or N], nonnegative.
        double viscous_friction = 0;             ///< Output viscous coefficient [N*m*s/rad or N*s/m], nonnegative.
        double friction_smoothing_speed = 1e-3;  ///< Positive Coulomb smoothing speed [rad/s or m/s].
        double rotor_inertia = 0;                ///< Motor-side rotor inertia [kg*m^2], nonnegative.
    };

    /// Construct with the default SI parameters.
    ChPMSMActuatorModel() : ChPMSMActuatorModel(Parameters{}) {}

    /// Copy and validate parameters; throw std::invalid_argument for invalid values.
    explicit ChPMSMActuatorModel(const Parameters& parameters) : m_p(parameters) { ValidateParameters(); }

    /// Quasi-static effort: ideal current tracking with i_d = 0, i_q* reduced to the largest magnitude whose
    /// steady terminal voltage fits the inverter limit (minimum-voltage current beyond that). Disabled: the steady
    /// short-circuit braking effort. Exact steady state of the closed loop whenever the voltage is not saturated.
    double ComputeEffort(const ActuatorCommand& command, const ActuatorState& state) const override {
        double id, iq;
        SteadyCurrents(command, state.velocity, id, iq, nullptr);
        return EffortFromTorque(MotorTorque(id, iq), state.velocity);
    }

    std::size_t GetNumStates() const override { return 4; }

    void SetInitialConditions(std::vector<double>& y0) const override { y0.assign(4, 0); }

    void CalculateRHS(double, const std::vector<double>& y, const ActuatorCommand& command, const ActuatorState& state, std::vector<double>& rhs) const override {
        CheckStates(y);
        double we = ElectricalSpeed(state.velocity);
        Drive d = Controller(y, command, we);
        rhs.assign(4, 0);
        rhs[I_D] = (d.vd - m_p.resistance * y[I_D] + we * m_p.inductance_q * y[I_Q]) / m_p.inductance_d;
        rhs[I_Q] = (d.vq - m_p.resistance * y[I_Q] - we * (m_p.inductance_d * y[I_D] + m_p.flux_linkage)) / m_p.inductance_q;
        if (command.enabled) {
            rhs[Z_D] = d.ed + (d.vd - d.ud) / (m_p.current_bandwidth * m_p.inductance_d);
            rhs[Z_Q] = d.eq + (d.vq - d.uq) / (m_p.current_bandwidth * m_p.inductance_q);
        } else {
            rhs[Z_D] = -m_p.current_bandwidth * y[Z_D];
            rhs[Z_Q] = -m_p.current_bandwidth * y[Z_Q];
        }
    }

    /// Current dynamics are stiff (time constants L/R and 1/w_c).
    bool IsStiff() const override { return true; }

    /// Rotor inertia reflected through the gear: gear_ratio^2 * rotor_inertia.
    double GetReflectedInertia() const override { return m_p.gear_ratio * m_p.gear_ratio * m_p.rotor_inertia; }

    double EffortFromStates(const std::vector<double>& y, const ActuatorCommand&, const ActuatorState& state) const override {
        CheckStates(y);
        return EffortFromTorque(MotorTorque(y[I_D], y[I_Q]), state.velocity);
    }

    /// Quasi-static telemetry; electrical_power is the inverter output 1.5 (v_d i_d + v_q i_q) and current the
    /// phase-current amplitude.
    ActuatorTelemetry ComputeTelemetry(const ActuatorCommand& command, const ActuatorState& state) const override {
        double id, iq;
        bool voltage_limited = false;
        SteadyCurrents(command, state.velocity, id, iq, &voltage_limited);
        double we = ElectricalSpeed(state.velocity);
        double vd = m_p.resistance * id - we * m_p.inductance_q * iq;
        double vq = m_p.resistance * iq + we * (m_p.inductance_d * id + m_p.flux_linkage);
        if (!command.enabled)
            vd = vq = 0;
        return Telemetry(id, iq, vd, vq, voltage_limited, command, state);
    }

    /// Telemetry from integrated currents and the controller's applied voltage.
    ActuatorTelemetry ComputeTelemetryFromStates(const std::vector<double>& y, const ActuatorCommand& command, const ActuatorState& state) const override {
        CheckStates(y);
        Drive d = Controller(y, command, ElectricalSpeed(state.velocity));
        return Telemetry(y[I_D], y[I_Q], d.vd, d.vq, d.saturated, command, state);
    }

    /// Motor torque per q-axis current amplitude, 1.5 p psi [N*m/A].
    double GetTorqueConstant() const { return 1.5 * m_p.pole_pairs * m_p.flux_linkage; }

    /// Phase-voltage amplitude limit, bus_voltage / sqrt(3) [V].
    double GetVoltageLimit() const { return m_p.bus_voltage / std::sqrt(3.0); }

    /// Ideal output no-load speed (back-EMF amplitude equal to the voltage limit) [rad/s or m/s].
    double GetNoLoadSpeed() const { return GetVoltageLimit() / (m_p.pole_pairs * m_p.flux_linkage * m_p.gear_ratio); }

    /// Ideal output peak effort at max_current, before transmission loss and friction [N*m or N].
    double GetPeakEffort() const { return m_p.gear_ratio * GetTorqueConstant() * m_p.max_current; }

    /// Magnetic energy stored in the windings, 0.75 (Ld i_d^2 + Lq i_q^2) [J].
    double GetMagneticEnergy(const std::vector<double>& y) const {
        CheckStates(y);
        return 0.75 * (m_p.inductance_d * y[I_D] * y[I_D] + m_p.inductance_q * y[I_Q] * y[I_Q]);
    }

    const Parameters& GetParameters() const { return m_p; }

  private:
    struct Drive {
        double vd = 0, vq = 0;  // applied voltage
        double ud = 0, uq = 0;  // requested (unsaturated) voltage
        double ed = 0, eq = 0;  // current errors
        bool saturated = false;
    };

    void ValidateParameters() const {
        auto positive = [](double x) { return std::isfinite(x) && x > 0; };
        auto nonnegative = [](double x) { return std::isfinite(x) && x >= 0; };
        if (m_p.pole_pairs < 1 || !positive(m_p.resistance) || !positive(m_p.inductance_d) || !positive(m_p.inductance_q) ||
            !positive(m_p.flux_linkage) || !positive(m_p.bus_voltage) || !positive(m_p.max_current) || !positive(m_p.current_bandwidth) ||
            !positive(m_p.gear_ratio) || !positive(m_p.loss_smoothing_speed) || !positive(m_p.friction_smoothing_speed) ||
            !nonnegative(m_p.coulomb_friction) || !nonnegative(m_p.viscous_friction) || !nonnegative(m_p.rotor_inertia) ||
            !(m_p.motoring_efficiency > 0 && m_p.motoring_efficiency <= 1) || !(m_p.generating_efficiency > 0 && m_p.generating_efficiency <= 1))
            throw std::invalid_argument("Invalid PMSM actuator model parameters");
    }

    void CheckStates(const std::vector<double>& y) const {
        if (y.size() != 4)
            throw std::invalid_argument("PMSM actuator model requires 4 states");
    }

    double ElectricalSpeed(double output_speed) const { return m_p.pole_pairs * m_p.gear_ratio * output_speed; }

    double CurrentReference(const ActuatorCommand& command) const {
        return command.enabled ? std::clamp(command.value, -1.0, 1.0) * m_p.max_current : 0;
    }

    Drive Controller(const std::vector<double>& y, const ActuatorCommand& command, double we) const {
        Drive d;
        if (!command.enabled)
            return d;
        double ki = m_p.current_bandwidth * m_p.resistance;
        d.ed = -y[I_D];
        d.eq = CurrentReference(command) - y[I_Q];
        d.ud = m_p.current_bandwidth * m_p.inductance_d * d.ed + ki * y[Z_D] - we * m_p.inductance_q * y[I_Q];
        d.uq = m_p.current_bandwidth * m_p.inductance_q * d.eq + ki * y[Z_Q] + we * (m_p.inductance_d * y[I_D] + m_p.flux_linkage);
        double magnitude = std::hypot(d.ud, d.uq);
        double limit = GetVoltageLimit();
        double scale = magnitude > limit ? limit / magnitude : 1;
        d.vd = scale * d.ud;
        d.vq = scale * d.uq;
        d.saturated = scale < 1;
        return d;
    }

    void SteadyCurrents(const ActuatorCommand& command, double output_speed, double& id, double& iq, bool* voltage_limited) const {
        double we = ElectricalSpeed(output_speed);
        double R = m_p.resistance, Ld = m_p.inductance_d, Lq = m_p.inductance_q, psi = m_p.flux_linkage;
        if (voltage_limited)
            *voltage_limited = false;
        if (!command.enabled) {
            double det = R * R + we * we * Ld * Lq;
            id = -we * we * Lq * psi / det;
            iq = -R * we * psi / det;
            return;
        }
        // With i_d = 0: |v|^2 = (we Lq iq)^2 + (R iq + we psi)^2 <= Vmax^2, a quadratic in iq.
        double vmax = GetVoltageLimit();
        double a = R * R + we * we * Lq * Lq;
        double b = 2 * R * we * psi;
        double c = we * we * psi * psi - vmax * vmax;
        double discriminant = b * b - 4 * a * c;
        double reference = CurrentReference(command);
        id = 0;
        if (discriminant < 0) {
            iq = -b / (2 * a);
        } else {
            double root = std::sqrt(discriminant);
            iq = std::clamp(reference, (-b - root) / (2 * a), (-b + root) / (2 * a));
        }
        if (voltage_limited)
            *voltage_limited = iq != reference;
    }

    double MotorTorque(double id, double iq) const {
        return 1.5 * m_p.pole_pairs * (m_p.flux_linkage * iq + (m_p.inductance_d - m_p.inductance_q) * id * iq);
    }

    double GearLoss(double ideal_effort, double output_speed) const {
        return ActuatorTransmissionLoss(ideal_effort, output_speed, m_p.loss_smoothing_speed, 1 - m_p.motoring_efficiency, 1 / m_p.generating_efficiency - 1);
    }

    double EffortFromTorque(double motor_torque, double output_speed) const {
        double ideal = m_p.gear_ratio * motor_torque;
        return ideal - GearLoss(ideal, output_speed) - m_p.coulomb_friction * std::tanh(output_speed / m_p.friction_smoothing_speed) -
               m_p.viscous_friction * output_speed;
    }

    ActuatorTelemetry Telemetry(double id, double iq, double vd, double vq, bool voltage_saturated, const ActuatorCommand& command, const ActuatorState& state) const {
        ActuatorTelemetry t;
        double ideal = m_p.gear_ratio * MotorTorque(id, iq);
        t.applied_effort = EffortFromTorque(MotorTorque(id, iq), state.velocity);
        t.mechanical_power = t.applied_effort * state.velocity;
        t.electrical_power = 1.5 * (vd * id + vq * iq);
        t.current = std::hypot(id, iq);
        t.copper_loss = 1.5 * m_p.resistance * (id * id + iq * iq);
        t.gear_loss = GearLoss(ideal, state.velocity) * state.velocity;
        if (t.mechanical_power > 0 && t.electrical_power > 0)
            t.efficiency = std::clamp(t.mechanical_power / t.electrical_power, 0.0, 1.0);
        else if (t.mechanical_power < 0 && t.electrical_power < 0)
            t.efficiency = std::clamp(t.electrical_power / t.mechanical_power, 0.0, 1.0);
        t.command_saturated = command.enabled && std::abs(command.value) > 1;
        t.current_saturated = command.enabled && std::abs(command.value) >= 1;
        t.voltage_saturated = voltage_saturated;
        return t;
    }

    Parameters m_p;
};

}  // namespace actuator_models
}  // namespace chrono

#endif
