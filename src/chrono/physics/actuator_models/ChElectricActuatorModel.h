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

#ifndef CH_ELECTRIC_ACTUATOR_MODEL_H
#define CH_ELECTRIC_ACTUATOR_MODEL_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "ChActuatorModel.h"

namespace chrono {
namespace actuator_models {

/// DC-equivalent electric actuator with optional winding-current dynamics and a passive, bidirectional gearbox.
/// Positive terminal power is drawn from the bus. Disabling selects zero voltage (shorted winding), not an open circuit.
/// All evaluations are pure; winding current is owned and integrated by the caller when inductance is positive.
class ChElectricActuatorModel : public ChActuatorModel {
  public:
    /// SI parameters: output-side velocity/effort and motor-side electrical constants.
    /// An ideal SI motor requires Kt == Ke: electromagnetic power Kt*i*omega must equal Ke*omega*i.
    /// The explicit opt-out is only for externally converted/non-SI data; it forfeits the SI passivity guarantee.
    struct Parameters {
        double resistance = 1;               ///< Winding resistance [ohm], finite and positive.
        double inductance = 0.01;            ///< Winding inductance [H], finite and nonnegative; zero is quasi-static.
        double torque_constant = 1;          ///< Motor torque/current [N*m/A] (or linear force/current [N/A]).
        double back_emf_constant = 1;        ///< Back EMF/speed [V*s/rad] (or [V*s/m]), finite and positive.
        double bus_voltage = 24;             ///< Maximum terminal voltage magnitude [V], finite and positive.
        double gear_ratio = 1;               ///< Motor speed/output speed, finite and positive.
        double motoring_efficiency = 0.9;    ///< Gearbox motoring efficiency in (0, 1].
        double generating_efficiency = 0.8;  ///< Gearbox generating efficiency in (0, 1].
        double max_current = 10;             ///< Symmetric winding-current limit [A], finite and positive.
        /// Electromagnetic output effort ceiling [N*m or N], positive; passive friction is additional.
        /// Enforced conservatively by a fixed current cap max_effort*generating_efficiency/(gear_ratio*Kt).
        double max_effort = std::numeric_limits<double>::infinity();
        double coulomb_friction = 0;                ///< Output Coulomb friction magnitude [N*m or N], finite and nonnegative.
        double viscous_friction = 0;                ///< Output viscous coefficient [N*m*s/rad or N*s/m], nonnegative.
        double friction_smoothing_speed = 1e-3;     ///< Positive, finite Coulomb smoothing speed [rad/s or m/s].
        bool allow_inconsistent_constants = false;  ///< Explicitly waive the relative 1e-6 Kt/Ke consistency check.
        double quadrant_smoothing = 1e-3;           ///< Positive, finite electromagnetic-power crossover width [W].
    };

    /// Construct with the default SI parameters.
    ChElectricActuatorModel() : ChElectricActuatorModel(Parameters{}) {}

    /// Copy and validate parameters; throw std::invalid_argument for invalid values or inconsistent SI constants.
    explicit ChElectricActuatorModel(const Parameters& parameters) : m_parameters(parameters) { ValidateParameters(); }

    /// Evaluate the quasi-static current solution, even if the configured model also exposes a winding state.
    double ComputeEffort(const ActuatorCommand& command, const ActuatorState& state) const override {
        double current = SteadyStateCurrent(command, state.velocity);
        return EffortFromCurrent(current, state.velocity);
    }

    /// One winding-current state [A] for positive inductance; otherwise zero states.
    std::size_t GetNumStates() const override { return m_parameters.inductance > 0 ? 1 : 0; }

    /// Initialize winding current to zero, or return an empty vector for the quasi-static model.
    void SetInitialConditions(std::vector<double>& y0) const override { y0.assign(GetNumStates(), 0); }

    /// Current limiting blocks outward derivatives at either bound; inward derivatives remain unconstrained.
    /// This represents a dissipative ideal current limiter, not an open circuit or a regenerative bus controller.
    /// An integrator must resolve limit crossings; effort and telemetry project numerical overshoots onto the bounds.
    void CalculateRHS(double, const std::vector<double>& y, const ActuatorCommand& command, const ActuatorState& state, std::vector<double>& rhs) const override {
        rhs.assign(GetNumStates(), 0);
        if (GetNumStates() == 0)
            return;

        double current = y.at(0);
        double voltage = AppliedVoltage(command, nullptr);
        double motor_speed = m_parameters.gear_ratio * state.velocity;
        double current_rate = (voltage - m_parameters.resistance * current - m_parameters.back_emf_constant * motor_speed) / m_parameters.inductance;
        double limit = CurrentLimit();
        if ((current >= limit && current_rate > 0) || (current <= -limit && current_rate < 0))
            current_rate = 0;
        rhs[0] = current_rate;
    }

    /// Exact state derivative on each smooth branch; the active outward clamp has a zero Jacobian row.
    bool CalculateJacobian(double time,
                           const std::vector<double>& y,
                           const ActuatorCommand& command,
                           const ActuatorState& state,
                           std::vector<std::vector<double>>& jacobian) const override {
        if (GetNumStates() == 0)
            return false;
        std::vector<double> rhs;
        CalculateRHS(time, y, command, state, rhs);
        bool clamped = std::abs(y.at(0)) >= CurrentLimit() && rhs[0] == 0;
        jacobian.assign(1, std::vector<double>(1, clamped ? 0 : -m_parameters.resistance / m_parameters.inductance));
        return true;
    }

    /// Winding dynamics require a stiff solver whenever inductance is positive.
    bool IsStiff() const override { return GetNumStates() > 0; }

    /// Evaluate effort using the bounded integrated winding current; a zero-state model uses the algebraic solution.
    double EffortFromStates(const std::vector<double>& y, const ActuatorCommand& command, const ActuatorState& state) const override {
        if (GetNumStates() == 0)
            return ComputeEffort(command, state);
        double current = std::clamp(y.at(0), -CurrentLimit(), CurrentLimit());
        return EffortFromCurrent(current, state.velocity);
    }

    /// Return quasi-static terminal power V*i, copper loss R*i^2, effort, efficiency, and active-limit flags.
    ActuatorTelemetry ComputeTelemetry(const ActuatorCommand& command, const ActuatorState& state) const override {
        return TelemetryFromCurrent(SteadyStateCurrent(command, state.velocity), command, state);
    }

    /// Use the integrated (bounded) winding current, not its quasi-static solution.
    ActuatorTelemetry ComputeTelemetryFromStates(const std::vector<double>& y, const ActuatorCommand& command, const ActuatorState& state) const override {
        return GetNumStates() == 0 ? ComputeTelemetry(command, state) : TelemetryFromCurrent(y.at(0), command, state);
    }

    /// Output stall effort at full voltage, limited by both winding resistance and the current bound.
    /// At zero power the smooth gearbox is lossless; far from crossover its motoring efficiency applies.
    double GetStallEffort() const { return ComputeEffort({1, true}, {}); }

    /// Ideal, frictionless output no-load speed at full bus voltage [rad/s or m/s].
    double GetNoLoadSpeed() const { return m_parameters.bus_voltage / (m_parameters.back_emf_constant * m_parameters.gear_ratio); }

  private:
    ActuatorTelemetry TelemetryFromCurrent(double requested_current, const ActuatorCommand& command, const ActuatorState& state) const {
        ActuatorTelemetry telemetry;
        bool voltage_saturated = false;
        double voltage = AppliedVoltage(command, &voltage_saturated);
        double limit = CurrentLimit();
        double current = std::clamp(requested_current, -limit, limit);
        telemetry.applied_effort = EffortFromCurrent(current, state.velocity);
        telemetry.mechanical_power = telemetry.applied_effort * state.velocity;
        telemetry.electrical_power = voltage * current;
        telemetry.current = current;
        telemetry.copper_loss = current * current * m_parameters.resistance;
        if (telemetry.mechanical_power > 0 && telemetry.electrical_power > 0)
            telemetry.efficiency = std::clamp(telemetry.mechanical_power / telemetry.electrical_power, 0.0, 1.0);
        else if (telemetry.mechanical_power < 0 && telemetry.electrical_power < 0)
            telemetry.efficiency = std::clamp(telemetry.electrical_power / telemetry.mechanical_power, 0.0, 1.0);
        telemetry.command_saturated = command.enabled && std::abs(command.value) > 1;
        telemetry.current_saturated = std::abs(requested_current) >= limit;
        telemetry.effort_saturated = telemetry.current_saturated && EffortCurrentLimit() < m_parameters.max_current;
        telemetry.voltage_saturated = voltage_saturated;
        return telemetry;
    }

    void ValidateParameters() const {
        if (!std::isfinite(m_parameters.resistance) || !(m_parameters.resistance > 0) || !std::isfinite(m_parameters.inductance) || !(m_parameters.inductance >= 0) ||
            !std::isfinite(m_parameters.torque_constant) || !(m_parameters.torque_constant > 0) || !std::isfinite(m_parameters.back_emf_constant) ||
            !(m_parameters.back_emf_constant > 0) || !std::isfinite(m_parameters.bus_voltage) || !(m_parameters.bus_voltage > 0) || !std::isfinite(m_parameters.gear_ratio) ||
            !(m_parameters.gear_ratio > 0) || !std::isfinite(m_parameters.max_current) || !(m_parameters.max_current > 0) || !(m_parameters.max_effort > 0) ||
            !std::isfinite(m_parameters.friction_smoothing_speed) || !(m_parameters.friction_smoothing_speed > 0) || !std::isfinite(m_parameters.coulomb_friction) ||
            m_parameters.coulomb_friction < 0 || !std::isfinite(m_parameters.viscous_friction) || m_parameters.viscous_friction < 0 ||
            !(m_parameters.motoring_efficiency > 0 && m_parameters.motoring_efficiency <= 1) ||
            !(m_parameters.generating_efficiency > 0 && m_parameters.generating_efficiency <= 1))
            throw std::invalid_argument("Invalid electric actuator model parameters");
        if (!std::isfinite(m_parameters.quadrant_smoothing) || m_parameters.quadrant_smoothing <= 0)
            throw std::invalid_argument("Quadrant smoothing power must be positive and finite");
        if (!m_parameters.allow_inconsistent_constants &&
            std::abs(m_parameters.torque_constant - m_parameters.back_emf_constant) > 1e-6 * std::max(m_parameters.torque_constant, m_parameters.back_emf_constant))
            throw std::invalid_argument("SI torque and back-EMF constants must match (or explicitly opt out)");
    }

    double AppliedVoltage(const ActuatorCommand& command, bool* saturated) const {
        double requested_voltage = command.enabled ? command.value * m_parameters.bus_voltage : 0;
        double voltage = std::clamp(requested_voltage, -m_parameters.bus_voltage, m_parameters.bus_voltage);
        if (saturated)
            *saturated = voltage != requested_voltage;
        return voltage;
    }

    double EffortCurrentLimit() const { return m_parameters.max_effort * m_parameters.generating_efficiency / (m_parameters.gear_ratio * m_parameters.torque_constant); }

    double CurrentLimit() const { return std::min(m_parameters.max_current, EffortCurrentLimit()); }

    double SteadyStateCurrent(const ActuatorCommand& command, double output_speed) const {
        double voltage = AppliedVoltage(command, nullptr);
        double motor_speed = m_parameters.gear_ratio * output_speed;
        double requested_current = (voltage - m_parameters.back_emf_constant * motor_speed) / m_parameters.resistance;
        return std::clamp(requested_current, -CurrentLimit(), CurrentLimit());
    }

    double EffortFromCurrent(double current, double output_speed) const {
        double motor_effort = m_parameters.torque_constant * current;
        double motor_speed = m_parameters.gear_ratio * output_speed;
        double power = motor_effort * motor_speed;
        // Blend through unity at zero power: a direct blend between the two efficiencies can create energy.
        double q = std::min(std::abs(power) / m_parameters.quadrant_smoothing, 1.0);
        double blend = q * q * (3 - 2 * q);
        double factor = power >= 0 ? 1 + (m_parameters.motoring_efficiency - 1) * blend : 1 + (1 / m_parameters.generating_efficiency - 1) * blend;
        double requested_effort = motor_effort * m_parameters.gear_ratio * factor;
        requested_effort -= m_parameters.coulomb_friction * std::tanh(output_speed / m_parameters.friction_smoothing_speed);
        requested_effort -= m_parameters.viscous_friction * output_speed;
        return requested_effort;
    }

    Parameters m_parameters;
};

}  // namespace actuator_models
}  // namespace chrono

#endif
