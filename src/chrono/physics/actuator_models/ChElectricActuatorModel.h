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

/// DC electric actuator model with winding-current dynamics or a quasi-static current solution.
class ChElectricActuatorModel : public ChActuatorModel {
  public:
    /// State velocity and effort are output-side quantities; torque and back-EMF constants are motor-side values.
    struct Parameters {
        double resistance = 1;
        double inductance = 0.01;
        double torque_constant = 1;
        double back_emf_constant = 1;
        double bus_voltage = 24;
        double gear_ratio = 1;
        double motoring_efficiency = 0.9;
        double generating_efficiency = 0.8;
        double max_current = 10;
        double max_effort = std::numeric_limits<double>::infinity();
        double coulomb_friction = 0;
        double viscous_friction = 0;
        double friction_smoothing_speed = 1e-3;
    };

    ChElectricActuatorModel() : ChElectricActuatorModel(Parameters{}) {}

    explicit ChElectricActuatorModel(const Parameters& parameters) : m_parameters(parameters) { ValidateParameters(); }

    double ComputeEffort(const ActuatorCommand& command, const ActuatorState& state) const override {
        double current = SteadyStateCurrent(command, state.velocity);
        return EffortFromCurrent(current, state.velocity, nullptr, nullptr);
    }

    std::size_t GetNumStates() const override { return m_parameters.inductance > 0 ? 1 : 0; }

    void SetInitialConditions(std::vector<double>& y0) const override { y0.assign(GetNumStates(), 0); }

    void CalculateRHS(double, const std::vector<double>& y, const ActuatorCommand& command, const ActuatorState& state, std::vector<double>& rhs) const override {
        rhs.assign(GetNumStates(), 0);
        if (GetNumStates() == 0)
            return;

        double current = y.at(0);
        double voltage = AppliedVoltage(command, nullptr);
        double motor_speed = m_parameters.gear_ratio * state.velocity;
        double current_rate = (voltage - m_parameters.resistance * current - m_parameters.back_emf_constant * motor_speed) / m_parameters.inductance;
        if ((current >= m_parameters.max_current && current_rate > 0) || (current <= -m_parameters.max_current && current_rate < 0))
            current_rate = 0;
        rhs[0] = current_rate;
    }

    bool CalculateJacobian(double, const std::vector<double>&, const ActuatorCommand&, const ActuatorState&, std::vector<std::vector<double>>& jacobian) const override {
        if (GetNumStates() == 0)
            return false;
        jacobian.assign(1, std::vector<double>(1, -m_parameters.resistance / m_parameters.inductance));
        return true;
    }

    bool IsStiff() const override { return GetNumStates() > 0; }

    double EffortFromStates(const std::vector<double>& y, const ActuatorCommand& command, const ActuatorState& state) const override {
        if (GetNumStates() == 0)
            return ComputeEffort(command, state);
        double current = std::clamp(y.at(0), -m_parameters.max_current, m_parameters.max_current);
        return EffortFromCurrent(current, state.velocity, nullptr, nullptr);
    }

    ActuatorTelemetry ComputeTelemetry(const ActuatorCommand& command, const ActuatorState& state) const override {
        ActuatorTelemetry telemetry;
        bool voltage_saturated = false;
        bool current_saturated = false;
        bool effort_saturated = false;
        bool command_saturated = !command.enabled && command.value != 0;
        double voltage = AppliedVoltage(command, &voltage_saturated);
        double current = SteadyStateCurrent(command, state.velocity, &current_saturated, &command_saturated);
        telemetry.applied_effort = EffortFromCurrent(current, state.velocity, &effort_saturated, &current_saturated);
        telemetry.mechanical_power = telemetry.applied_effort * state.velocity;
        telemetry.electrical_power = voltage * current;
        if (std::abs(telemetry.electrical_power) > std::numeric_limits<double>::epsilon()) {
            telemetry.efficiency = std::min(1.0, std::abs(telemetry.mechanical_power / telemetry.electrical_power));
        }
        telemetry.command_saturated = command_saturated;
        telemetry.current_saturated = current_saturated;
        telemetry.effort_saturated = effort_saturated;
        telemetry.voltage_saturated = voltage_saturated;
        return telemetry;
    }

    /// Maximum output effort at zero speed and current limit.
    double GetStallEffort() const {
        return std::min(m_parameters.max_effort, m_parameters.gear_ratio * m_parameters.motoring_efficiency * m_parameters.torque_constant * m_parameters.max_current);
    }

    /// Ideal output no-load speed at full bus voltage.
    double GetNoLoadSpeed() const { return m_parameters.bus_voltage / (m_parameters.back_emf_constant * m_parameters.gear_ratio); }

  private:
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
    }

    double AppliedVoltage(const ActuatorCommand& command, bool* saturated) const {
        double requested_voltage = command.enabled ? command.value * m_parameters.bus_voltage : 0;
        double voltage = std::clamp(requested_voltage, -m_parameters.bus_voltage, m_parameters.bus_voltage);
        if (saturated)
            *saturated = voltage != requested_voltage;
        return voltage;
    }

    double SteadyStateCurrent(const ActuatorCommand& command, double output_speed, bool* saturated = nullptr, bool* command_saturated = nullptr) const {
        if (command_saturated)
            *command_saturated = !command.enabled && command.value != 0;
        bool voltage_saturated = false;
        double voltage = AppliedVoltage(command, &voltage_saturated);
        if (command_saturated)
            *command_saturated = *command_saturated || voltage_saturated;
        double motor_speed = m_parameters.gear_ratio * output_speed;
        double requested_current = (voltage - m_parameters.back_emf_constant * motor_speed) / m_parameters.resistance;
        double current = std::clamp(requested_current, -m_parameters.max_current, m_parameters.max_current);
        if (saturated)
            *saturated = current != requested_current;
        return current;
    }

    double EffortFromCurrent(double current, double output_speed, bool* effort_saturated, bool* current_saturated) const {
        double requested_current = current;
        current = std::clamp(current, -m_parameters.max_current, m_parameters.max_current);
        if (current_saturated)
            *current_saturated = *current_saturated || current != requested_current;

        double motor_effort = m_parameters.torque_constant * current;
        double motor_speed = m_parameters.gear_ratio * output_speed;
        double efficiency = motor_effort * motor_speed >= 0 ? m_parameters.motoring_efficiency : m_parameters.generating_efficiency;
        double requested_effort = motor_effort * m_parameters.gear_ratio * (motor_effort * motor_speed >= 0 ? efficiency : 1 / efficiency);
        requested_effort -= m_parameters.coulomb_friction * std::tanh(output_speed / m_parameters.friction_smoothing_speed);
        requested_effort -= m_parameters.viscous_friction * output_speed;
        double effort = std::clamp(requested_effort, -m_parameters.max_effort, m_parameters.max_effort);
        if (effort_saturated)
            *effort_saturated = effort != requested_effort;
        return effort;
    }

    Parameters m_parameters;
};

}  // namespace actuator_models
}  // namespace chrono

#endif
