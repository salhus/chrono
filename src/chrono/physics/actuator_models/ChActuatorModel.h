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

#ifndef CH_ACTUATOR_MODEL_H
#define CH_ACTUATOR_MODEL_H

#include <cstddef>
#include <vector>

namespace chrono {
namespace actuator_models {

/// Actuator command. The default enabled value makes a newly constructed command simulation-ready.
/// The model interprets value as normalized input unless otherwise documented; disabled-command behavior is model-specific.
struct ActuatorCommand {
    double value = 0;
    bool enabled = true;
};

/// Kinematic state supplied by the containing simulation.
struct ActuatorState {
    double position = 0;
    double velocity = 0;
    double time = 0;
};

/// Optional actuator output and power telemetry.
struct ActuatorTelemetry {
    double applied_effort = 0;
    double mechanical_power = 0;
    double electrical_power = 0;
    double efficiency = 0;
    bool command_saturated = false;
    bool current_saturated = false;
    bool effort_saturated = false;
    bool voltage_saturated = false;
};

/// Engine-independent actuator model interface.
/// Positive effort extends a linear actuator or increases a rotational coordinate. Resistive effort has the opposite
/// sign to velocity.
class ChActuatorModel {
  public:
    virtual ~ChActuatorModel() = default;

    /// Compute effort for models without internal state. This function must not mutate model state.
    virtual double ComputeEffort(const ActuatorCommand& command, const ActuatorState& state) const = 0;

    /// Number of internal states (zero for stateless models).
    virtual std::size_t GetNumStates() const { return 0; }

    /// Set default initial conditions for internal states.
    virtual void SetInitialConditions(std::vector<double>& y0) const { y0.assign(GetNumStates(), 0); }

    /// Evaluate internal state derivatives.
    virtual void CalculateRHS(double, const std::vector<double>&, const ActuatorCommand&, const ActuatorState&, std::vector<double>& rhs) const { rhs.assign(GetNumStates(), 0); }

    /// Optionally provide the internal-state Jacobian, row-major.
    virtual bool CalculateJacobian(double, const std::vector<double>&, const ActuatorCommand&, const ActuatorState&, std::vector<std::vector<double>>&) const { return false; }

    /// Declare whether internal dynamics are stiff.
    virtual bool IsStiff() const { return false; }

    /// Compute effort from internal states; stateless models use ComputeEffort.
    virtual double EffortFromStates(const std::vector<double>&, const ActuatorCommand& command, const ActuatorState& state) const { return ComputeEffort(command, state); }

    /// Compute optional telemetry for the current model output.
    virtual ActuatorTelemetry ComputeTelemetry(const ActuatorCommand& command, const ActuatorState& state) const {
        ActuatorTelemetry telemetry;
        telemetry.applied_effort = ComputeEffort(command, state);
        telemetry.mechanical_power = telemetry.applied_effort * state.velocity;
        return telemetry;
    }
};

}  // namespace actuator_models
}  // namespace chrono

#endif
