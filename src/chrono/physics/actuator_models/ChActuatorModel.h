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

/// Normalized actuator input, clamped to [-1, 1] by the supplied models.
/// Commands are enabled by default. Disabling removes active input, not passive friction or back-driving.
struct ActuatorCommand {
    double value = 0;     ///< Fraction of bus voltage (electric) or available effort (tabulated).
    bool enabled = true;  ///< Apply the input; false selects zero voltage/commanded effort.
};

/// Output-side kinematics in SI units, supplied by the containing simulation.
/// Position is displacement from the reference configuration, not absolute length or a wrapped angle.
/// Rotational position includes accumulated turns. Positive velocity increases position.
struct ActuatorState {
    double position = 0;  ///< Reference displacement [m] or relative, unwrapped angle [rad].
    double velocity = 0;  ///< Displacement rate [m/s] or angular velocity [rad/s].
    double time = 0;      ///< Simulation time [s].
};

/// Output-side telemetry; power is in watts and effort in N or N*m.
/// Electrical power is positive when drawn from the bus, negative when returned.
/// Efficiency is mechanical/electrical power when both are positive, electrical/mechanical power when both
/// are negative, and zero otherwise (including zero power). It is bounded to [0, 1].
/// A stateful model can temporarily deliver stored energy: instantaneous efficiency is not an energy balance.
struct ActuatorTelemetry {
    double applied_effort = 0;       ///< Signed effort delivered to the mechanical system.
    double mechanical_power = 0;     ///< applied_effort * velocity; positive means mechanical output.
    double electrical_power = 0;     ///< Terminal voltage times winding current (electric model).
    double efficiency = 0;           ///< Instantaneous terminal conversion ratio, or table estimate.
    bool command_saturated = false;  ///< An enabled normalized command was clipped to [-1, 1].
    bool current_saturated = false;  ///< The current constraint is active.
    bool effort_saturated = false;   ///< An output effort constraint is active.
    bool voltage_saturated = false;  ///< Requested enabled voltage exceeded the bus voltage.
    double current = 0;              ///< Winding current [A]; zero if not provided by the model.
    double copper_loss = 0;          ///< Winding Joule loss i^2 R [W]; zero if not provided.
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

    /// Evaluate internal state derivatives without mutating the model or the supplied states.
    virtual void CalculateRHS(double, const std::vector<double>&, const ActuatorCommand&, const ActuatorState&, std::vector<double>& rhs) const { rhs.assign(GetNumStates(), 0); }

    /// Optionally provide d(rhs)/d(y), row-major, with GetNumStates() rows and columns.
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

    /// Compute telemetry using integrated states. Stateful models should override this along with EffortFromStates.
    /// The default preserves model-specific stateless telemetry and updates the mechanical output from states.
    virtual ActuatorTelemetry ComputeTelemetryFromStates(const std::vector<double>& y, const ActuatorCommand& command, const ActuatorState& state) const {
        auto telemetry = ComputeTelemetry(command, state);
        telemetry.applied_effort = EffortFromStates(y, command, state);
        telemetry.mechanical_power = telemetry.applied_effort * state.velocity;
        return telemetry;
    }
};

}  // namespace actuator_models
}  // namespace chrono

#endif
