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

#ifndef CH_ACTUATOR_MODEL_LINEAR_H
#define CH_ACTUATOR_MODEL_LINEAR_H

#include <memory>

#include "chrono/physics/ChActuator.h"
#include "chrono/physics/actuator_models/ChActuatorModel.h"

namespace chrono {

/// Chrono binding for an engine-independent linear actuator model.
/// The model receives displacement from the initial length (s - s_0), extension speed, and evaluation time.
/// SetActuatorInitialLength() defines the reference for stand-alone use; Connect() uses the initial body separation.
/// Input is clamped to [-1,1] and defaults to zero. Disabling the command is distinct from disabling the physics item:
/// the model determines disabled-command behavior, including passive forces and internal-state dynamics.
/// Call Initialize() for stand-alone use, or Connect() for attached use, before evaluating a stateful model.
class ChApi ChActuatorModelLinear : public ChLinearActuator {
  public:
    /// Construct a binding sharing the specified immutable model. Throws std::invalid_argument for a null model.
    explicit ChActuatorModelLinear(std::shared_ptr<const actuator_models::ChActuatorModel> model);
    virtual ~ChActuatorModelLinear() {}

    /// Get the actuator model.
    std::shared_ptr<const actuator_models::ChActuatorModel> GetModel() const { return m_model; }

    /// Enable or disable the model command (enabled by default), without disabling the physics item.
    void SetEnabled(bool enabled) { m_enabled = enabled; }

    /// Return whether the model command is enabled.
    bool IsEnabled() const { return m_enabled; }

    /// Get force [N] using the current kinematics, integrated internal states, and input at the specified time.
    virtual double GetActuatorForce(double time) override;

    /// Get model telemetry using the current kinematics and integrated internal states at the specified time.
    /// This does not advance the actuator or the containing system.
    actuator_models::ActuatorTelemetry GetTelemetry(double time) const;

    /// Forward the model's solver stiffness declaration.
    virtual bool IsStiff() const override;

  protected:
    /// Return the model's internal-state dimension.
    virtual unsigned int GetNumStates() const override;
    /// Copy model initial conditions, validating the state dimension required by ChExternalDynamicsODE.
    virtual void SetInitialConditions(ChVectorDynamic<>& y0) override;
    /// Evaluate internal-state derivatives with the current actuator kinematics and command.
    virtual void CalculateRHS(double time, const ChVectorDynamic<>& y, ChVectorDynamic<>& rhs) override;
    /// Resize and load a validated square model Jacobian, or return false to request finite differences.
    virtual bool CalculateJac(double time,
                              const ChVectorDynamic<>& y,
                              const ChVectorDynamic<>& rhs,
                              ChMatrixDynamic<>& J) override;

  private:
    std::shared_ptr<const actuator_models::ChActuatorModel> m_model;  ///< shared immutable model
    bool m_enabled = true;                                         ///< model command enable flag
};

}  // namespace chrono

#endif
