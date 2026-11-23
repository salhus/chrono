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
//
// Shared Newton-matrix coupling and reflected inertia for actuator model bindings.
//
// =============================================================================

#ifndef CH_ACTUATOR_MODEL_COUPLING_H
#define CH_ACTUATOR_MODEL_COUPLING_H

#include <vector>

#include "chrono/core/ChMatrix.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/actuator_models/ChActuatorModel.h"

namespace chrono {

/// Newton-matrix coupling policy for actuator model bindings.
enum class ChActuatorCouplingMode {
    AUTO,     ///< couple if the model is stiff or has reflected inertia
    ENABLED,  ///< always couple (requires a solver that supports KRM blocks, e.g. a direct sparse solver)
    DISABLED  ///< explicit effort only (not allowed with reflected inertia)
};

namespace actuator_coupling {

/// Resolve a coupling mode for a model. Throws std::invalid_argument if reflected inertia is requested with DISABLED.
ChApi bool ResolveCoupling(ChActuatorCouplingMode mode, const actuator_models::ChActuatorModel& model);

/// Copy integrated ODE states into a model state vector.
ChApi std::vector<double> CopyStates(const ChVectorDynamic<>& states, std::size_t expected);

/// Sensitivities of a model's effort E(q, v, y) and internal-state derivatives f(q, v, y) at one evaluation point.
struct ModelSensitivities {
    double effort_position = 0;         ///< dE/dq
    double effort_velocity = 0;         ///< dE/dv
    std::vector<double> effort_states;  ///< dE/dy
    std::vector<double> rhs_position;   ///< df/dq
    std::vector<double> rhs_velocity;   ///< df/dv
};

/// Effort from the integrated states (stateful models) or the algebraic evaluation (stateless models).
ChApi double EvaluateEffort(const actuator_models::ChActuatorModel& model,
                            const std::vector<double>& y,
                            const actuator_models::ActuatorCommand& command,
                            const actuator_models::ActuatorState& state);

/// Central finite differences of the model's pure evaluations (no model state is mutated).
ChApi ModelSensitivities ComputeSensitivities(const actuator_models::ChActuatorModel& model,
                                              const std::vector<double>& y,
                                              const actuator_models::ActuatorCommand& command,
                                              const actuator_models::ActuatorState& state);

/// Load the coupling block ordered [body 1 (6), body 2 (6), internal states (n)]:
///   H = Mfactor * Jr * g g^T - Rfactor * dQ/dv - Kfactor * dQ/dx,
/// where Q = E * g on the bodies, f on the states, and g = d(rate)/d(body velocities).
/// The state-state block is left to ChExternalDynamicsODE. Changes of g with configuration are neglected.
ChApi void LoadCouplingKRM(ChMatrixRef H,
                           const ChVectorDynamic<>& g,
                           const ModelSensitivities& sensitivities,
                           double reflected_inertia,
                           double Kfactor,
                           double Rfactor,
                           double Mfactor);

/// R += c * Jr * g (g^T w), restricted to active bodies.
ChApi void LoadReflectedInertiaMv(ChBody* body1,
                                  ChBody* body2,
                                  const ChVectorDynamic<>& g,
                                  double reflected_inertia,
                                  ChVectorDynamic<>& R,
                                  const ChVectorDynamic<>& w,
                                  double c);

/// Md += c * Jr * diag(g g^T) on active bodies; err accumulates the neglected off-diagonal magnitude.
ChApi void LoadReflectedInertiaLumped(ChBody* body1,
                                      ChBody* body2,
                                      const ChVectorDynamic<>& g,
                                      double reflected_inertia,
                                      ChVectorDynamic<>& Md,
                                      double& err,
                                      double c);

}  // namespace actuator_coupling
}  // namespace chrono

#endif
