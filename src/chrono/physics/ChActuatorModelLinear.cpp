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

#include <stdexcept>
#include <utility>
#include <vector>

#include "chrono/physics/ChActuatorModelLinear.h"

namespace chrono {

ChActuatorModelLinear::ChActuatorModelLinear(std::shared_ptr<const actuator_models::ChActuatorModel> model) : m_model(std::move(model)) {
    if (!m_model)
        throw std::invalid_argument("Linear actuator model must not be null");
}

double ChActuatorModelLinear::GetActuatorForce(double time) {
    actuator_models::ActuatorCommand command{GetInput(time), true};
    actuator_models::ActuatorState state{s, sd, time};
    if (m_model->GetNumStates() == 0)
        return m_model->ComputeEffort(command, state);

    std::vector<double> y(GetStates().size());
    for (std::size_t i = 0; i < y.size(); ++i)
        y[i] = GetStates()(static_cast<Eigen::Index>(i));
    return m_model->EffortFromStates(y, command, state);
}

bool ChActuatorModelLinear::IsStiff() const {
    return m_model->IsStiff();
}

unsigned int ChActuatorModelLinear::GetNumStates() const {
    return static_cast<unsigned int>(m_model->GetNumStates());
}

void ChActuatorModelLinear::SetInitialConditions(ChVectorDynamic<>& y0) {
    std::vector<double> states;
    m_model->SetInitialConditions(states);
    y0.resize(static_cast<Eigen::Index>(states.size()));
    for (std::size_t i = 0; i < states.size(); ++i)
        y0(static_cast<Eigen::Index>(i)) = states[i];
}

void ChActuatorModelLinear::CalculateRHS(double time, const ChVectorDynamic<>& y, ChVectorDynamic<>& rhs) {
    std::vector<double> states(static_cast<std::size_t>(y.size()));
    for (Eigen::Index i = 0; i < y.size(); ++i)
        states[static_cast<std::size_t>(i)] = y(i);
    std::vector<double> derivatives;
    actuator_models::ActuatorCommand command{GetInput(time), true};
    actuator_models::ActuatorState state{s, sd, time};
    m_model->CalculateRHS(time, states, command, state, derivatives);
    rhs.resize(static_cast<Eigen::Index>(derivatives.size()));
    for (std::size_t i = 0; i < derivatives.size(); ++i)
        rhs(static_cast<Eigen::Index>(i)) = derivatives[i];
}

bool ChActuatorModelLinear::CalculateJac(double time, const ChVectorDynamic<>& y, const ChVectorDynamic<>&, ChMatrixDynamic<>& J) {
    std::vector<double> states(static_cast<std::size_t>(y.size()));
    for (Eigen::Index i = 0; i < y.size(); ++i)
        states[static_cast<std::size_t>(i)] = y(i);
    std::vector<std::vector<double>> jacobian;
    actuator_models::ActuatorCommand command{GetInput(time), true};
    actuator_models::ActuatorState state{s, sd, time};
    if (!m_model->CalculateJacobian(time, states, command, state, jacobian))
        return false;
    for (std::size_t i = 0; i < jacobian.size(); ++i) {
        for (std::size_t j = 0; j < jacobian[i].size(); ++j)
            J(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = jacobian[i][j];
    }
    return true;
}

}  // namespace chrono
