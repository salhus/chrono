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
    actuator_models::ActuatorCommand command{GetInput(time), m_enabled};
    actuator_models::ActuatorState state{s - s_0, sd, time};
    if (m_model->GetNumStates() == 0)
        return m_model->ComputeEffort(command, state);

    if (static_cast<std::size_t>(GetStates().size()) != m_model->GetNumStates())
        throw std::runtime_error("Linear actuator internal states are not initialized");
    std::vector<double> y(GetStates().size());
    for (std::size_t i = 0; i < y.size(); ++i)
        y[i] = GetStates()(static_cast<Eigen::Index>(i));
    return m_model->EffortFromStates(y, command, state);
}

actuator_models::ActuatorTelemetry ChActuatorModelLinear::GetTelemetry(double time) const {
    if (static_cast<std::size_t>(GetStates().size()) != m_model->GetNumStates())
        throw std::runtime_error("Linear actuator internal states are not initialized");
    std::vector<double> states(static_cast<std::size_t>(GetStates().size()));
    for (std::size_t i = 0; i < states.size(); ++i)
        states[i] = GetStates()(static_cast<Eigen::Index>(i));
    return m_model->ComputeTelemetryFromStates(states, {GetInput(time), m_enabled}, {s - s_0, sd, time});
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
    if (states.size() != m_model->GetNumStates())
        throw std::runtime_error("Linear actuator model returned invalid initial-state dimensions");
    y0.resize(static_cast<Eigen::Index>(states.size()));
    for (std::size_t i = 0; i < states.size(); ++i)
        y0(static_cast<Eigen::Index>(i)) = states[i];
}

void ChActuatorModelLinear::CalculateRHS(double time, const ChVectorDynamic<>& y, ChVectorDynamic<>& rhs) {
    if (static_cast<std::size_t>(y.size()) != m_model->GetNumStates())
        throw std::runtime_error("Linear actuator model received invalid state dimensions");
    std::vector<double> states(static_cast<std::size_t>(y.size()));
    for (Eigen::Index i = 0; i < y.size(); ++i)
        states[static_cast<std::size_t>(i)] = y(i);
    std::vector<double> derivatives;
    actuator_models::ActuatorCommand command{GetInput(time), m_enabled};
    actuator_models::ActuatorState state{s - s_0, sd, time};
    m_model->CalculateRHS(time, states, command, state, derivatives);
    if (derivatives.size() != m_model->GetNumStates())
        throw std::runtime_error("Linear actuator model returned invalid RHS dimensions");
    rhs.resize(static_cast<Eigen::Index>(derivatives.size()));
    for (std::size_t i = 0; i < derivatives.size(); ++i)
        rhs(static_cast<Eigen::Index>(i)) = derivatives[i];
}

bool ChActuatorModelLinear::CalculateJac(double time, const ChVectorDynamic<>& y, const ChVectorDynamic<>&, ChMatrixDynamic<>& J) {
    if (static_cast<std::size_t>(y.size()) != m_model->GetNumStates())
        throw std::runtime_error("Linear actuator model received invalid state dimensions");
    std::vector<double> states(static_cast<std::size_t>(y.size()));
    for (Eigen::Index i = 0; i < y.size(); ++i)
        states[static_cast<std::size_t>(i)] = y(i);
    std::vector<std::vector<double>> jacobian;
    actuator_models::ActuatorCommand command{GetInput(time), m_enabled};
    actuator_models::ActuatorState state{s - s_0, sd, time};
    if (!m_model->CalculateJacobian(time, states, command, state, jacobian))
        return false;
    if (jacobian.size() != states.size())
        throw std::runtime_error("Linear actuator model returned invalid Jacobian row count");
    for (const auto& row : jacobian) {
        if (row.size() != states.size())
            throw std::runtime_error("Linear actuator model returned invalid Jacobian column count");
    }
    J.resize(y.size(), y.size());
    for (std::size_t i = 0; i < jacobian.size(); ++i) {
        for (std::size_t j = 0; j < jacobian[i].size(); ++j)
            J(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = jacobian[i][j];
    }
    return true;
}

}  // namespace chrono
