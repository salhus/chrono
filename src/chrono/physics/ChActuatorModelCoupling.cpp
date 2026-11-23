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

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "chrono/physics/ChActuatorModelCoupling.h"

namespace chrono {
namespace actuator_coupling {

namespace {

double Perturbation(double x) {
    return 1e-7 * std::max(1.0, std::abs(x));
}

void EvaluateRHS(const actuator_models::ChActuatorModel& model,
                 const std::vector<double>& y,
                 const actuator_models::ActuatorCommand& command,
                 const actuator_models::ActuatorState& state,
                 std::vector<double>& rhs) {
    model.CalculateRHS(state.time, y, command, state, rhs);
    if (rhs.size() != model.GetNumStates())
        throw std::runtime_error("Actuator model returned invalid RHS dimensions");
}

// Accumulate the contribution of one connected body (block index 0 or 1) if its variables are active.
template <typename Fn>
void ForActiveBodies(ChBody* body1, ChBody* body2, Fn&& fn) {
    if (body1->Variables().IsActive())
        fn(body1->Variables().GetOffset(), 0);
    if (body2->Variables().IsActive())
        fn(body2->Variables().GetOffset(), 6);
}

}  // namespace

bool ResolveCoupling(ChActuatorCouplingMode mode, const actuator_models::ChActuatorModel& model) {
    double inertia = model.GetReflectedInertia();
    if (!std::isfinite(inertia) || inertia < 0)
        throw std::invalid_argument("Actuator model reflected inertia must be finite and nonnegative");
    switch (mode) {
        case ChActuatorCouplingMode::ENABLED:
            return true;
        case ChActuatorCouplingMode::DISABLED:
            if (inertia > 0)
                throw std::invalid_argument("Actuator reflected inertia requires Newton-matrix coupling");
            return false;
        default:
            return model.IsStiff() || inertia > 0;
    }
}

std::vector<double> CopyStates(const ChVectorDynamic<>& states, std::size_t expected) {
    if (static_cast<std::size_t>(states.size()) != expected)
        throw std::runtime_error("Actuator internal states are not initialized");
    return std::vector<double>(states.data(), states.data() + states.size());
}

double EvaluateEffort(const actuator_models::ChActuatorModel& model,
                      const std::vector<double>& y,
                      const actuator_models::ActuatorCommand& command,
                      const actuator_models::ActuatorState& state) {
    return model.GetNumStates() == 0 ? model.ComputeEffort(command, state) : model.EffortFromStates(y, command, state);
}

ModelSensitivities ComputeSensitivities(const actuator_models::ChActuatorModel& model,
                                        const std::vector<double>& y,
                                        const actuator_models::ActuatorCommand& command,
                                        const actuator_models::ActuatorState& state) {
    const std::size_t n = model.GetNumStates();
    ModelSensitivities s;
    s.effort_states.assign(n, 0);
    s.rhs_position.assign(n, 0);
    s.rhs_velocity.assign(n, 0);

    std::vector<double> rhs_plus, rhs_minus;

    // Position and velocity sensitivities
    for (int k = 0; k < 2; ++k) {
        auto plus = state;
        auto minus = state;
        double& p = k == 0 ? plus.position : plus.velocity;
        double& m = k == 0 ? minus.position : minus.velocity;
        double h = Perturbation(p);
        p += h;
        m -= h;
        double dE = (EvaluateEffort(model, y, command, plus) - EvaluateEffort(model, y, command, minus)) / (2 * h);
        (k == 0 ? s.effort_position : s.effort_velocity) = dE;
        if (n > 0) {
            EvaluateRHS(model, y, command, plus, rhs_plus);
            EvaluateRHS(model, y, command, minus, rhs_minus);
            auto& target = k == 0 ? s.rhs_position : s.rhs_velocity;
            for (std::size_t i = 0; i < n; ++i)
                target[i] = (rhs_plus[i] - rhs_minus[i]) / (2 * h);
        }
    }

    // Internal-state sensitivities of the effort
    auto yp = y;
    for (std::size_t j = 0; j < n; ++j) {
        double h = Perturbation(y[j]);
        yp[j] = y[j] + h;
        double plus = EvaluateEffort(model, yp, command, state);
        yp[j] = y[j] - h;
        double minus = EvaluateEffort(model, yp, command, state);
        yp[j] = y[j];
        s.effort_states[j] = (plus - minus) / (2 * h);
    }

    return s;
}

void LoadCouplingKRM(ChMatrixRef H,
                     const ChVectorDynamic<>& g,
                     const ModelSensitivities& s,
                     double reflected_inertia,
                     double Kfactor,
                     double Rfactor,
                     double Mfactor) {
    const auto n = static_cast<Eigen::Index>(s.effort_states.size());
    if (H.rows() != 12 + n || H.cols() != 12 + n || g.size() != 12)
        throw std::runtime_error("Actuator coupling block has unexpected dimensions");

    H.setZero();

    // Body-body: reflected inertia, damping dE/dv and stiffness dE/dq along the actuator coordinate
    double coefficient = Mfactor * reflected_inertia - Rfactor * s.effort_velocity - Kfactor * s.effort_position;
    H.topLeftCorner(12, 12) = coefficient * g * g.transpose();

    // Body-state: effort sensitivity to internal states (states are velocity-level variables)
    for (Eigen::Index j = 0; j < n; ++j)
        H.block(0, 12 + j, 12, 1) = (-Rfactor * s.effort_states[j]) * g;

    // State-body: internal dynamics sensitivity to actuator kinematics
    for (Eigen::Index i = 0; i < n; ++i)
        H.block(12 + i, 0, 1, 12) = (-Rfactor * s.rhs_velocity[i] - Kfactor * s.rhs_position[i]) * g.transpose();
}

void LoadReflectedInertiaMv(ChBody* body1,
                            ChBody* body2,
                            const ChVectorDynamic<>& g,
                            double reflected_inertia,
                            ChVectorDynamic<>& R,
                            const ChVectorDynamic<>& w,
                            double c) {
    if (reflected_inertia == 0)
        return;
    double rate = 0;
    ForActiveBodies(body1, body2, [&](unsigned int offset, Eigen::Index block) { rate += g.segment(block, 6).dot(w.segment(offset, 6)); });
    ForActiveBodies(body1, body2, [&](unsigned int offset, Eigen::Index block) { R.segment(offset, 6) += (c * reflected_inertia * rate) * g.segment(block, 6); });
}

void LoadReflectedInertiaLumped(ChBody* body1,
                                ChBody* body2,
                                const ChVectorDynamic<>& g,
                                double reflected_inertia,
                                ChVectorDynamic<>& Md,
                                double& err,
                                double c) {
    if (reflected_inertia == 0)
        return;
    double sum_abs = 0;
    double sum_sq = 0;
    ForActiveBodies(body1, body2, [&](unsigned int offset, Eigen::Index block) {
        auto gb = g.segment(block, 6);
        Md.segment(offset, 6) += (c * reflected_inertia) * gb.cwiseAbs2();
        sum_abs += gb.cwiseAbs().sum();
        sum_sq += gb.squaredNorm();
    });
    err += c * reflected_inertia * (sum_abs * sum_abs - sum_sq);
}

}  // namespace actuator_coupling
}  // namespace chrono
