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

#include <cmath>
#include <memory>

#include "gtest/gtest.h"

#include "chrono/functions/ChFunctionConst.h"
#include "chrono/physics/ChActuatorModelLinear.h"
#include "chrono/physics/ChActuatorModelRotational.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/actuator_models/ChElectricActuatorModel.h"
#include "chrono/physics/actuator_models/ChTabulatedActuatorModel.h"

using namespace chrono;
using namespace chrono::actuator_models;

namespace {

class CommandEffortModel : public ChActuatorModel {
  public:
    double ComputeEffort(const ActuatorCommand& command, const ActuatorState&) const override { return 10 * command.value; }
};

ChElectricActuatorModel::Parameters MakeElectricParameters() {
    ChElectricActuatorModel::Parameters parameters;
    parameters.inductance = 0;
    parameters.resistance = 2;
    parameters.torque_constant = 0.5;
    parameters.back_emf_constant = 0.2;
    parameters.bus_voltage = 12;
    parameters.gear_ratio = 3;
    parameters.motoring_efficiency = 0.8;
    parameters.generating_efficiency = 0.5;
    parameters.max_current = 10;
    return parameters;
}

}  // namespace

TEST(ActuatorModels, ElectricSteadyStateAndQuadrants) {
    ChElectricActuatorModel model(MakeElectricParameters());
    ActuatorState state{0, 0, 0};
    EXPECT_NEAR(model.ComputeEffort({0.5, true}, state), 3.6, 1e-12);

    auto parameters = MakeElectricParameters();
    parameters.gear_ratio = 2;
    parameters.torque_constant = 1;
    parameters.back_emf_constant = 1;
    parameters.resistance = 1;
    parameters.bus_voltage = 10;
    parameters.motoring_efficiency = 0.8;
    parameters.generating_efficiency = 0.5;
    ChElectricActuatorModel quadrant_model(parameters);
    state.velocity = 1;
    EXPECT_NEAR(quadrant_model.ComputeEffort({1, true}, state), 12.8, 1e-12);
    EXPECT_NEAR(quadrant_model.ComputeEffort({0, true}, state), -8, 1e-12);
}

TEST(ActuatorModels, ElectricVoltageAndCurrentLimitsAndNoLoadSpeed) {
    auto parameters = MakeElectricParameters();
    parameters.bus_voltage = 10;
    parameters.resistance = 1;
    parameters.max_current = 5;
    parameters.max_effort = 4;
    ChElectricActuatorModel model(parameters);

    auto telemetry = model.ComputeTelemetry({2, true}, {0, 0, 0});
    EXPECT_TRUE(telemetry.voltage_saturated);
    EXPECT_TRUE(telemetry.current_saturated);
    EXPECT_TRUE(telemetry.command_saturated);
    EXPECT_TRUE(telemetry.effort_saturated);
    EXPECT_NEAR(model.GetStallEffort(), 4, 1e-12);
    EXPECT_NEAR(model.GetNoLoadSpeed(), 10.0 / 0.2 / 3.0, 1e-12);
    EXPECT_NEAR(model.ComputeEffort({1, true}, {0, model.GetNoLoadSpeed(), 0}), 0, 1e-12);
}

TEST(ActuatorModels, ElectricStateJacobianMatchesFiniteDifference) {
    auto parameters = MakeElectricParameters();
    parameters.inductance = 0.2;
    parameters.resistance = 4;
    parameters.max_current = 100;
    ChElectricActuatorModel model(parameters);
    std::vector<std::vector<double>> jacobian;
    std::vector<double> y{0.7};
    std::vector<double> rhs;
    ActuatorCommand command{0.4, true};
    ActuatorState state{0, 0.1, 0};

    ASSERT_TRUE(model.CalculateJacobian(0, y, command, state, jacobian));
    constexpr double delta = 1e-6;
    std::vector<double> y_plus{y[0] + delta};
    std::vector<double> y_minus{y[0] - delta};
    std::vector<double> rhs_plus;
    std::vector<double> rhs_minus;
    model.CalculateRHS(0, y_plus, command, state, rhs_plus);
    model.CalculateRHS(0, y_minus, command, state, rhs_minus);
    EXPECT_NEAR(jacobian[0][0], (rhs_plus[0] - rhs_minus[0]) / (2 * delta), 1e-8);
    EXPECT_NEAR(jacobian[0][0], -20, 1e-12);
}

TEST(ActuatorModels, TabulatedInterpolationClampingAndPurity) {
    ChTabulatedActuatorModel model({0, 10}, {100, 50}, {0, 10}, {0, 100}, {0.5, 0.7, 0.6, 0.8}, {0.2, 0.4, 0.3, 0.5}, {0, 10}, {0, 4});
    ActuatorCommand command{0.5, true};
    ActuatorState state{0, 5, 0};
    double first = model.ComputeEffort(command, state);
    EXPECT_TRUE(ActuatorCommand{}.enabled);
    EXPECT_DOUBLE_EQ(first, model.ComputeEffort(command, state));
    EXPECT_NEAR(first, 35.5, 1e-12);
    EXPECT_NEAR(model.GetEfficiency(command, state), 0.621, 1e-12);

    state.velocity = -5;
    EXPECT_NEAR(model.GetEfficiency(command, state), 0.329, 1e-12);
    EXPECT_NEAR(model.ComputeEffort({0, true}, state), 2, 1e-12);
    state.velocity = 20;
    EXPECT_NEAR(model.ComputeEffort(command, state), 21, 1e-12);
}

TEST(ActuatorModels, TabulatedRejectsNonmonotonicAxes) {
    EXPECT_THROW(ChTabulatedActuatorModel({0, 0}, {1, 1}, {0, 1}, {0, 1}, {0.5, 0.5, 0.5, 0.5}, {0.5, 0.5, 0.5, 0.5}, {0, 1}, {0, 0}), std::invalid_argument);
}

TEST(ActuatorBindings, LinearActuatorDrivesAttachedBodyInExtensionDirection) {
    ChSystemNSC system;
    system.SetGravitationalAcceleration(VNULL);
    auto ground = chrono_types::make_shared<ChBody>();
    ground->SetFixed(true);
    system.AddBody(ground);
    auto body = chrono_types::make_shared<ChBody>();
    body->SetPos(ChVector3d(1, 0, 0));
    body->SetMass(1);
    body->SetInertiaXX(ChVector3d(1, 1, 1));
    system.AddBody(body);

    auto actuator = chrono_types::make_shared<ChActuatorModelLinear>(std::make_shared<CommandEffortModel>());
    actuator->SetInputFunction(chrono_types::make_shared<ChFunctionConst>(1));
    actuator->Connect(ground, body, true, VNULL, VNULL);
    system.Add(actuator);
    system.DoStepDynamics(1e-3);

    EXPECT_GT(body->GetPosDt().x(), 0);
}

TEST(ActuatorBindings, RotationalActuatorDrivesAttachedBodyInPositiveDirection) {
    ChSystemNSC system;
    system.SetGravitationalAcceleration(VNULL);
    auto ground = chrono_types::make_shared<ChBody>();
    ground->SetFixed(true);
    system.AddBody(ground);
    auto body = chrono_types::make_shared<ChBody>();
    body->SetMass(1);
    body->SetInertiaXX(ChVector3d(1, 1, 1));
    system.AddBody(body);

    auto actuator = chrono_types::make_shared<ChActuatorModelRotational>(std::make_shared<CommandEffortModel>());
    actuator->SetInputFunction(chrono_types::make_shared<ChFunctionConst>(1));
    actuator->Connect(ground, body, ChFrame<>(VNULL, QUNIT));
    system.Add(actuator);
    system.DoStepDynamics(1e-3);

    EXPECT_GT(body->GetAngVelParent().z(), 0);
}

TEST(ActuatorBindings, ElectricStateIntegratesWithImplicitStepper) {
    auto parameters = MakeElectricParameters();
    parameters.inductance = 0.05;
    parameters.resistance = 2;
    auto model = std::make_shared<ChElectricActuatorModel>(parameters);
    auto actuator = chrono_types::make_shared<ChActuatorModelRotational>(model);
    actuator->SetInputFunction(chrono_types::make_shared<ChFunctionConst>(1));
    actuator->SetActuatorInitialAngle(0);
    actuator->Initialize();

    ChSystemNSC system;
    system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT);
    system.SetSolverType(ChSolver::Type::SPARSE_QR);
    system.Add(actuator);
    system.DoStepDynamics(1e-3);

    ASSERT_EQ(actuator->GetStates().size(), 1);
    EXPECT_GT(actuator->GetStates()(0), 0);
}
