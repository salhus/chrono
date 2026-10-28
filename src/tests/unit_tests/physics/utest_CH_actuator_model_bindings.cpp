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

#include <memory>
#include <stdexcept>

#include "gtest/gtest.h"

#include "chrono/functions/ChFunctionConst.h"
#include "chrono/physics/ChActuatorModelLinear.h"
#include "chrono/physics/ChActuatorModelRotational.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/actuator_models/ChElectricActuatorModel.h"

using namespace chrono;
using namespace chrono::actuator_models;

namespace {

class PositionModel : public ChActuatorModel {
  public:
    double ComputeEffort(const ActuatorCommand&, const ActuatorState& state) const override { return state.position; }
};

class CommandModel : public ChActuatorModel {
  public:
    double ComputeEffort(const ActuatorCommand& command, const ActuatorState&) const override {
        return command.enabled ? command.value : -2;
    }
};

class InspectModel : public PositionModel {
  public:
    std::size_t GetNumStates() const override { return 1; }

    void CalculateRHS(double,
                      const std::vector<double>&,
                      const ActuatorCommand& command,
                      const ActuatorState& state,
                      std::vector<double>& rhs) const override {
        rhs = {state.position + state.velocity + (command.enabled ? command.value : -1)};
    }

    bool CalculateJacobian(double,
                           const std::vector<double>&,
                           const ActuatorCommand& command,
                           const ActuatorState& state,
                           std::vector<std::vector<double>>& jacobian) const override {
        jacobian = {{state.position + (command.enabled ? command.value : -1)}};
        return true;
    }
};

class MalformedModel : public InspectModel {
  public:
    enum class Output { INITIAL, RHS, JAC_ROWS, JAC_COLUMNS };
    explicit MalformedModel(Output output) : m_output(output) {}

    void SetInitialConditions(std::vector<double>& states) const override {
        states.assign(m_output == Output::INITIAL ? 2 : 1, 0);
    }

    void CalculateRHS(double,
                      const std::vector<double>&,
                      const ActuatorCommand&,
                      const ActuatorState&,
                      std::vector<double>& rhs) const override {
        rhs.assign(m_output == Output::RHS ? 2 : 1, 0);
    }

    bool CalculateJacobian(double,
                           const std::vector<double>&,
                           const ActuatorCommand&,
                           const ActuatorState&,
                           std::vector<std::vector<double>>& jacobian) const override {
        jacobian.assign(m_output == Output::JAC_ROWS ? 2 : 1,
                        std::vector<double>(m_output == Output::JAC_COLUMNS ? 2 : 1, 0));
        return true;
    }

  private:
    Output m_output;
};

template <class Binding>
class ExposedBinding : public Binding {
  public:
    using Binding::Binding;
    using Binding::CalculateJac;
    using Binding::CalculateRHS;
    using Binding::SetInitialConditions;

    void Refresh(double time) { this->Update(time, UpdateFlags::UPDATE_ALL); }
};

using LinearBinding = ExposedBinding<ChActuatorModelLinear>;
using RotationalBinding = ExposedBinding<ChActuatorModelRotational>;

template <class Binding>
void CheckCommandAndJacobian(Binding& actuator) {
    actuator.Initialize();
    ChVectorDynamic<> y(1);
    y.setZero();
    ChVectorDynamic<> rhs;
    ChMatrixDynamic<> jacobian;
    EXPECT_TRUE(actuator.IsEnabled());
    EXPECT_DOUBLE_EQ(actuator.GetInput(0), 0);
    actuator.CalculateRHS(0, y, rhs);
    EXPECT_DOUBLE_EQ(rhs(0), 0);
    actuator.SetInputFunction(chrono_types::make_shared<ChFunctionConst>(2));
    actuator.CalculateRHS(0, y, rhs);
    EXPECT_DOUBLE_EQ(rhs(0), 1);
    EXPECT_TRUE(actuator.CalculateJac(0, y, rhs, jacobian));
    ASSERT_EQ(jacobian.rows(), 1);
    ASSERT_EQ(jacobian.cols(), 1);
    EXPECT_DOUBLE_EQ(jacobian(0, 0), 1);
    actuator.SetEnabled(false);
    EXPECT_FALSE(actuator.IsEnabled());
    actuator.CalculateRHS(0, y, rhs);
    EXPECT_DOUBLE_EQ(rhs(0), -1);
    actuator.CalculateJac(0, y, rhs, jacobian);
    EXPECT_DOUBLE_EQ(jacobian(0, 0), -1);
    actuator.SetEnabled(true);
    actuator.SetInputFunction(nullptr);
    EXPECT_DOUBLE_EQ(actuator.GetInput(0), 0);
}

template <class Binding>
void CheckMalformedOutputs() {
    for (auto output : {MalformedModel::Output::INITIAL, MalformedModel::Output::RHS,
                        MalformedModel::Output::JAC_ROWS, MalformedModel::Output::JAC_COLUMNS}) {
        Binding actuator(std::make_shared<MalformedModel>(output));
        // Initialize with valid storage before testing the malformed model callbacks directly.
        if (output != MalformedModel::Output::INITIAL)
            actuator.Initialize();
        ChVectorDynamic<> y(1);
        y.setZero();
        ChVectorDynamic<> rhs;
        ChMatrixDynamic<> jacobian;
        if (output == MalformedModel::Output::INITIAL) {
            EXPECT_THROW(actuator.SetInitialConditions(y), std::runtime_error);
        } else if (output == MalformedModel::Output::RHS) {
            EXPECT_THROW(actuator.CalculateRHS(0, y, rhs), std::runtime_error);
        } else {
            EXPECT_THROW(actuator.CalculateJac(0, y, rhs, jacobian), std::runtime_error);
            EXPECT_EQ(jacobian.size(), 0);
        }
    }
}

template <class Binding>
void CheckIntegratedTelemetry() {
    ChElectricActuatorModel::Parameters parameters;
    parameters.inductance = 0.2;
    parameters.max_current = 100;
    auto model = std::make_shared<ChElectricActuatorModel>(parameters);
    auto actuator = std::make_shared<Binding>(model);
    actuator->SetInputFunction(chrono_types::make_shared<ChFunctionConst>(0.5));
    actuator->Initialize();
    const Binding& read_only = *actuator;
    EXPECT_DOUBLE_EQ(read_only.GetTelemetry(0).applied_effort, 0);
    EXPECT_DOUBLE_EQ(read_only.GetTelemetry(0).electrical_power, 0);
    ChSystemNSC system;
    system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT);
    system.SetSolverType(ChSolver::Type::SPARSE_QR);
    system.Add(actuator);
    system.DoStepDynamics(0.001);
    ASSERT_GT(actuator->GetStates()(0), 0);
    auto expected = model->ComputeTelemetryFromStates({actuator->GetStates()(0)}, {0.5, true}, {0, 0, 0.001});
    auto actual = read_only.GetTelemetry(0.001);
    EXPECT_DOUBLE_EQ(actual.applied_effort, expected.applied_effort);
    EXPECT_DOUBLE_EQ(actual.electrical_power, expected.electrical_power);
    EXPECT_LT(actual.electrical_power, model->ComputeTelemetry({0.5, true}, {}).electrical_power);
    actuator->SetEnabled(false);
    EXPECT_DOUBLE_EQ(read_only.GetTelemetry(0.001).electrical_power, 0);
}

}  // namespace

TEST(ActuatorModelBindings, ZeroInitialPositionAndReferences) {
    LinearBinding linear(std::make_shared<PositionModel>());
    RotationalBinding rotational(std::make_shared<PositionModel>());
    linear.Initialize();
    rotational.Initialize();
    EXPECT_DOUBLE_EQ(linear.GetActuatorForce(0), 0);
    EXPECT_DOUBLE_EQ(rotational.GetActuatorTorque(0), 0);
    linear.SetActuatorInitialLength(2);
    rotational.SetActuatorInitialAngle(3);
    EXPECT_DOUBLE_EQ(linear.GetActuatorForce(0), 0);
    EXPECT_DOUBLE_EQ(rotational.GetActuatorTorque(0), 0);
    linear.SetActuatorLength(2.5, 0.2);
    rotational.SetActuatorAngle(3 + 4 * CH_2PI, 0.3);
    EXPECT_DOUBLE_EQ(linear.GetActuatorForce(0), 0.5);
    EXPECT_NEAR(rotational.GetActuatorTorque(0), 4 * CH_2PI, 1e-12);
    EXPECT_DOUBLE_EQ(linear.GetTelemetry(0).mechanical_power, 0.1);
}

TEST(ActuatorModelBindings, CommandsAndEmptyJacobians) {
    LinearBinding linear(std::make_shared<InspectModel>());
    RotationalBinding rotational(std::make_shared<InspectModel>());
    CheckCommandAndJacobian(linear);
    CheckCommandAndJacobian(rotational);
}

TEST(ActuatorModelBindings, DisabledCommandReachesEffortAndTelemetry) {
    LinearBinding linear(std::make_shared<CommandModel>());
    RotationalBinding rotational(std::make_shared<CommandModel>());
    linear.Initialize();
    rotational.Initialize();
    EXPECT_DOUBLE_EQ(linear.GetActuatorForce(0), 0);
    EXPECT_DOUBLE_EQ(rotational.GetActuatorTorque(0), 0);
    linear.SetInputFunction(chrono_types::make_shared<ChFunctionConst>(0.5));
    rotational.SetInputFunction(chrono_types::make_shared<ChFunctionConst>(0.5));
    EXPECT_DOUBLE_EQ(linear.GetActuatorForce(0), 0.5);
    EXPECT_DOUBLE_EQ(rotational.GetActuatorTorque(0), 0.5);
    linear.SetEnabled(false);
    rotational.SetEnabled(false);
    EXPECT_DOUBLE_EQ(linear.GetActuatorForce(0), -2);
    EXPECT_DOUBLE_EQ(rotational.GetActuatorTorque(0), -2);
    EXPECT_DOUBLE_EQ(linear.GetTelemetry(0).applied_effort, -2);
    EXPECT_DOUBLE_EQ(rotational.GetTelemetry(0).applied_effort, -2);
}

TEST(ActuatorModelBindings, RejectMalformedDimensions) {
    CheckMalformedOutputs<LinearBinding>();
    CheckMalformedOutputs<RotationalBinding>();
}

TEST(ActuatorModelBindings, ConnectedLinearReferenceAndCurrentRHS) {
    auto body1 = std::make_shared<ChBody>();
    auto body2 = std::make_shared<ChBody>();
    body2->SetPos(ChVector3d(2, 0, 0));
    LinearBinding actuator(std::make_shared<InspectModel>());
    actuator.Connect(body1, body2, true, VNULL, VNULL);
    EXPECT_DOUBLE_EQ(actuator.GetActuatorForce(0), 0);
    body2->SetPos(ChVector3d(2.5, 0, 0));
    body2->SetPosDt(ChVector3d(0.2, 0, 0));
    actuator.Refresh(1);
    EXPECT_DOUBLE_EQ(actuator.GetActuatorForce(1), 0.5);
    EXPECT_NEAR(actuator.GetRHS()(0), 0.7, 1e-12);
}

TEST(ActuatorModelBindings, ConnectedMultipleTurnsBothDirections) {
    for (double direction : {-1.0, 1.0}) {
        auto body1 = std::make_shared<ChBody>();
        auto body2 = std::make_shared<ChBody>();
        RotationalBinding actuator(std::make_shared<InspectModel>());
        actuator.Connect(body1, body2, ChFrame<>());
        EXPECT_DOUBLE_EQ(actuator.GetActuatorTorque(0), 0);
        for (int step = 1; step <= 96; ++step) {
            double angle = direction * step * CH_PI / 8;
            body2->SetRot(QuatFromAngleZ(angle));
            actuator.Refresh(step);
            EXPECT_NEAR(actuator.GetActuatorTorque(step), angle, 1e-11);
            EXPECT_NEAR(actuator.GetRHS()(0), angle, 1e-11);
        }
    }
}

TEST(ActuatorModelBindings, TelemetryUsesIntegratedCurrent) {
    CheckIntegratedTelemetry<LinearBinding>();
    CheckIntegratedTelemetry<RotationalBinding>();
}
