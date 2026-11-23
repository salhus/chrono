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
// Monolithic Newton-matrix coupling, reflected inertia, and reaction loads of actuator model bindings.
//
// =============================================================================

#include <cmath>
#include <memory>
#include <stdexcept>

#include "gtest/gtest.h"

#include "chrono/functions/ChFunctionConst.h"
#include "chrono/physics/ChActuatorModelCoupling.h"
#include "chrono/physics/ChActuatorModelLinear.h"
#include "chrono/physics/ChActuatorModelRotational.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/actuator_models/ChElectricActuatorModel.h"

using namespace chrono;
using namespace chrono::actuator_models;

namespace {

class ConstantEffortModel : public ChActuatorModel {
  public:
    ConstantEffortModel(double effort, double inertia) : m_effort(effort), m_inertia(inertia) {}
    double ComputeEffort(const ActuatorCommand&, const ActuatorState&) const override { return m_effort; }
    double GetReflectedInertia() const override { return m_inertia; }

  private:
    double m_effort;
    double m_inertia;
};

// E = 2 q + 3 v + 5 y, dy/dt = 7 q + 11 v - y
class AffineModel : public ChActuatorModel {
  public:
    double ComputeEffort(const ActuatorCommand&, const ActuatorState& s) const override { return 2 * s.position + 3 * s.velocity; }
    std::size_t GetNumStates() const override { return 1; }
    void CalculateRHS(double, const std::vector<double>& y, const ActuatorCommand&, const ActuatorState& s, std::vector<double>& rhs) const override {
        rhs = {7 * s.position + 11 * s.velocity - y[0]};
    }
    double EffortFromStates(const std::vector<double>& y, const ActuatorCommand& c, const ActuatorState& s) const override {
        return ComputeEffort(c, s) + 5 * y[0];
    }
};

ChElectricActuatorModel::Parameters DemoParameters(double inductance) {
    ChElectricActuatorModel::Parameters p;
    p.torque_constant = 0.2;
    p.back_emf_constant = 0.2;
    p.gear_ratio = 4;
    p.resistance = 2;
    p.bus_voltage = 24;
    p.motoring_efficiency = 0.9;
    p.max_current = 100;
    p.inductance = inductance;
    return p;
}

std::shared_ptr<ChBody> MakeBody(ChSystem& system, double inertia_z, bool fixed) {
    auto body = chrono_types::make_shared<ChBody>();
    body->SetMass(1);
    body->SetInertiaXX(ChVector3d(1, 1, inertia_z));
    body->SetFixed(fixed);
    system.AddBody(body);
    return body;
}

void ConfigureSystem(ChSystem& system) {
    system.SetGravitationalAcceleration(VNULL);
    // Velocity-dependent coupling terms need a Newton-type integrator; EULER_IMPLICIT_LINEARIZED applies them to the
    // absolute velocity and biases the equilibrium.
    system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT);
    system.SetSolverType(ChSolver::Type::SPARSE_QR);
}

// Drive a load of inertia J_load from rest with a rotational electric actuator; return the final output speed.
double SpinUp(double inductance, ChActuatorCouplingMode mode, double step, double duration, bool* active = nullptr) {
    ChSystemNSC system;
    ConfigureSystem(system);
    auto ground = MakeBody(system, 1, true);
    auto rotor = MakeBody(system, 0.001, false);
    auto actuator = chrono_types::make_shared<ChActuatorModelRotational>(std::make_shared<ChElectricActuatorModel>(DemoParameters(inductance)));
    actuator->SetCouplingMode(mode);
    actuator->SetInputFunction(chrono_types::make_shared<ChFunctionConst>(1));
    actuator->Connect(ground, rotor, ChFrame<>());
    system.Add(actuator);
    if (active)
        *active = actuator->IsCouplingActive();
    while (system.GetChTime() < duration - step / 2)
        system.DoStepDynamics(step);
    return rotor->GetAngVelParent().z();
}

}  // namespace

TEST(ActuatorCoupling, FiniteDifferenceSensitivitiesAndBlock) {
    AffineModel model;
    ActuatorCommand command{0.3, true};
    ActuatorState state{0.4, -0.5, 0};
    auto s = actuator_coupling::ComputeSensitivities(model, {0.6}, command, state);
    EXPECT_NEAR(s.effort_position, 2, 1e-6);
    EXPECT_NEAR(s.effort_velocity, 3, 1e-6);
    ASSERT_EQ(s.effort_states.size(), 1u);
    EXPECT_NEAR(s.effort_states[0], 5, 1e-6);
    EXPECT_NEAR(s.rhs_position[0], 7, 1e-6);
    EXPECT_NEAR(s.rhs_velocity[0], 11, 1e-6);

    ChVectorDynamic<> g(12);
    g.setZero();
    g(5) = -1;
    g(11) = 1;
    ChMatrixDynamic<> H(13, 13);
    double Kf = -0.01, Rf = -0.1, Mf = 1;
    actuator_coupling::LoadCouplingKRM(H, g, s, 0.25, Kf, Rf, Mf);
    double body = Mf * 0.25 - Rf * 3 - Kf * 2;
    EXPECT_NEAR(H(11, 11), body, 1e-6);
    EXPECT_NEAR(H(5, 11), -body, 1e-6);
    EXPECT_NEAR(H(11, 12), -Rf * 5, 1e-6);
    EXPECT_NEAR(H(5, 12), Rf * 5, 1e-6);
    EXPECT_NEAR(H(12, 11), -Rf * 11 - Kf * 7, 1e-6);
    EXPECT_DOUBLE_EQ(H(12, 12), 0);
    EXPECT_DOUBLE_EQ(H(0, 0), 0);

    ChMatrixDynamic<> wrong(12, 12);
    EXPECT_THROW(actuator_coupling::LoadCouplingKRM(wrong, g, s, 0, Kf, Rf, Mf), std::runtime_error);
}

TEST(ActuatorCoupling, ResolveModes) {
    ConstantEffortModel plain(1, 0), inertial(1, 0.5);
    ChElectricActuatorModel stiff(DemoParameters(1e-3));
    EXPECT_FALSE(actuator_coupling::ResolveCoupling(ChActuatorCouplingMode::AUTO, plain));
    EXPECT_TRUE(actuator_coupling::ResolveCoupling(ChActuatorCouplingMode::AUTO, inertial));
    EXPECT_TRUE(actuator_coupling::ResolveCoupling(ChActuatorCouplingMode::AUTO, stiff));
    EXPECT_TRUE(actuator_coupling::ResolveCoupling(ChActuatorCouplingMode::ENABLED, plain));
    EXPECT_FALSE(actuator_coupling::ResolveCoupling(ChActuatorCouplingMode::DISABLED, stiff));
    EXPECT_THROW(actuator_coupling::ResolveCoupling(ChActuatorCouplingMode::DISABLED, inertial), std::invalid_argument);
    ConstantEffortModel bad(1, -1);
    EXPECT_THROW(actuator_coupling::ResolveCoupling(ChActuatorCouplingMode::AUTO, bad), std::invalid_argument);
}

TEST(ActuatorCoupling, RotationalReactionConservesMomentum) {
    ChSystemNSC system;
    ConfigureSystem(system);
    auto body1 = MakeBody(system, 2, false);
    auto body2 = MakeBody(system, 2, false);
    auto actuator = chrono_types::make_shared<ChActuatorModelRotational>(std::make_shared<ConstantEffortModel>(1, 0));
    actuator->Connect(body1, body2, ChFrame<>());
    system.Add(actuator);
    for (int i = 0; i < 100; ++i)
        system.DoStepDynamics(0.01);
    double w1 = body1->GetAngVelParent().z();
    double w2 = body2->GetAngVelParent().z();
    EXPECT_NEAR(w1 + w2, 0, 1e-9);
    EXPECT_NEAR(w2, 0.5, 1e-6);
    EXPECT_NEAR(w1, -0.5, 1e-6);
}

TEST(ActuatorCoupling, RotationalReflectedInertia) {
    const double torque = 2, load = 0.5, reflected = 1.5, duration = 1;
    for (bool body1_fixed : {true, false}) {
        ChSystemNSC system;
        ConfigureSystem(system);
        auto body1 = MakeBody(system, load, body1_fixed);
        auto body2 = MakeBody(system, load, false);
        auto actuator = chrono_types::make_shared<ChActuatorModelRotational>(std::make_shared<ConstantEffortModel>(torque, reflected));
        actuator->Connect(body1, body2, ChFrame<>());
        system.Add(actuator);
        EXPECT_TRUE(actuator->IsCouplingActive());
        for (int i = 0; i < 100; ++i)
            system.DoStepDynamics(duration / 100);
        double w1 = body1->GetAngVelParent().z();
        double w2 = body2->GetAngVelParent().z();
        double reduced = body1_fixed ? load : load / 2;
        EXPECT_NEAR(w2 - w1, torque * duration / (reduced + reflected), 1e-9);
        if (!body1_fixed)
            EXPECT_NEAR(w1 + w2, 0, 1e-9);
    }
}

TEST(ActuatorCoupling, LinearReflectedMass) {
    const double force = 3, mass = 1, reflected = 2, duration = 1;
    ChSystemNSC system;
    ConfigureSystem(system);
    auto ground = MakeBody(system, 1, true);
    auto slider = MakeBody(system, 1, false);
    slider->SetPos(ChVector3d(1, 0, 0));
    auto actuator = chrono_types::make_shared<ChActuatorModelLinear>(std::make_shared<ConstantEffortModel>(force, reflected));
    actuator->Connect(ground, slider, false, VNULL, ChVector3d(1, 0, 0));
    system.Add(actuator);
    ASSERT_TRUE(actuator->IsCouplingActive());
    for (int i = 0; i < 100; ++i)
        system.DoStepDynamics(duration / 100);
    EXPECT_NEAR(slider->GetPosDt().x(), force * duration / (mass + reflected), 1e-9);
    EXPECT_NEAR(slider->GetPosDt().y(), 0, 1e-12);
}

TEST(ActuatorCoupling, DisabledWithReflectedInertiaThrows) {
    ChSystemNSC system;
    auto ground = MakeBody(system, 1, true);
    auto rotor = MakeBody(system, 1, false);
    ChActuatorModelRotational actuator(std::make_shared<ConstantEffortModel>(1, 0.1));
    actuator.SetCouplingMode(ChActuatorCouplingMode::DISABLED);
    EXPECT_THROW(actuator.Connect(ground, rotor, ChFrame<>()), std::invalid_argument);
}

TEST(ActuatorCoupling, ReconnectClearsCouplingBlock) {
    ChSystemNSC system;
    system.SetGravitationalAcceleration(VNULL);
    system.SetSolverType(ChSolver::Type::PSOR);  // PSOR rejects any KRM block
    auto ground = MakeBody(system, 1, true);
    auto rotor = MakeBody(system, 1, false);
    auto actuator = chrono_types::make_shared<ChActuatorModelRotational>(std::make_shared<ConstantEffortModel>(1, 0));
    actuator->SetCouplingMode(ChActuatorCouplingMode::ENABLED);
    actuator->Connect(ground, rotor, ChFrame<>());
    system.Add(actuator);
    EXPECT_TRUE(actuator->IsCouplingActive());
    EXPECT_THROW(system.DoStepDynamics(0.01), std::runtime_error);
    actuator->SetCouplingMode(ChActuatorCouplingMode::DISABLED);
    actuator->Connect(ground, rotor, ChFrame<>());
    EXPECT_FALSE(actuator->IsCouplingActive());
    EXPECT_NO_THROW(system.DoStepDynamics(0.01));
}

TEST(ActuatorCoupling, QuasiStaticMotorLargeStep) {
    const double no_load_speed = 24 / (0.2 * 4);
    bool active = true;
    double coupled = SpinUp(0, ChActuatorCouplingMode::ENABLED, 0.02, 1, &active);
    EXPECT_TRUE(active);
    EXPECT_NEAR(coupled, no_load_speed, 1e-3);
    double explicit_effort = SpinUp(0, ChActuatorCouplingMode::AUTO, 0.02, 1, &active);
    EXPECT_FALSE(active);
    EXPECT_FALSE(std::isfinite(explicit_effort) && std::abs(explicit_effort - no_load_speed) < 1);
}

TEST(ActuatorCoupling, StatefulMotorLargeStep) {
    const double no_load_speed = 24 / (0.2 * 4);
    bool active = false;
    double coupled = SpinUp(1e-4, ChActuatorCouplingMode::AUTO, 0.02, 1, &active);
    EXPECT_TRUE(active);
    EXPECT_NEAR(coupled, no_load_speed, 1e-3);
}

TEST(ActuatorCoupling, ElectricRotorInertiaAndRatedStall) {
    auto parameters = DemoParameters(0);
    parameters.rotor_inertia = 0.01;
    ChElectricActuatorModel model(parameters);
    EXPECT_DOUBLE_EQ(model.GetReflectedInertia(), 16 * 0.01);
    EXPECT_NEAR(model.GetRatedStallEffort(), 0.9 * 4 * 0.2 * 12, 1e-12);
    parameters.max_current = 10;
    EXPECT_NEAR(ChElectricActuatorModel(parameters).GetRatedStallEffort(), 0.9 * 4 * 0.2 * 10, 1e-12);
    parameters.max_effort = 5;
    parameters.coulomb_friction = 0.5;
    EXPECT_NEAR(ChElectricActuatorModel(parameters).GetRatedStallEffort(), 4.5, 1e-12);
    for (double bad : {-1.0, std::nan(""), double(INFINITY)}) {
        parameters.rotor_inertia = bad;
        EXPECT_THROW(ChElectricActuatorModel{parameters}, std::invalid_argument);
    }
}

TEST(ActuatorCoupling, ElectricRotorInertiaSpinUp) {
    // Quasi-static motoring branch: E(w) = 8.64 - 0.288 w, so w(t) = 30 (1 - exp(-0.288 t / (J_load + N^2 J_rotor))).
    auto parameters = DemoParameters(0);
    const double load = 0.001, duration = 0.005, step = 1e-5;
    for (double rotor_inertia : {0.0, load / 16}) {
        parameters.rotor_inertia = rotor_inertia;
        ChSystemNSC system;
        ConfigureSystem(system);
        auto ground = MakeBody(system, 1, true);
        auto body = MakeBody(system, load, false);
        auto actuator = chrono_types::make_shared<ChActuatorModelRotational>(std::make_shared<ChElectricActuatorModel>(parameters));
        actuator->SetCouplingMode(ChActuatorCouplingMode::ENABLED);
        actuator->SetInputFunction(chrono_types::make_shared<ChFunctionConst>(1));
        actuator->Connect(ground, body, ChFrame<>());
        system.Add(actuator);
        for (int i = 0; i < 500; ++i)
            system.DoStepDynamics(step);
        double inertia = load + 16 * rotor_inertia;
        EXPECT_NEAR(body->GetAngVelParent().z(), 30 * (1 - std::exp(-0.288 * duration / inertia)), 0.05);
    }
}
