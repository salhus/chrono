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
// PMSM/FOC and valve-controlled hydraulic cylinder actuator models: closed-form checks and Chrono bindings.
//
// =============================================================================

#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

#include "gtest/gtest.h"

#include "chrono/functions/ChFunctionConst.h"
#include "chrono/physics/ChActuatorModelLinear.h"
#include "chrono/physics/ChActuatorModelRotational.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/actuator_models/ChHydraulicCylinderModel.h"
#include "chrono/physics/actuator_models/ChPMSMActuatorModel.h"

using namespace chrono;
using namespace chrono::actuator_models;

namespace {

// Classical RK4 on the model states at fixed command and kinematics.
std::vector<double> Integrate(const ChActuatorModel& model, std::vector<double> y, const ActuatorCommand& c, const ActuatorState& s, double dt, double duration) {
    std::vector<double> k1, k2, k3, k4, tmp(y.size());
    auto axpy = [&](const std::vector<double>& k, double a) {
        for (std::size_t i = 0; i < y.size(); ++i)
            tmp[i] = y[i] + a * k[i];
        return tmp;
    };
    int steps = static_cast<int>(std::round(duration / dt));
    for (int n = 0; n < steps; ++n) {
        model.CalculateRHS(0, y, c, s, k1);
        model.CalculateRHS(0, axpy(k1, dt / 2), c, s, k2);
        model.CalculateRHS(0, axpy(k2, dt / 2), c, s, k3);
        model.CalculateRHS(0, axpy(k3, dt), c, s, k4);
        for (std::size_t i = 0; i < y.size(); ++i)
            y[i] += dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
    return y;
}

std::shared_ptr<ChBody> MakeBody(ChSystem& system, double mass, double inertia_z, bool fixed) {
    auto body = chrono_types::make_shared<ChBody>();
    body->SetMass(mass);
    body->SetInertiaXX(ChVector3d(1, 1, inertia_z));
    body->SetFixed(fixed);
    system.AddBody(body);
    return body;
}

void ConfigureSystem(ChSystem& system) {
    system.SetGravitationalAcceleration(VNULL);
    system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT);
    system.SetSolverType(ChSolver::Type::SPARSE_QR);
}

ChHydraulicCylinderModel::Parameters SymmetricCylinder() {
    ChHydraulicCylinderModel::Parameters p;
    p.annulus_area = p.cap_area;
    return p;
}

}  // namespace

// ----------------------------------------------------------------------------- PMSM

TEST(PMSMActuatorModel, ValidatesParameters) {
    EXPECT_NO_THROW(ChPMSMActuatorModel{});
    auto check = [](auto modify) {
        ChPMSMActuatorModel::Parameters p;
        modify(p);
        EXPECT_THROW(ChPMSMActuatorModel{p}, std::invalid_argument);
    };
    check([](auto& p) { p.pole_pairs = 0; });
    check([](auto& p) { p.inductance_q = 0; });
    check([](auto& p) { p.flux_linkage = -1; });
    check([](auto& p) { p.current_bandwidth = std::nan(""); });
    check([](auto& p) { p.motoring_efficiency = 1.1; });
    check([](auto& p) { p.rotor_inertia = -1; });
    ChPMSMActuatorModel model;
    EXPECT_EQ(model.GetNumStates(), 4u);
    EXPECT_TRUE(model.IsStiff());
    EXPECT_NEAR(model.GetTorqueConstant(), 0.12, 1e-15);
    EXPECT_NEAR(model.GetNoLoadSpeed(), 48 / std::sqrt(3.0) / 0.08, 1e-12);
    EXPECT_THROW(model.EffortFromStates({0, 0}, {}, {}), std::invalid_argument);
}

TEST(PMSMActuatorModel, ElectromechanicalPowerIdentity) {
    // Inverter power = copper loss + d(magnetic energy)/dt + T_motor * w_motor, including saliency and saturation.
    ChPMSMActuatorModel::Parameters p;
    p.inductance_d = 0.8e-3;
    p.inductance_q = 1.4e-3;
    ChPMSMActuatorModel model(p);
    const std::vector<std::vector<double>> states = {{0.3, 4, 0.01, -0.02}, {-2, -7, 0.1, 0.3}, {1, 9, -0.2, 0.5}};
    for (const auto& y : states) {
        for (double speed : {-400.0, -20.0, 0.0, 150.0, 500.0}) {
            for (ActuatorCommand c : {ActuatorCommand{0.7, true}, ActuatorCommand{-1.5, true}, ActuatorCommand{0.4, false}}) {
                ActuatorState s{0, speed, 0};
                std::vector<double> rhs;
                model.CalculateRHS(0, y, c, s, rhs);
                auto t = model.ComputeTelemetryFromStates(y, c, s);
                double magnetic_rate = 1.5 * (p.inductance_d * y[0] * rhs[0] + p.inductance_q * y[1] * rhs[1]);
                EXPECT_NEAR(t.electrical_power - t.copper_loss - magnetic_rate, t.mechanical_power, 1e-9 * (1 + std::abs(t.electrical_power)));
            }
        }
    }
}

TEST(PMSMActuatorModel, CurrentLoopIsFirstOrderAtBandwidth) {
    ChPMSMActuatorModel::Parameters p;
    p.inductance_q = 1.5e-3;
    ChPMSMActuatorModel model(p);
    ActuatorCommand c{0.5, true};
    ActuatorState s{0, 50, 0};
    double tau = 1 / p.current_bandwidth;
    auto y = Integrate(model, {0, 0, 0, 0}, c, s, tau / 200, tau);
    EXPECT_NEAR(y[ChPMSMActuatorModel::I_Q], 5 * (1 - std::exp(-1.0)), 1e-6);
    EXPECT_NEAR(y[ChPMSMActuatorModel::I_D], 0, 1e-9);
    y = Integrate(model, y, c, s, tau / 200, 20 * tau);
    EXPECT_NEAR(model.EffortFromStates(y, c, s), model.ComputeEffort(c, s), 1e-6);
    EXPECT_NEAR(model.ComputeEffort(c, s), 0.12 * 5, 1e-12);
    EXPECT_FALSE(model.ComputeTelemetry(c, s).voltage_saturated);
}

TEST(PMSMActuatorModel, VoltageLimitedTorqueSpeedCurve) {
    ChPMSMActuatorModel model;
    double no_load = model.GetNoLoadSpeed();
    EXPECT_NEAR(model.ComputeEffort({1, true}, {0, no_load, 0}), 0, 1e-12);
    EXPECT_LT(model.ComputeEffort({1, true}, {0, 1.1 * no_load, 0}), 0);
    EXPECT_NEAR(model.ComputeEffort({1, true}, {0, 0, 0}), model.GetPeakEffort(), 1e-12);
    double previous = model.GetPeakEffort();
    for (double fraction : {0.5, 0.8, 0.9, 0.95, 0.99}) {
        double effort = model.ComputeEffort({1, true}, {0, fraction * no_load, 0});
        EXPECT_LE(effort, previous + 1e-12);
        previous = effort;
    }
    EXPECT_TRUE(model.ComputeTelemetry({1, true}, {0, 0.99 * no_load, 0}).voltage_saturated);
}

TEST(PMSMActuatorModel, ShortCircuitBraking) {
    ChPMSMActuatorModel model;
    const auto& p = model.GetParameters();
    ActuatorCommand off{1, false};
    ActuatorState s{0, 100, 0};
    double we = 400;
    double det = p.resistance * p.resistance + we * we * p.inductance_d * p.inductance_q;
    double expected = 1.5 * p.pole_pairs * p.flux_linkage * (-p.resistance * we * p.flux_linkage / det);
    EXPECT_NEAR(model.ComputeEffort(off, s), expected, 1e-12);
    EXPECT_LT(expected, 0);
    auto y = Integrate(model, {3, 8, 0.01, 0.02}, off, s, 1e-6, 0.05);
    EXPECT_NEAR(model.EffortFromStates(y, off, s), expected, 1e-6);
    EXPECT_NEAR(y[ChPMSMActuatorModel::Z_Q], 0, 1e-12);
    EXPECT_DOUBLE_EQ(model.ComputeTelemetryFromStates(y, off, s).electrical_power, 0);
}

namespace {

struct SpinUpResult {
    double speed = 0;
    double relative_energy_residual = 0;
};

SpinUpResult PMSMSpinUp(double step, double duration) {
    ChPMSMActuatorModel::Parameters p;
    p.gear_ratio = 2;
    p.rotor_inertia = 2e-5;
    p.motoring_efficiency = 0.95;
    p.generating_efficiency = 0.9;
    auto model = std::make_shared<ChPMSMActuatorModel>(p);
    ChSystemNSC system;
    ConfigureSystem(system);
    auto ground = MakeBody(system, 1, 1, true);
    auto load = MakeBody(system, 1, 4e-3, false);
    auto actuator = chrono_types::make_shared<ChActuatorModelRotational>(model);
    actuator->SetInputFunction(chrono_types::make_shared<ChFunctionConst>(1));
    actuator->Connect(ground, load, ChFrame<>());
    system.Add(actuator);
    EXPECT_TRUE(actuator->IsCouplingActive());

    auto states = [&]() { return std::vector<double>(actuator->GetStates().data(), actuator->GetStates().data() + 4); };
    auto losses = [&]() {
        auto t = actuator->GetTelemetry(system.GetChTime());
        return std::array<double, 2>{t.electrical_power, t.copper_loss + t.gear_loss};
    };
    double input = 0, dissipated = 0;
    auto previous = losses();
    while (system.GetChTime() < duration - step / 2) {
        system.DoStepDynamics(step);
        auto current = losses();
        input += 0.5 * step * (previous[0] + current[0]);
        dissipated += 0.5 * step * (previous[1] + current[1]);
        previous = current;
    }
    double w = load->GetAngVelParent().z();
    double kinetic = 0.5 * (4e-3 + model->GetReflectedInertia()) * w * w;
    SpinUpResult r;
    r.speed = w;
    r.relative_energy_residual = (input - dissipated - kinetic - model->GetMagneticEnergy(states())) / input;
    return r;
}

}  // namespace

TEST(PMSMActuatorModel, ChronoSpinUpEnergyBalance) {
    ChPMSMActuatorModel::Parameters p;
    p.gear_ratio = 2;
    const double no_load = ChPMSMActuatorModel(p).GetNoLoadSpeed();  // output side
    // Current-limited acceleration, then a slow voltage-limited approach to the no-load speed.
    auto fine = PMSMSpinUp(1e-4, 1.5);
    EXPECT_NEAR(fine.speed, no_load, 5e-3 * no_load);
    EXPECT_LT(std::abs(fine.relative_energy_residual), 0.01);
    // Monolithic coupling keeps a step 4x the current-loop time constant stable and accurate at steady state.
    auto coarse = PMSMSpinUp(2e-3, 1.5);
    EXPECT_NEAR(coarse.speed, fine.speed, 0.01 * no_load);
    EXPECT_LT(std::abs(coarse.relative_energy_residual), 0.03);
}

// ----------------------------------------------------------------------------- Hydraulic cylinder

TEST(HydraulicCylinderModel, ValidatesParameters) {
    EXPECT_NO_THROW(ChHydraulicCylinderModel{});
    auto check = [](auto modify) {
        ChHydraulicCylinderModel::Parameters p;
        modify(p);
        EXPECT_THROW(ChHydraulicCylinderModel{p}, std::invalid_argument);
    };
    check([](auto& p) { p.cap_area = 0; });
    check([](auto& p) { p.initial_position = 0.6; });
    check([](auto& p) { p.bulk_modulus = -1; });
    check([](auto& p) { p.relief_pressure = 10e6; });
    check([](auto& p) { p.tank_pressure = 30e6; });
    check([](auto& p) { p.initial_pressure = -1; });
    check([](auto& p) { p.transition_pressure = 0; });
    ChHydraulicCylinderModel model;
    std::vector<double> y0;
    model.SetInitialConditions(y0);
    EXPECT_EQ(y0, std::vector<double>(2, 0));
    EXPECT_TRUE(model.IsStiff());
}

TEST(HydraulicCylinderModel, StallForcesAndOrifice) {
    ChHydraulicCylinderModel::Parameters p;
    p.tank_pressure = 0.3e6;
    p.initial_pressure = p.tank_pressure;
    ChHydraulicCylinderModel model(p);
    EXPECT_NEAR(model.ComputeEffort({1, true}, {}), 20e6 * p.cap_area - 0.3e6 * p.annulus_area, 1e-6);
    EXPECT_NEAR(model.ComputeEffort({-1, true}, {}), 0.3e6 * p.cap_area - 20e6 * p.annulus_area, 1e-6);
    EXPECT_NEAR(model.GetMaxExtendForce(), 20e6 * p.cap_area - 0.3e6 * p.annulus_area, 1e-6);
    EXPECT_NEAR(model.OrificeFlow(p.nominal_pressure_drop), p.nominal_flow, 1e-3 * p.nominal_flow);
    EXPECT_DOUBLE_EQ(model.OrificeFlow(-1e6), -model.OrificeFlow(1e6));
    // Closed valve: a moving load meets relief and tank pressures.
    EXPECT_NEAR(model.ComputeEffort({0, true}, {0, 0.1, 0}), 0.3e6 * p.cap_area - 25e6 * p.annulus_area, 1e-6);
    EXPECT_NEAR(model.ComputeEffort({0.8, false}, {0, -0.1, 0}), model.ComputeEffort({0, true}, {0, -0.1, 0}), 1e-9);
}

TEST(HydraulicCylinderModel, QuasiStaticNoLoadSpeed) {
    ChHydraulicCylinderModel model(SymmetricCylinder());
    const auto& p = model.GetParameters();
    // Bisection for zero effort.
    double lo = 0, hi = 2;
    for (int i = 0; i < 200; ++i) {
        double mid = 0.5 * (lo + hi);
        (model.ComputeEffort({1, true}, {0, mid, 0}) > 0 ? lo : hi) = mid;
    }
    double expected = model.OrificeFlow(0.5 * p.supply_pressure) / p.cap_area;
    EXPECT_NEAR(lo, expected, 1e-9);
    EXPECT_NEAR(expected, p.nominal_flow * std::sqrt(10e6 / 3.5e6) / p.cap_area, 1e-4 * expected);
}

TEST(HydraulicCylinderModel, PressureBoundsAndSteadyFlows) {
    ChHydraulicCylinderModel model;
    const auto& p = model.GetParameters();
    std::vector<double> rhs;
    // Closed valve, retracting load: cap chamber compressed at relief, rod chamber expanded at tank.
    model.CalculateRHS(0, {p.relief_pressure, p.tank_pressure}, {0, true}, {0, -0.05, 0}, rhs);
    EXPECT_DOUBLE_EQ(rhs[0], 0);  // relief blocks further rise
    EXPECT_DOUBLE_EQ(rhs[1], 0);  // anti-cavitation blocks further fall
    EXPECT_NEAR(model.EffortFromStates({30e6, -1e6}, {}, {}), p.relief_pressure * p.cap_area, 1e-6);
    // Quasi-static pressures are a dynamic equilibrium (leak-free).
    ActuatorCommand c{0.6, true};
    ActuatorState s{0.1, 0.12, 0};
    auto t = model.ComputeTelemetry(c, s);
    auto y = Integrate(model, {10e6, 5e6}, c, s, 1e-6, 0.3);
    EXPECT_NEAR(model.EffortFromStates(y, c, s), t.applied_effort, 1e-3 * std::abs(t.applied_effort));
    EXPECT_GT(t.electrical_power, t.mechanical_power);
}

TEST(HydraulicCylinderModel, AnalyticJacobianMatchesFiniteDifferences) {
    ChHydraulicCylinderModel::Parameters p;
    p.leakage_coefficient = 1e-12;
    ChHydraulicCylinderModel model(p);
    for (double u : {0.7, -0.4, 0.0}) {
        ActuatorCommand c{u, true};
        ActuatorState s{0.05, 0.03, 0};
        std::vector<double> y{12e6, 6e6};
        std::vector<std::vector<double>> J;
        ASSERT_TRUE(model.CalculateJacobian(0, y, c, s, J));
        for (int j = 0; j < 2; ++j) {
            double d = 1e-6 * y[j];
            auto yp = y, ym = y;
            yp[j] += d;
            ym[j] -= d;
            std::vector<double> rp, rm;
            model.CalculateRHS(0, yp, c, s, rp);
            model.CalculateRHS(0, ym, c, s, rm);
            for (int i = 0; i < 2; ++i) {
                double fd = (rp[i] - rm[i]) / (2 * d);
                EXPECT_NEAR(J[i][j], fd, 1e-5 * std::abs(J[i][i]) + 1e-6) << "u=" << u << " i=" << i << " j=" << j;
            }
        }
    }
}

// Large linearized implicit steps from a pressure bound (the default initial state is at tank pressure) must settle,
// not oscillate between tank and relief: the Jacobian keeps the interior slope on and past a bound.
TEST(HydraulicCylinderModel, LargeImplicitStepsFromBoundSettle) {
    ChHydraulicCylinderModel model;
    const auto& p = model.GetParameters();
    for (double h : {0.02, 0.05, 0.1}) {
        std::vector<double> y{p.tank_pressure, p.tank_pressure}, rhs;
        std::vector<std::vector<double>> J;
        for (int k = 0; k < 40; ++k) {
            model.CalculateRHS(0, y, {1, true}, {}, rhs);
            ASSERT_TRUE(model.CalculateJacobian(0, y, {1, true}, {}, J));
            // (I - h J) dy = h rhs
            double a = 1 - h * J[0][0], b = -h * J[0][1], c = -h * J[1][0], d = 1 - h * J[1][1];
            double det = a * d - b * c;
            y[0] += (d * h * rhs[0] - b * h * rhs[1]) / det;
            y[1] += (a * h * rhs[1] - c * h * rhs[0]) / det;
        }
        EXPECT_NEAR(y[0], p.supply_pressure, 1e-3 * p.supply_pressure) << "h=" << h;
        EXPECT_NEAR(y[1], p.tank_pressure, 1e-3 * p.supply_pressure) << "h=" << h;
    }
}

// Regression: ChExternalDynamicsODE must not add its identity mass twice (variables + KRM block). With the double
// count, Newton contracted by only ~1/2 per iteration and the fixed-iteration Euler step lagged the model rate by ~6%.
TEST(HydraulicCylinderModel, ImplicitStepAdvancesPressureAtModelRate) {
    for (auto type : {ChTimestepper::Type::EULER_IMPLICIT, ChTimestepper::Type::HHT}) {
        auto parameters = SymmetricCylinder();
        parameters.initial_pressure = 10e6;
        auto model = std::make_shared<ChHydraulicCylinderModel>(parameters);
        ChSystemNSC system;
        ConfigureSystem(system);
        system.SetTimestepperType(type);
        auto ground = MakeBody(system, 1, 1, true);
        auto rod = MakeBody(system, 100, 1, false);
        rod->SetPos(ChVector3d(1, 0, 0));
        rod->SetPosDt(ChVector3d(0.05, 0, 0));
        auto actuator = chrono_types::make_shared<ChActuatorModelLinear>(model);
        actuator->Connect(ground, rod, false, VNULL, ChVector3d(1, 0, 0));
        system.Add(actuator);
        std::vector<double> rhs;
        model->CalculateRHS(0, {10e6, 10e6}, {0, true}, {0, 0.05, 0}, rhs);
        const double h = 1e-6;
        system.DoStepDynamics(h);
        EXPECT_NEAR((actuator->GetStates()(0) - 10e6) / h, rhs[0], 1e-3 * std::abs(rhs[0]));
        EXPECT_NEAR((actuator->GetStates()(1) - 10e6) / h, rhs[1], 1e-3 * std::abs(rhs[1]));
    }
}

TEST(HydraulicCylinderModel, ChronoNoLoadExtension) {
    auto model = std::make_shared<ChHydraulicCylinderModel>(SymmetricCylinder());
    const auto& p = model->GetParameters();
    ChSystemNSC system;
    ConfigureSystem(system);
    auto ground = MakeBody(system, 1, 1, true);
    auto rod = MakeBody(system, 20, 1, false);
    rod->SetPos(ChVector3d(1, 0, 0));
    auto actuator = chrono_types::make_shared<ChActuatorModelLinear>(model);
    actuator->SetInputFunction(chrono_types::make_shared<ChFunctionConst>(1));
    actuator->Connect(ground, rod, false, VNULL, ChVector3d(1, 0, 0));
    system.Add(actuator);
    EXPECT_TRUE(actuator->IsCouplingActive());
    for (int i = 0; i < 200; ++i)
        system.DoStepDynamics(1e-3);
    double expected = model->OrificeFlow(0.5 * p.supply_pressure) / p.cap_area;
    EXPECT_NEAR(rod->GetPosDt().x(), expected, 5e-3 * expected);
    EXPECT_NEAR(actuator->GetStates()(0) + actuator->GetStates()(1), p.supply_pressure, 0.01 * p.supply_pressure);
}

namespace {

// Period of a mass on a locked (centered-spool) cylinder, from velocity zero crossings.
double LockedPeriod(double mass, double step, ChActuatorCouplingMode mode, double* peak_speed) {
    auto parameters = SymmetricCylinder();
    parameters.initial_pressure = 10e6;
    auto model = std::make_shared<ChHydraulicCylinderModel>(parameters);
    ChSystemNSC system;
    ConfigureSystem(system);
    auto ground = MakeBody(system, 1, 1, true);
    auto rod = MakeBody(system, mass, 1, false);
    rod->SetPos(ChVector3d(1, 0, 0));
    rod->SetPosDt(ChVector3d(0.05, 0, 0));
    auto actuator = chrono_types::make_shared<ChActuatorModelLinear>(model);
    actuator->SetCouplingMode(mode);
    actuator->Connect(ground, rod, false, VNULL, ChVector3d(1, 0, 0));
    system.Add(actuator);
    std::vector<double> crossings;
    double previous = rod->GetPosDt().x();
    *peak_speed = 0;
    while (crossings.size() < 3 && system.GetChTime() < 0.2) {
        system.DoStepDynamics(step);
        double v = rod->GetPosDt().x();
        *peak_speed = std::max(*peak_speed, std::abs(v));
        if (!std::isfinite(v))
            return NAN;
        if ((previous > 0) != (v > 0))
            crossings.push_back(system.GetChTime() - step * v / (v - previous));
        previous = v;
    }
    return crossings.size() == 3 ? crossings[2] - crossings[0] : NAN;
}

}  // namespace

TEST(HydraulicCylinderModel, ChronoLockedNaturalFrequency) {
    ChHydraulicCylinderModel model(SymmetricCylinder());
    const double mass = 100;
    double expected = 2 * CH_PI / std::sqrt(model.GetHydraulicStiffness(0) / mass);
    double peak = 0;
    EXPECT_NEAR(LockedPeriod(mass, 1e-5, ChActuatorCouplingMode::AUTO, &peak), expected, 0.01 * expected);
    EXPECT_LE(peak, 0.05 + 1e-9);
    // A step of an eighth of the period stays stable with monolithic coupling.
    double coarse = LockedPeriod(mass, expected / 8, ChActuatorCouplingMode::AUTO, &peak);
    EXPECT_TRUE(std::isfinite(coarse));
    EXPECT_LE(peak, 0.05 + 1e-9);
}
