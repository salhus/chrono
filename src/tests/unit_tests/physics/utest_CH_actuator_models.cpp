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
#include <limits>
#include <memory>
#include <utility>

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
    parameters.back_emf_constant = 0.5;
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
    // At stall, i = V/R and the zero-power gearbox factor is unity.
    EXPECT_NEAR(model.ComputeEffort({0.5, true}, state), 4.5, 1e-12);

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
    // The effort ceiling sets a conservative current cap Ieff=max_effort*eta_g/(N Kt).
    // At zero power the gearbox factor is unity, so stall effort is max_effort*eta_g.
    EXPECT_NEAR(model.GetStallEffort(), 2, 1e-12);
    EXPECT_NEAR(model.GetNoLoadSpeed(), 10.0 / 0.5 / 3.0, 1e-12);
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

TEST(ActuatorModels, ElectricSIConstantsAndSmoothingValidation) {
    auto parameters = MakeElectricParameters();
    // SI electromechanical conversion conserves power only when Kt and Ke agree.
    parameters.back_emf_constant = 0.2;
    EXPECT_THROW(ChElectricActuatorModel{parameters}, std::invalid_argument);
    parameters.allow_inconsistent_constants = true;
    EXPECT_NO_THROW(ChElectricActuatorModel{parameters});
    parameters.allow_inconsistent_constants = false;
    parameters.back_emf_constant = parameters.torque_constant * (1 + 0.5e-6);
    EXPECT_NO_THROW(ChElectricActuatorModel{parameters});
    parameters.back_emf_constant = parameters.torque_constant * (1 - 0.5e-6);
    EXPECT_NO_THROW(ChElectricActuatorModel{parameters});
    parameters.back_emf_constant = parameters.torque_constant * (1 + 2e-6);
    EXPECT_THROW(ChElectricActuatorModel{parameters}, std::invalid_argument);
    parameters.back_emf_constant = parameters.torque_constant;
    for (double invalid : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::quiet_NaN()}) {
        parameters.quadrant_smoothing = invalid;
        EXPECT_THROW(ChElectricActuatorModel{parameters}, std::invalid_argument);
        parameters.quadrant_smoothing = 1e-3;
        parameters.friction_smoothing_speed = invalid;
        EXPECT_THROW(ChElectricActuatorModel{parameters}, std::invalid_argument);
        parameters.friction_smoothing_speed = 1e-3;
    }
}

TEST(ActuatorModels, ElectricAnalyticCurrentStallAndNoLoad) {
    auto parameters = MakeElectricParameters();
    parameters.max_current = 100;
    ChElectricActuatorModel model(parameters);
    // Ohm's law: i = (u Vbus - Ke N omega)/R, including both signed quadrants.
    for (double command : {-1.0, -0.3, 0.0, 0.3, 1.0}) {
        for (double speed : {-4.0, 0.0, 4.0}) {
            double expected_current = (command * parameters.bus_voltage -
                                       parameters.back_emf_constant * parameters.gear_ratio * speed) /
                                      parameters.resistance;
            auto telemetry = model.ComputeTelemetry({command, true}, {0, speed, 0});
            EXPECT_NEAR(telemetry.current, expected_current, 1e-12);
            EXPECT_NEAR(telemetry.copper_loss, expected_current * expected_current * parameters.resistance, 1e-12);
        }
    }
    // Full-voltage stall is N Kt min(Vbus/R, Imax), not N Kt Imax unconditionally.
    EXPECT_NEAR(model.GetStallEffort(), 9, 1e-12);
    EXPECT_NEAR(model.GetStallEffort(), model.ComputeEffort({1, true}, {}), 1e-12);
    parameters.max_current = 2;
    ChElectricActuatorModel current_limited(parameters);
    EXPECT_NEAR(current_limited.GetStallEffort(), 3, 1e-12);
    for (double sign : {-1.0, 1.0}) {
        auto telemetry = model.ComputeTelemetry({sign, true}, {0, sign * model.GetNoLoadSpeed(), 0});
        EXPECT_NEAR(telemetry.current, 0, 1e-12);
        EXPECT_NEAR(telemetry.applied_effort, 0, 1e-12);
        EXPECT_NEAR(telemetry.electrical_power, 0, 1e-12);
    }
}

TEST(ActuatorModels, ElectricCompactQuadrantSmoothingIsC1AndPassive) {
    auto parameters = MakeElectricParameters();
    parameters.inductance = 0.2;
    parameters.quadrant_smoothing = 2;
    ChElectricActuatorModel model(parameters);
    auto effort = [&](double current, double speed) {
        return model.EffortFromStates({current}, {}, {0, speed, 0});
    };
    // The compact cubic h(q)=q^2(3-2q) joins unity at p=0 to the exact far-quadrant factors.
    for (double current : {-2.0, 2.0}) {
        for (double normalized_power : {-3.0, -1.0, -0.5, -0.01, 0.0, 0.01, 0.5, 1.0, 3.0}) {
            double speed = normalized_power * parameters.quadrant_smoothing /
                           (parameters.gear_ratio * parameters.torque_constant * current);
            double power = parameters.gear_ratio * parameters.torque_constant * current * speed;
            double q = std::min(std::abs(normalized_power), 1.0);
            double h = q * q * (3 - 2 * q);
            double factor = power >= 0 ? 1 + (parameters.motoring_efficiency - 1) * h
                                       : 1 + (1 / parameters.generating_efficiency - 1) * h;
            double actual = effort(current, speed);
            EXPECT_NEAR(actual, parameters.gear_ratio * parameters.torque_constant * current * factor, 1e-12);
            EXPECT_GE(power - actual * speed, -1e-12);
        }
    }
    // h'(0)=h'(1)=0 guarantees matching first derivatives at all three branch joins,
    // whether electromagnetic power crosses the join by changing current or output speed.
    constexpr double delta = 1e-6;
    for (bool vary_current : {false, true}) {
        for (double join : {-1.0, 0.0, 1.0}) {
            double crossover = join * parameters.quadrant_smoothing /
                               (parameters.gear_ratio * parameters.torque_constant);
            auto sample = [&](double x) { return vary_current ? effort(x, 1) : effort(1, x); };
            double left = (sample(crossover) - sample(crossover - delta)) / delta;
            double right = (sample(crossover + delta) - sample(crossover)) / delta;
            EXPECT_NEAR(left, right, 2e-5);
        }
    }
    EXPECT_DOUBLE_EQ(effort(2, 0), 3);
    EXPECT_DOUBLE_EQ(effort(0, 2), 0);
}

TEST(ActuatorModels, ElectricCurrentClampJacobianAndInwardRelease) {
    auto parameters = MakeElectricParameters();
    parameters.inductance = 0.2;
    parameters.max_current = 2;
    ChElectricActuatorModel model(parameters);
    // An outward active current constraint has di/dt=0 and Jacobian zero;
    // an inward voltage releases it and recovers the winding slope -R/L.
    for (double sign : {-1.0, 1.0}) {
        for (bool outward : {false, true}) {
            ActuatorCommand command{(outward ? sign : -sign), true};
            std::vector<double> current{sign * parameters.max_current};
            std::vector<double> rhs;
            std::vector<std::vector<double>> jacobian;
            model.CalculateRHS(0, current, command, {}, rhs);
            ASSERT_TRUE(model.CalculateJacobian(0, current, command, {}, jacobian));
            ASSERT_EQ(rhs.size(), 1);
            ASSERT_EQ(jacobian.size(), 1);
            ASSERT_EQ(jacobian[0].size(), 1);
            if (outward) {
                EXPECT_DOUBLE_EQ(rhs[0], 0);
                EXPECT_DOUBLE_EQ(jacobian[0][0], 0);
            } else {
                EXPECT_LT(sign * rhs[0], 0);
                EXPECT_DOUBLE_EQ(jacobian[0][0], -parameters.resistance / parameters.inductance);
            }
        }
    }
}

TEST(ActuatorModels, ElectricTelemetryFlagsAndDisabledShortCircuit) {
    auto parameters = MakeElectricParameters();
    parameters.inductance = 0.2;
    parameters.max_current = 2;
    ChElectricActuatorModel model(parameters);
    // Command clipping and voltage clipping describe only enabled |u|>1, not a current constraint.
    for (double command : {-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0}) {
        for (bool enabled : {false, true}) {
            auto telemetry = model.ComputeTelemetryFromStates({0.5}, {command, enabled}, {});
            EXPECT_EQ(telemetry.command_saturated, enabled && std::abs(command) > 1);
            EXPECT_EQ(telemetry.voltage_saturated, enabled && std::abs(command) > 1);
            EXPECT_FALSE(telemetry.current_saturated);
            EXPECT_FALSE(telemetry.effort_saturated);
            EXPECT_DOUBLE_EQ(telemetry.current, 0.5);
            EXPECT_DOUBLE_EQ(telemetry.copper_loss, 0.5);
            double voltage = enabled ? parameters.bus_voltage * std::clamp(command, -1.0, 1.0) : 0;
            EXPECT_DOUBLE_EQ(telemetry.electrical_power, voltage * telemetry.current);
        }
    }
    for (double current : {-3.0, -2.0, -1.9, 1.9, 2.0, 3.0}) {
        auto telemetry = model.ComputeTelemetryFromStates({current}, {0, true}, {});
        EXPECT_EQ(telemetry.current_saturated, std::abs(current) >= parameters.max_current);
        EXPECT_DOUBLE_EQ(telemetry.current, std::clamp(current, -2.0, 2.0));
        EXPECT_FALSE(telemetry.command_saturated);
        EXPECT_FALSE(telemetry.voltage_saturated);
        EXPECT_FALSE(telemetry.effort_saturated);
    }
    // Disable shorts the winding (V=0): back EMF still drives braking current and copper dissipation.
    for (double speed : {-1.0, 1.0}) {
        auto disabled = model.ComputeTelemetry({5, false}, {0, speed, 0});
        auto zero_voltage = model.ComputeTelemetry({0, true}, {0, speed, 0});
        EXPECT_DOUBLE_EQ(disabled.current, zero_voltage.current);
        EXPECT_DOUBLE_EQ(disabled.applied_effort, zero_voltage.applied_effort);
        EXPECT_DOUBLE_EQ(disabled.electrical_power, 0);
        EXPECT_LT(disabled.mechanical_power, 0);
        EXPECT_GT(disabled.copper_loss, 0);
        EXPECT_FALSE(disabled.command_saturated);
        EXPECT_FALSE(disabled.voltage_saturated);
    }
    // Back-driving beyond no-load speed returns power to an enabled bus, unlike a shorted winding.
    auto generating = model.ComputeTelemetry({1, true}, {0, 10, 0});
    EXPECT_LT(generating.electrical_power, 0);
    EXPECT_LT(generating.mechanical_power, 0);
    // Motoring reports Pmech/Pbus; generating reports Pbus/Pmech, in either signed direction.
    for (double sign : {-1.0, 1.0}) {
        auto motoring = model.ComputeTelemetryFromStates({sign}, {sign, true}, {0, sign, 0});
        EXPECT_NEAR(motoring.mechanical_power, 1.2, 1e-12);
        EXPECT_NEAR(motoring.electrical_power, 12, 1e-12);
        EXPECT_NEAR(motoring.efficiency, 0.1, 1e-12);
        auto regenerating = model.ComputeTelemetryFromStates({-sign}, {0.1 * sign, true}, {0, sign, 0});
        EXPECT_NEAR(regenerating.mechanical_power, -3, 1e-12);
        EXPECT_NEAR(regenerating.electrical_power, -1.2, 1e-12);
        EXPECT_NEAR(regenerating.efficiency, 0.4, 1e-12);
        auto mixed_power = model.ComputeTelemetryFromStates({sign}, {-sign, true}, {0, sign, 0});
        EXPECT_DOUBLE_EQ(mixed_power.efficiency, 0);
    }
}

TEST(ActuatorModels, ElectricFourQuadrantEnergyBalance) {
    auto parameters = MakeElectricParameters();
    parameters.max_current = 2;
    parameters.quadrant_smoothing = 1;
    parameters.coulomb_friction = 0.3;
    parameters.viscous_friction = 0.1;
    parameters.friction_smoothing_speed = 0.2;
    // Kirchhoff plus electromagnetic conversion gives:
    // V i - tau omega = R i^2 + L i di/dt + gearbox_loss + friction_loss + limiter_loss.
    for (auto configuration : {std::make_pair(0.0, std::numeric_limits<double>::infinity()),
                               std::make_pair(0.2, std::numeric_limits<double>::infinity()),
                               std::make_pair(0.0, 1.5), std::make_pair(0.2, 1.5)}) {
        double inductance = configuration.first;
        double max_effort = configuration.second;
        parameters.inductance = inductance;
        parameters.max_effort = max_effort;
        double current_limit = std::min(parameters.max_current,
                                        max_effort * parameters.generating_efficiency /
                                            (parameters.gear_ratio * parameters.torque_constant));
        ChElectricActuatorModel model(parameters);
        for (double command : {-2.0, -1.0, -0.1, 0.0, 0.1, 1.0, 2.0}) {
            for (double speed : {-12.0, -2.0, -0.001, 0.0, 0.001, 2.0, 12.0}) {
                for (bool enabled : {false, true}) {
                    for (double fraction : {-1.0, -0.5, -0.005, 0.0, 0.005, 0.5, 1.0}) {
                        double initial_current = fraction * current_limit;
                        ActuatorCommand input{command, enabled};
                        ActuatorState state{0, speed, 0};
                        auto telemetry = inductance > 0
                                             ? model.ComputeTelemetryFromStates({initial_current}, input, state)
                                             : model.ComputeTelemetry(input, state);
                        double current = telemetry.current;
                        double voltage = enabled ? parameters.bus_voltage * std::clamp(command, -1.0, 1.0) : 0;
                        double electromagnetic_power =
                            parameters.gear_ratio * parameters.torque_constant * current * speed;
                        double friction = parameters.coulomb_friction *
                                              std::tanh(speed / parameters.friction_smoothing_speed) +
                                          parameters.viscous_friction * speed;
                        double friction_loss = friction * speed;
                        double gearbox_loss = electromagnetic_power -
                                              (telemetry.applied_effort + friction) * speed;
                        double current_rate = 0;
                        if (inductance > 0) {
                            std::vector<double> rhs;
                            model.CalculateRHS(0, {initial_current}, input, state, rhs);
                            current_rate = rhs[0];
                        }
                        double stored_power = inductance * current * current_rate;
                        double limiter_loss =
                            (voltage - parameters.resistance * current -
                             parameters.back_emf_constant * parameters.gear_ratio * speed -
                             inductance * current_rate) *
                            current;
                        SCOPED_TRACE(::testing::Message() << "max_effort=" << max_effort
                                                         << " L=" << inductance << " u=" << command
                                                         << " enabled=" << enabled << " w=" << speed
                                                         << " i=" << current);
                        EXPECT_GE(telemetry.copper_loss, 0);
                        EXPECT_GE(friction_loss, -1e-12);
                        EXPECT_GE(gearbox_loss, -1e-12);
                        EXPECT_GE(limiter_loss, -1e-12);
                        EXPECT_NEAR(telemetry.electrical_power - telemetry.mechanical_power - telemetry.copper_loss -
                                        stored_power - gearbox_loss - friction_loss - limiter_loss,
                                    0, 2e-12);
                        EXPECT_GE(telemetry.efficiency, 0);
                        EXPECT_LE(telemetry.efficiency, 1);
                }
            }
        }
      }
    }
}

TEST(ActuatorModels, ElectricEffortCeilingUsesPassiveCurrentLimiting) {
    auto parameters = MakeElectricParameters();
    parameters.inductance = 0.2;
    parameters.max_effort = 3;
    parameters.coulomb_friction = 0.4;
    parameters.viscous_friction = 1;
    ChElectricActuatorModel model(parameters);
    double current_limit = parameters.max_effort * parameters.generating_efficiency /
                           (parameters.gear_ratio * parameters.torque_constant);
    // A clipped command need not activate the effort constraint: integrated current is independent of voltage.
    auto command_only = model.ComputeTelemetryFromStates({0}, {2, true}, {});
    EXPECT_TRUE(command_only.command_saturated);
    EXPECT_TRUE(command_only.voltage_saturated);
    EXPECT_FALSE(command_only.effort_saturated);
    // Conversely, an effort-current constraint can be active without clipping the command or terminal voltage.
    auto effort_only = model.ComputeTelemetryFromStates({2 * current_limit}, {0, true}, {});
    EXPECT_FALSE(effort_only.command_saturated);
    EXPECT_FALSE(effort_only.voltage_saturated);
    EXPECT_TRUE(effort_only.effort_saturated);
    EXPECT_NEAR(effort_only.current, current_limit, 1e-12);
    EXPECT_NEAR(effort_only.applied_effort, parameters.max_effort * parameters.generating_efficiency, 1e-12);
    // At |omega| above no-load speed, the finite current cap retains regenerative braking and bus return.
    for (double sign : {-1.0, 1.0}) {
        double speed = sign * 2 * model.GetNoLoadSpeed();
        auto generating = model.ComputeTelemetry({sign, true}, {0, speed, 0});
        EXPECT_TRUE(generating.effort_saturated);
        EXPECT_NEAR(generating.current, -sign * current_limit, 1e-12);
        EXPECT_LT(generating.electrical_power, 0);
        EXPECT_LT(generating.mechanical_power, 0);
        EXPECT_GE(generating.electrical_power - generating.mechanical_power - generating.copper_loss, 0);
        double friction = parameters.coulomb_friction *
                              std::tanh(speed / parameters.friction_smoothing_speed) +
                          parameters.viscous_friction * speed;
        EXPECT_NEAR(std::abs(generating.applied_effort + friction), parameters.max_effort, 1e-12);
    }
    // The electromagnetic ceiling is enforced through the winding equation, never by erasing regenerative
    // braking effort. Passive friction remains additional and may exceed the electromagnetic rating.
    for (double speed : {-10.0, 10.0}) {
        auto telemetry = model.ComputeTelemetryFromStates({0}, {0, true}, {0, speed, 0});
        EXPECT_FALSE(telemetry.effort_saturated);
        EXPECT_GT(std::abs(telemetry.applied_effort), parameters.max_effort);
        EXPECT_LT(telemetry.mechanical_power, 0);
    }
    // The same fixed effort-current boundary must block outward state evolution and release inward motion.
    for (double sign : {-1.0, 1.0}) {
        for (bool outward : {false, true}) {
            std::vector<double> rhs;
            std::vector<std::vector<double>> jacobian;
            ActuatorCommand command{outward ? sign : -sign, true};
            model.CalculateRHS(0, {sign * current_limit}, command, {}, rhs);
            ASSERT_TRUE(model.CalculateJacobian(0, {sign * current_limit}, command, {}, jacobian));
            EXPECT_DOUBLE_EQ(jacobian[0][0], outward ? 0 : -parameters.resistance / parameters.inductance);
            if (outward)
                EXPECT_DOUBLE_EQ(rhs[0], 0);
            else
                EXPECT_LT(sign * rhs[0], 0);
        }
    }
}

TEST(ActuatorModels, ElectricAnalyticRotorIntegrationAndIntegratedEnergy) {
    auto parameters = MakeElectricParameters();
    parameters.max_current = 100;
    parameters.quadrant_smoothing = 1e-6;
    ChElectricActuatorModel model(parameters);
    constexpr double inertia = 2;
    double no_load = model.GetNoLoadSpeed();
    // Away from crossover J omega_dot = N Kt eta_m (V-Ke N omega)/R in motoring,
    // or N Kt (V-Ke N omega)/(R eta_g) while generating.
    for (bool generating : {false, true}) {
        double initial_speed = generating ? 1.5 * no_load : 0.25 * no_load;
        double time_constant = inertia * parameters.resistance /
                               (parameters.gear_ratio * parameters.gear_ratio *
                                parameters.torque_constant * parameters.back_emf_constant) *
                               (generating ? parameters.generating_efficiency : 1 / parameters.motoring_efficiency);
        double step = time_constant / 2000;
        double speed = initial_speed;
        double bus_energy = 0;
        double dissipated_energy = 0;
        auto rhs = [&](double velocity) {
            auto telemetry = model.ComputeTelemetry({1, true}, {0, velocity, 0});
            return telemetry.applied_effort / inertia;
        };
        auto powers = [&](double velocity) {
            auto telemetry = model.ComputeTelemetry({1, true}, {0, velocity, 0});
            double motor_power = parameters.gear_ratio * parameters.torque_constant * telemetry.current * velocity;
            EXPECT_GT(std::abs(motor_power), parameters.quadrant_smoothing);
            double gear_loss = motor_power - telemetry.mechanical_power;
            EXPECT_GE(gear_loss, 0);
            return std::pair<double, double>{telemetry.electrical_power, telemetry.copper_loss + gear_loss};
        };
        // RK4 integrates both J omega_dot=tau and the bus/loss energy quadratures.
        for (int i = 0; i < 4000; ++i) {
            double k1 = rhs(speed);
            double v2 = speed + step * k1 / 2;
            double k2 = rhs(v2);
            double v3 = speed + step * k2 / 2;
            double k3 = rhs(v3);
            double v4 = speed + step * k3;
            double k4 = rhs(v4);
            auto p1 = powers(speed);
            auto p2 = powers(v2);
            auto p3 = powers(v3);
            auto p4 = powers(v4);
            bus_energy += step * (p1.first + 2 * p2.first + 2 * p3.first + p4.first) / 6;
            dissipated_energy += step * (p1.second + 2 * p2.second + 2 * p3.second + p4.second) / 6;
            speed += step * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
            double exact = no_load + (initial_speed - no_load) * std::exp(-(i + 1) * step / time_constant);
            EXPECT_NEAR(speed, exact, 2e-12);
        }
        // Integrated bus work equals kinetic-energy gain plus nonnegative Joule and gear dissipation.
        double kinetic_energy_change = inertia * (speed * speed - initial_speed * initial_speed) / 2;
        EXPECT_GE(dissipated_energy, 0);
        EXPECT_NEAR(bus_energy - kinetic_energy_change - dissipated_energy, 0, 2e-9);
    }
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

TEST(ActuatorModels, TabulatedFrictionSmoothingValidationAndC1Passivity) {
    auto make_model = [](double smoothing) {
        return ChTabulatedActuatorModel({0}, {10}, {0}, {0}, {0.8}, {0.5}, {0}, {2}, smoothing);
    };
    auto model = make_model(0.2);
    // The regularized Coulomb law is tau_f=-F_c tanh(omega/s), with slope -F_c/s at zero.
    // Its power is always nonpositive, including when input is disabled.
    EXPECT_DOUBLE_EQ(model.GetFrictionSmoothingSpeed(), 0.2);
    EXPECT_DOUBLE_EQ(model.ComputeEffort({0, true}, {}), 0);
    for (double speed : {-10.0, -0.2, -1e-8, 0.0, 1e-8, 0.2, 10.0}) {
        double expected = -2 * std::tanh(speed / 0.2);
        auto telemetry = model.ComputeTelemetry({4, false}, {0, speed, 0});
        EXPECT_NEAR(telemetry.applied_effort, expected, 1e-12);
        EXPECT_LE(telemetry.mechanical_power, 0);
        EXPECT_FALSE(telemetry.command_saturated);
    }
    constexpr double delta = 1e-7;
    double left_slope = (model.ComputeEffort({0, true}, {}) -
                         model.ComputeEffort({0, true}, {0, -delta, 0})) /
                        delta;
    double right_slope = (model.ComputeEffort({0, true}, {0, delta, 0}) -
                          model.ComputeEffort({0, true}, {})) /
                         delta;
    EXPECT_NEAR(left_slope, -10, 1e-10);
    EXPECT_NEAR(right_slope, -10, 1e-10);
    EXPECT_NEAR(left_slope, right_slope, 1e-12);
    model.SetFrictionSmoothingSpeed(0.4);
    EXPECT_DOUBLE_EQ(model.GetFrictionSmoothingSpeed(), 0.4);
    EXPECT_NEAR(model.ComputeEffort({0, true}, {0, 0.2, 0}), -2 * std::tanh(0.5), 1e-12);
    for (double invalid : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::quiet_NaN()}) {
        EXPECT_THROW(make_model(invalid), std::invalid_argument);
        EXPECT_THROW(model.SetFrictionSmoothingSpeed(invalid), std::invalid_argument);
        EXPECT_DOUBLE_EQ(model.GetFrictionSmoothingSpeed(), 0.4);
    }
}

TEST(ActuatorModels, TabulatedTelemetryFlagsAndBidirectionalPower) {
    ChTabulatedActuatorModel model({0}, {10}, {0}, {0}, {0.8}, {0.5}, {0}, {2}, 0.2);
    // A table model clips only enabled normalized input, not a voltage/current/effort constraint.
    // Its terminal power is Pm/eta_m in motoring and Pm*eta_g in generating.
    for (double command : {-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0}) {
        for (double speed : {-2.0, 0.0, 2.0}) {
            for (bool enabled : {false, true}) {
                ActuatorCommand input{command, enabled};
                ActuatorState state{0, speed, 0};
                auto telemetry = model.ComputeTelemetry(input, state);
                double active_effort = enabled ? 10 * std::clamp(command, -1.0, 1.0) : 0;
                double effort = active_effort - 2 * std::tanh(speed / 0.2);
                double power = effort * speed;
                double expected_electrical = power >= 0 ? power / 0.8 : power * 0.5;
                EXPECT_EQ(telemetry.command_saturated, enabled && std::abs(command) > 1);
                EXPECT_FALSE(telemetry.current_saturated);
                EXPECT_FALSE(telemetry.voltage_saturated);
                EXPECT_FALSE(telemetry.effort_saturated);
                EXPECT_DOUBLE_EQ(telemetry.current, 0);
                EXPECT_DOUBLE_EQ(telemetry.copper_loss, 0);
                EXPECT_NEAR(telemetry.applied_effort, effort, 1e-12);
                EXPECT_NEAR(telemetry.mechanical_power, power, 1e-12);
                EXPECT_NEAR(telemetry.electrical_power, expected_electrical, 1e-12);
                EXPECT_GE(telemetry.electrical_power - telemetry.mechanical_power, 0);
                EXPECT_DOUBLE_EQ(telemetry.efficiency, power == 0 ? 0 : (power > 0 ? 0.8 : 0.5));
                EXPECT_DOUBLE_EQ(model.ComputeEffort(input, state), model.ComputeEffort(input, state));
            }
        }
    }
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
