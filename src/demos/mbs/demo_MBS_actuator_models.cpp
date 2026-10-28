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
//
// A DC-equivalent actuator drives a frictionless rotor in two independent runs:
// regeneration from above no-load speed and motoring from rest. Kt = Ke in SI.
// The quasi-static winding (L = 0) gives omega_0 = V/(N*Ke) = 30 rad/s.
// Outside the small power-smoothing band, omega - omega_0 decays exponentially:
// tau_m = J*R/(N^2*Kt*Ke*eta_m) = 0.347222 s; tau_g = J*R*eta_g/(N^2*Kt*Ke) = 0.25 s.
// Bus power is negative during regeneration. The energy residual should tend to zero
// as the integration step decreases; copper and gearbox losses remain nonnegative.
// =============================================================================

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <memory>

#include "chrono/functions/ChFunctionConst.h"
#include "chrono/physics/ChActuatorModelRotational.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/actuator_models/ChElectricActuatorModel.h"

using namespace chrono;
using namespace chrono::actuator_models;

void Run(double initial_speed, const char* phase) {
    ChSystemNSC system;
    system.SetGravitationalAcceleration(VNULL);
    system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT);
    system.SetSolverType(ChSolver::Type::SPARSE_QR);

    auto ground = chrono_types::make_shared<ChBody>();
    ground->SetFixed(true);
    system.AddBody(ground);

    auto rotor = chrono_types::make_shared<ChBody>();
    rotor->SetMass(1);
    const double inertia = 0.1;  // [kg*m^2]
    rotor->SetInertiaXX(ChVector3d(inertia, inertia, inertia));
    rotor->SetAngVelParent(ChVector3d(0, 0, initial_speed));
    system.AddBody(rotor);

    ChElectricActuatorModel::Parameters parameters;
    parameters.resistance = 2;           // [ohm]
    parameters.inductance = 0;           // [H], quasi-static electrical limit
    parameters.torque_constant = 0.2;    // [N*m/A]
    parameters.back_emf_constant = 0.2;  // [V*s/rad], equal to Kt in SI
    parameters.bus_voltage = 24;         // [V]
    parameters.gear_ratio = 4;
    parameters.motoring_efficiency = 0.9;
    parameters.generating_efficiency = 0.8;
    parameters.max_current = 20;  // [A], not reached in either run
    auto model = std::make_shared<ChElectricActuatorModel>(parameters);

    auto actuator = chrono_types::make_shared<ChActuatorModelRotational>(model);
    actuator->SetInputFunction(chrono_types::make_shared<ChFunctionConst>(1));
    actuator->Connect(ground, rotor, ChFrame<>(VNULL, QUNIT));
    system.Add(actuator);

    std::cout << "\n" << phase << ": no-load speed = " << model->GetNoLoadSpeed() << " rad/s\n";
    std::cout << "time[s],speed[rad/s],torque[N*m],mechanical_power[W],bus_power[W],copper_loss[W],current[A],efficiency[1]\n";
    const double step = 1e-4;
    const int steps = 20000;
    double bus_energy = 0;
    double copper_energy = 0;
    double gear_energy = 0;
    auto previous = actuator->GetTelemetry(0);
    auto gear_loss = [&](const ActuatorTelemetry& telemetry, double speed) {
        double electromagnetic_power = parameters.gear_ratio * parameters.torque_constant * telemetry.current * speed;
        return std::max(0.0, electromagnetic_power - telemetry.mechanical_power);
    };
    double previous_gear_loss = gear_loss(previous, initial_speed);
    for (int i = 0; i <= steps; ++i) {
        double time = system.GetChTime();
        double speed = rotor->GetAngVelParent().z();
        auto telemetry = actuator->GetTelemetry(time);
        double loss = gear_loss(telemetry, speed);
        if (i > 0) {
            bus_energy += 0.5 * step * (previous.electrical_power + telemetry.electrical_power);
            copper_energy += 0.5 * step * (previous.copper_loss + telemetry.copper_loss);
            gear_energy += 0.5 * step * (previous_gear_loss + loss);
        }
        if (i % 2000 == 0) {
            std::cout << time << "," << speed << "," << telemetry.applied_effort << "," << telemetry.mechanical_power << "," << telemetry.electrical_power << ","
                      << telemetry.copper_loss << "," << telemetry.current << "," << telemetry.efficiency << "\n";
        }
        previous = telemetry;
        previous_gear_loss = loss;
        if (i < steps)
            system.DoStepDynamics(step);
    }

    double speed = rotor->GetAngVelParent().z();
    double kinetic_change = 0.5 * inertia * (speed * speed - initial_speed * initial_speed);
    std::cout << "Energy balance [J]: bus_drawn=" << bus_energy << ", delta_kinetic=" << kinetic_change << ", copper_loss=" << copper_energy << ", gear_loss=" << gear_energy
              << ", residual=" << bus_energy - kinetic_change - copper_energy - gear_energy << "\n";
}

int main() {
    std::cout << std::setprecision(9);
    Run(40, "Generating");
    Run(0, "Motoring");
    return 0;
}
