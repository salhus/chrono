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
// A DC-equivalent actuator drives a frictionless rotor: regeneration from above no-load
// speed, motoring from rest, then sinusoidal voltage with repeated speed reversals. Kt = Ke in SI.
// The quasi-static winding (L = 0) gives omega_0 = V/(N*Ke) = 30 rad/s.
// Generating at 40 rad/s: i=-4 A, torque=-4 N*m, Pbus=-96 W, Pcopper=32 W,
// Pmachine=-128 W, Pgear=32 W, Pmech=-160 W.
// Motoring at rest: i=12 A, torque=9.6 N*m, Pbus=Pcopper=288 W, other powers zero.
// Just outside the 0.01 rad/s smoothing band, motoring torque is approximately 8.64 N*m.
// Outside that speed band, omega - omega_0 decays exponentially:
// tau_m = J*R/(N^2*Kt*Ke*eta_m) = 0.347222 s; tau_g = J*R*eta_g/(N^2*Kt*Ke) = 0.25 s.
// Bus power is negative during regeneration. The energy residual should tend to zero
// as the integration step decreases; copper and gearbox losses remain nonnegative.
// =============================================================================

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>

#include "chrono/functions/ChFunctionConst.h"
#include "chrono/functions/ChFunctionSine.h"
#include "chrono/physics/ChActuatorModelRotational.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/actuator_models/ChElectricActuatorModel.h"

using namespace chrono;
using namespace chrono::actuator_models;

void Run(double initial_speed, const char* phase, bool oscillatory = false) {
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
    if (oscillatory)
        actuator->SetInputFunction(chrono_types::make_shared<ChFunctionSine>(0.1, 1));
    else
        actuator->SetInputFunction(chrono_types::make_shared<ChFunctionConst>(1));
    actuator->Connect(ground, rotor, ChFrame<>(VNULL, QUNIT));
    system.Add(actuator);

    std::cout << "\n" << phase << ": no-load speed = " << model->GetNoLoadSpeed() << " rad/s\n";
    std::cout << "time[s],speed[rad/s],torque[N*m],mechanical_power[W],bus_power[W],copper_loss[W],gear_loss[W],current[A],efficiency[1]\n";
    const double step = 1e-4;
    const int steps = oscillatory ? 40000 : 20000;
    double bus_energy = 0;
    double copper_energy = 0;
    double gear_energy = 0;
    auto previous = actuator->GetTelemetry(0);
    for (int i = 0; i <= steps; ++i) {
        double time = system.GetChTime();
        double speed = rotor->GetAngVelParent().z();
        auto telemetry = actuator->GetTelemetry(time);
        if (i > 0) {
            bus_energy += 0.5 * step * (previous.electrical_power + telemetry.electrical_power);
            copper_energy += 0.5 * step * (previous.copper_loss + telemetry.copper_loss);
            gear_energy += 0.5 * step * (previous.gear_loss + telemetry.gear_loss);
        }
        if (i % 2000 == 0 || (oscillatory && std::abs(speed) < 3 * parameters.loss_smoothing_speed && i % 10 == 0)) {
            std::cout << time << "," << speed << "," << telemetry.applied_effort << "," << telemetry.mechanical_power << "," << telemetry.electrical_power << ","
                      << telemetry.copper_loss << "," << telemetry.gear_loss << "," << telemetry.current << "," << telemetry.efficiency << "\n";
        }
        previous = telemetry;
        if (i < steps)
            system.DoStepDynamics(step);
    }

    double speed = rotor->GetAngVelParent().z();
    double kinetic_change = 0.5 * inertia * (speed * speed - initial_speed * initial_speed);
    double residual = bus_energy - kinetic_change - copper_energy - gear_energy;
    double scale = std::max({std::abs(bus_energy), std::abs(kinetic_change), copper_energy, gear_energy});
    std::cout << "Energy balance [J]: bus_drawn=" << bus_energy << ", delta_kinetic=" << kinetic_change << ", copper_loss=" << copper_energy << ", gear_loss=" << gear_energy
              << ", residual=" << residual << ", relative_residual=" << std::abs(residual) / scale << "\n";
}

int main() {
    std::cout << std::setprecision(9);
    Run(40, "Generating");
    Run(0, "Motoring");
    Run(0, "Reversals (0.1*sin(2*pi*t) voltage command)", true);
    return 0;
}
