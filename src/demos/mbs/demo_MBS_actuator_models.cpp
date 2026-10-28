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

#include <iostream>
#include <memory>

#include "chrono/functions/ChFunctionConst.h"
#include "chrono/physics/ChActuatorModelRotational.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/actuator_models/ChElectricActuatorModel.h"

using namespace chrono;
using namespace chrono::actuator_models;

int main() {
    ChSystemNSC system;
    system.SetGravitationalAcceleration(VNULL);
    system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT);
    system.SetSolverType(ChSolver::Type::SPARSE_QR);

    auto ground = chrono_types::make_shared<ChBody>();
    ground->SetFixed(true);
    system.AddBody(ground);

    auto rotor = chrono_types::make_shared<ChBody>();
    rotor->SetMass(1);
    rotor->SetInertiaXX(ChVector3d(0.1, 0.1, 0.1));
    rotor->SetAngVelParent(ChVector3d(0, 0, 40));
    system.AddBody(rotor);

    ChElectricActuatorModel::Parameters parameters;
    parameters.resistance = 2;
    parameters.inductance = 0.02;
    parameters.torque_constant = 0.1;
    parameters.back_emf_constant = 0.2;
    parameters.bus_voltage = 24;
    parameters.gear_ratio = 4;
    parameters.motoring_efficiency = 0.9;
    parameters.generating_efficiency = 0.8;
    parameters.max_current = 20;
    auto model = std::make_shared<ChElectricActuatorModel>(parameters);

    auto actuator = chrono_types::make_shared<ChActuatorModelRotational>(model);
    actuator->SetInputFunction(chrono_types::make_shared<ChFunctionConst>(1));
    actuator->Connect(ground, rotor, ChFrame<>(VNULL, QUNIT));
    system.Add(actuator);

    std::cout << "time, speed, torque, mechanical power, winding current\n";
    const double step = 1e-3;
    for (int i = 0; i <= 2000; ++i) {
        double time = system.GetChTime();
        if (i % 100 == 0) {
            double torque = actuator->GetActuatorTorque(time);
            std::cout << time << ", " << rotor->GetAngVelParent().z() << ", " << torque << ", " << torque * rotor->GetAngVelParent().z() << ", " << actuator->GetStates()(0)
                      << "\n";
        }
        system.DoStepDynamics(step);
    }

    return 0;
}
