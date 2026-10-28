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

#include "ChActuatorModel.h"
#include "ChElectricActuatorModel.h"
#include "ChTabulatedActuatorModel.h"

int main() {
    chrono::actuator_models::ChElectricActuatorModel electric_model;
    chrono::actuator_models::ActuatorCommand command;
    chrono::actuator_models::ActuatorState state;
    chrono::actuator_models::ChTabulatedActuatorModel tabulated_model({0, 10}, {10, 5}, {0, 10}, {0, 10}, {0.8, 0.8, 0.8, 0.8}, {0.6, 0.6, 0.6, 0.6}, {0, 10}, {0, 2});
    if (electric_model.ComputeEffort(command, state) != 0 || tabulated_model.ComputeEffort(command, state) != 0)
        return 1;
    command.value = 0.5;
    state.velocity = 5;
    auto electric = electric_model.ComputeTelemetry(command, state);
    auto tabulated = tabulated_model.ComputeTelemetry(command, state);
    // Both concrete models must instantiate, execute, and expose telemetry without Chrono or Eigen linkage.
    if (!(electric.current > 0 && electric.copper_loss > 0 && tabulated.applied_effort > 0 && tabulated.mechanical_power > 0 &&
          tabulated.electrical_power > tabulated.mechanical_power))
        return 2;
    command.enabled = false;
    return electric_model.ComputeEffort(command, state) < 0 && tabulated_model.ComputeEffort(command, state) < 0 ? 0 : 3;
}
