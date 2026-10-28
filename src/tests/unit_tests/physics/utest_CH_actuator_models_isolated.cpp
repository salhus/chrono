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
    return electric_model.ComputeEffort(command, state) == 0 ? 0 : 1;
}
