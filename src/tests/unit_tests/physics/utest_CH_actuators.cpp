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

#include "gtest/gtest.h"

#include "chrono/functions/ChFunctionConst.h"
#include "chrono/physics/ChActuatorRSDA.h"
#include "chrono/physics/ChActuatorTSDA.h"

using namespace chrono;

TEST(ActuatorFixes, TSDAConstantForce) {
    ChActuatorTSDA actuator;
    actuator.SetInputFunction(chrono_types::make_shared<ChFunctionConst>(0.5));
    actuator.SetConstantForce(3);
    actuator.SetBaseForce(2);

    EXPECT_DOUBLE_EQ(actuator.GetActuatorForce(0), 4);
}

TEST(ActuatorFixes, RSDATorqueAndRestAngle) {
    ChActuatorRSDA actuator;
    actuator.SetInputFunction(chrono_types::make_shared<ChFunctionConst>(0.5));
    actuator.SetActuatorAngle(2, 0);
    actuator.SetRestAngle(1);
    actuator.SetConstantTorque(4);
    actuator.SetBaseTorque(2);
    actuator.SetSpringCoefficient(3);

    EXPECT_DOUBLE_EQ(actuator.GetActuatorTorque(0), 2);
}
