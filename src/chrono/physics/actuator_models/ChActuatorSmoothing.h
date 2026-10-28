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

#ifndef CH_ACTUATOR_SMOOTHING_H
#define CH_ACTUATOR_SMOOTHING_H

#include <algorithm>
#include <cmath>

namespace chrono {
namespace actuator_models {

/// Compact C1 smoothstep: h(q)=q^2(3-2q) on [0,1], constant outside; max h'=3/2.
inline double ActuatorSmoothStep(double q) {
    q = std::clamp(q, 0.0, 1.0);
    return q * q * (3 - 2 * q);
}

/// Passive load-dependent loss torque (or force) at the output.
/// With x=velocity/smoothing_speed, s(x)=sign(x)*h(2|x|) reaches +/-1 at |x|=1/2.
/// The directional coefficient is cm*h(|x|) or cg*h(|x|), smoothly continued through zero at standstill.
/// This prevents the opposite quadrant's losses from reversing motoring effort inside the band.
/// |d(loss)/d(velocity)| <= 3*max(cm,cg)*|effort|/smoothing_speed (see manual for the product bound).
/// Both factors are C1 in velocity, loss*velocity>=0, and the directional branches are exact at |x|>=1.
/// The caller must supply finite smoothing_speed>0 and nonnegative loss coefficients cm and cg.
inline double ActuatorTransmissionLoss(double effort, double velocity, double smoothing_speed, double cm, double cg) {
    double x = velocity / smoothing_speed;
    double magnitude = std::abs(x);
    double sign = std::copysign(ActuatorSmoothStep(2 * magnitude), x);
    double directional_coefficient = std::signbit(effort) == std::signbit(velocity) ? cm : cg;
    double coefficient = directional_coefficient * ActuatorSmoothStep(magnitude);
    return coefficient * std::abs(effort) * sign;
}

}  // namespace actuator_models
}  // namespace chrono

#endif
