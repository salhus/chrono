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

#ifndef CH_TABULATED_ACTUATOR_MODEL_H
#define CH_TABULATED_ACTUATOR_MODEL_H

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "ChActuatorModel.h"

namespace chrono {
namespace actuator_models {

/// Stateless actuator model consuming user-supplied output-side effort-speed, efficiency, and friction tables.
/// Axes and effort magnitudes are nonnegative SI quantities; commands specify a fraction of the effort envelope.
/// Tables use linear/bilinear interpolation and clamp outside their axes. Data must describe a passive actuator.
class ChTabulatedActuatorModel : public ChActuatorModel {
  public:
    /// Efficiency tables are row-major: speed index first, then effort index.
    /// Each axis is nonempty, finite, and strictly increasing. Efficiencies lie in (0, 1].
    /// The optional positive, finite smoothing speed has units m/s or rad/s; existing constructors remain valid.
    ChTabulatedActuatorModel(std::vector<double> speed_axis,
                             std::vector<double> max_effort,
                             std::vector<double> efficiency_speed_axis,
                             std::vector<double> efficiency_effort_axis,
                             std::vector<double> motoring_efficiency,
                             std::vector<double> generating_efficiency,
                             std::vector<double> friction_speed_axis,
                             std::vector<double> friction_effort,
                             double friction_smoothing_speed = 1e-3)
        : m_speed_axis(std::move(speed_axis)),
          m_max_effort(std::move(max_effort)),
          m_efficiency_speed_axis(std::move(efficiency_speed_axis)),
          m_efficiency_effort_axis(std::move(efficiency_effort_axis)),
          m_motoring_efficiency(std::move(motoring_efficiency)),
          m_generating_efficiency(std::move(generating_efficiency)),
          m_friction_speed_axis(std::move(friction_speed_axis)),
          m_friction_effort(std::move(friction_effort)),
          m_friction_smoothing_speed(friction_smoothing_speed) {
        SetFrictionSmoothingSpeed(friction_smoothing_speed);
        ValidateAxis(m_speed_axis, m_max_effort.size(), "speed");
        ValidateAxis(m_efficiency_speed_axis, m_efficiency_speed_axis.size(), "efficiency speed");
        ValidateAxis(m_efficiency_effort_axis, m_efficiency_effort_axis.size(), "efficiency effort");
        ValidateAxis(m_friction_speed_axis, m_friction_effort.size(), "friction speed");
        ValidateNonnegativeAxis(m_speed_axis, "speed");
        ValidateNonnegativeAxis(m_efficiency_speed_axis, "efficiency speed");
        ValidateNonnegativeAxis(m_efficiency_effort_axis, "efficiency effort");
        ValidateNonnegativeAxis(m_friction_speed_axis, "friction speed");
        const auto map_size = m_efficiency_speed_axis.size() * m_efficiency_effort_axis.size();
        if (m_motoring_efficiency.size() != map_size || m_generating_efficiency.size() != map_size)
            throw std::invalid_argument("Efficiency map dimensions do not match their axes");
        ValidateEfficiency(m_motoring_efficiency);
        ValidateEfficiency(m_generating_efficiency);
        ValidateValues(m_max_effort);
        ValidateValues(m_friction_effort);
        for (double value : m_max_effort) {
            if (value < 0)
                throw std::invalid_argument("Maximum effort values must be nonnegative");
        }
        for (double value : m_friction_effort) {
            if (value < 0)
                throw std::invalid_argument("Friction values must be nonnegative");
        }
    }

    /// Set the velocity scale for the smooth Coulomb sign; changing it does not change table magnitudes.
    void SetFrictionSmoothingSpeed(double speed) {
        if (!std::isfinite(speed) || speed <= 0)
            throw std::invalid_argument("Friction smoothing speed must be positive and finite");
        m_friction_smoothing_speed = speed;
    }

    /// Get the friction smoothing speed [m/s or rad/s].
    double GetFrictionSmoothingSpeed() const { return m_friction_smoothing_speed; }

    /// Apply the signed command fraction of the speed envelope, then subtract smooth output-side friction.
    double ComputeEffort(const ActuatorCommand& command, const ActuatorState& state) const override {
        double speed = std::abs(state.velocity);
        double envelope = Interpolate1D(m_speed_axis, m_max_effort, speed);
        double commanded_effort = command.enabled ? std::clamp(command.value, -1.0, 1.0) * envelope : 0;
        double friction = Interpolate1D(m_friction_speed_axis, m_friction_effort, speed);
        double smooth_sign = std::tanh(state.velocity / m_friction_smoothing_speed);
        return commanded_effort - friction * smooth_sign;
    }

    /// Return bilinear efficiency using the sign of final mechanical power to select the quadrant.
    double GetEfficiency(const ActuatorCommand& command, const ActuatorState& state) const {
        double effort = ComputeEffort(command, state);
        const auto& map = effort * state.velocity >= 0 ? m_motoring_efficiency : m_generating_efficiency;
        return Interpolate2D(m_efficiency_speed_axis, m_efficiency_effort_axis, map, std::abs(state.velocity), std::abs(effort));
    }

    /// Estimate bus power as P_mech/eta for motoring or eta*P_mech for generating.
    /// Maps must include all device losses; no winding-current or separate copper-loss estimate is available.
    ActuatorTelemetry ComputeTelemetry(const ActuatorCommand& command, const ActuatorState& state) const override {
        ActuatorTelemetry telemetry;
        telemetry.applied_effort = ComputeEffort(command, state);
        telemetry.mechanical_power = telemetry.applied_effort * state.velocity;
        telemetry.efficiency = GetEfficiency(command, state);
        telemetry.electrical_power = telemetry.mechanical_power >= 0 ? telemetry.mechanical_power / telemetry.efficiency : telemetry.mechanical_power * telemetry.efficiency;
        if (telemetry.mechanical_power == 0)
            telemetry.efficiency = 0;
        telemetry.command_saturated = command.enabled && std::abs(command.value) > 1;
        return telemetry;
    }

  private:
    static void ValidateValues(const std::vector<double>& values) {
        for (double value : values) {
            if (!std::isfinite(value))
                throw std::invalid_argument("Table values must be finite");
        }
    }

    static void ValidateAxis(const std::vector<double>& axis, std::size_t values_size, const char* name) {
        if (axis.empty() || axis.size() != values_size)
            throw std::invalid_argument(std::string("Invalid ") + name + " axis size");
        ValidateValues(axis);
        for (std::size_t i = 1; i < axis.size(); ++i) {
            if (!(axis[i] > axis[i - 1]))
                throw std::invalid_argument(std::string(name) + " axis must be strictly increasing");
        }
    }

    static void ValidateNonnegativeAxis(const std::vector<double>& axis, const char* name) {
        if (axis.front() < 0)
            throw std::invalid_argument(std::string(name) + " axis must be nonnegative");
    }

    static void ValidateEfficiency(const std::vector<double>& values) {
        ValidateValues(values);
        for (double value : values) {
            if (!(value > 0 && value <= 1))
                throw std::invalid_argument("Efficiency values must be in (0, 1]");
        }
    }

    static std::size_t LowerIndex(const std::vector<double>& axis, double value) {
        if (axis.size() < 2 || value <= axis.front())
            return 0;
        if (value >= axis.back())
            return axis.size() - 2;
        return static_cast<std::size_t>(std::upper_bound(axis.begin(), axis.end(), value) - axis.begin() - 1);
    }

    static double Interpolate1D(const std::vector<double>& axis, const std::vector<double>& values, double x) {
        if (axis.size() == 1 || x <= axis.front())
            return values.front();
        if (x >= axis.back())
            return values.back();
        std::size_t i = LowerIndex(axis, x);
        double weight = (x - axis[i]) / (axis[i + 1] - axis[i]);
        return values[i] + weight * (values[i + 1] - values[i]);
    }

    static double Interpolate2D(const std::vector<double>& x_axis, const std::vector<double>& y_axis, const std::vector<double>& values, double x, double y) {
        std::size_t ix = LowerIndex(x_axis, x);
        std::size_t iy = LowerIndex(y_axis, y);
        std::size_t ix1 = std::min(ix + 1, x_axis.size() - 1);
        std::size_t iy1 = std::min(iy + 1, y_axis.size() - 1);
        double wx = ix1 == ix ? 0 : std::clamp((x - x_axis[ix]) / (x_axis[ix1] - x_axis[ix]), 0.0, 1.0);
        double wy = iy1 == iy ? 0 : std::clamp((y - y_axis[iy]) / (y_axis[iy1] - y_axis[iy]), 0.0, 1.0);
        auto at = [&](std::size_t i, std::size_t j) { return values[i * y_axis.size() + j]; };
        double low = at(ix, iy) + wy * (at(ix, iy1) - at(ix, iy));
        double high = at(ix1, iy) + wy * (at(ix1, iy1) - at(ix1, iy));
        return low + wx * (high - low);
    }

    std::vector<double> m_speed_axis;
    std::vector<double> m_max_effort;
    std::vector<double> m_efficiency_speed_axis;
    std::vector<double> m_efficiency_effort_axis;
    std::vector<double> m_motoring_efficiency;
    std::vector<double> m_generating_efficiency;
    std::vector<double> m_friction_speed_axis;
    std::vector<double> m_friction_effort;
    double m_friction_smoothing_speed;
};

}  // namespace actuator_models
}  // namespace chrono

#endif
