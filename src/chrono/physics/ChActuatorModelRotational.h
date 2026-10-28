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

#ifndef CH_ACTUATOR_MODEL_ROTATIONAL_H
#define CH_ACTUATOR_MODEL_ROTATIONAL_H

#include <memory>

#include "chrono/physics/ChActuator.h"
#include "chrono/physics/actuator_models/ChActuatorModel.h"

namespace chrono {

/// Chrono binding for an engine-independent rotational actuator model.
class ChApi ChActuatorModelRotational : public ChRotationalActuator {
  public:
    explicit ChActuatorModelRotational(std::shared_ptr<const actuator_models::ChActuatorModel> model);
    virtual ~ChActuatorModelRotational() {}

    /// Get the actuator model.
    std::shared_ptr<const actuator_models::ChActuatorModel> GetModel() const { return m_model; }

    /// Get model effort using the current actuator state and input.
    virtual double GetActuatorTorque(double time) override;

    /// Forward the model's solver stiffness declaration.
    virtual bool IsStiff() const override;

  protected:
    virtual unsigned int GetNumStates() const override;
    virtual void SetInitialConditions(ChVectorDynamic<>& y0) override;
    virtual void CalculateRHS(double time, const ChVectorDynamic<>& y, ChVectorDynamic<>& rhs) override;
    virtual bool CalculateJac(double time, const ChVectorDynamic<>& y, const ChVectorDynamic<>& rhs, ChMatrixDynamic<>& J) override;

  private:
    std::shared_ptr<const actuator_models::ChActuatorModel> m_model;
};

}  // namespace chrono

#endif
