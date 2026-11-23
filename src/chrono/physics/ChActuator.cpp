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
// Authors: Radu Serban
// =============================================================================

#include "chrono/utils/ChUtils.h"
#include "chrono/physics/ChActuator.h"

namespace chrono {

ChActuator::ChActuator() : is_attached(false), calculate_consistent_IC(false) {}

double ChActuator::GetInput(double t) const {
    return ref_fun ? ChClamp(ref_fun->GetVal(t), -1.0, +1.0) : 0;
}

void ChActuator::IntLoadResidual_F(const unsigned int off, ChVectorDynamic<>& R, const double c) {
    if (!IsActive())
        return;

    // Load dynamics of the actuator (if any)
    ChExternalDynamicsODE::IntLoadResidual_F(off, R, c);

    if (is_attached) {
        // Add forces to connected bodies (calculated in Update)
        if (m_body1->Variables().IsActive()) {
            R.segment(m_body1->Variables().GetOffset() + 0, 3) += c * m_Qforce.segment(0, 3);
            R.segment(m_body1->Variables().GetOffset() + 3, 3) += c * m_Qforce.segment(3, 3);
        }
        if (m_body2->Variables().IsActive()) {
            R.segment(m_body2->Variables().GetOffset() + 0, 3) += c * m_Qforce.segment(6, 3);
            R.segment(m_body2->Variables().GetOffset() + 3, 3) += c * m_Qforce.segment(9, 3);
        }
    }
}

void ChActuator::InjectKRMMatrices(ChSystemDescriptor& descriptor) {
    // Insert the optional coupling KRM block (if used)
    if (m_KRM_coupling.GetNumVariables() > 0)
        descriptor.InsertKRMBlock(&m_KRM_coupling);

    // Let the base class insert its own KRM block for actuator dynamics (if any)
    ChExternalDynamicsODE::InjectKRMMatrices(descriptor);
}

void ChActuator::LoadKRMMatrices(double Kfactor, double Rfactor, double Mfactor) {
    // Let derived classes load the coupling KRM block (if used)
    if (m_KRM_coupling.GetNumVariables() > 0)
        ComputeCouplingKRM(m_KRM_coupling.GetMatrix(), Kfactor, Rfactor, Mfactor);

    // Let the base class load its own KRM block for actuator dynamics (if any)
    ChExternalDynamicsODE::LoadKRMMatrices(Kfactor, Rfactor, Mfactor);
}

void ChActuator::SetupCouplingKRM() {
    // Drop any block from a previous Connect() (stale bodies or reallocated state variables).
    m_KRM_coupling = ChKRMBlock();
    if (!EnableCouplingKRM())
        return;
    std::vector<ChVariables*> vars;
    vars.push_back(&m_body1->Variables());
    vars.push_back(&m_body2->Variables());
    if (GetNumStates() > 0)
        vars.push_back(&Variables());
    m_KRM_coupling.SetVariables(vars);
}

// -----------------------------------------------------------------------------

ChLinearActuator::ChLinearActuator() : ChActuator(), s_0(0), s(0), sd(0), F0(0) {}

void ChLinearActuator::SetActuatorInitialLength(double len) {
    s_0 = len;
    s = len;
    sd = 0;
}

void ChLinearActuator::SetActuatorLength(double len, double vel) {
    // Do nothing if not attached to bodies
    if (is_attached)
        return;

    s = len;
    sd = vel;
}

void ChLinearActuator::SetInitialLoad(double initial_load) {
    calculate_consistent_IC = true;
    F0 = initial_load;
}

void ChLinearActuator::Connect(std::shared_ptr<ChBody> body1, std::shared_ptr<ChBody> body2, bool local, ChVector3d loc1, ChVector3d loc2) {
    is_attached = true;

    // Call base initialization (external dynamics)
    Initialize();

    // Cache connected bodies and body local connection points
    m_body1 = body1.get();
    m_body2 = body2.get();

    if (local) {
        m_loc1 = loc1;
        m_loc2 = loc2;
        m_aloc1 = body1->TransformPointLocalToParent(loc1);
        m_aloc2 = body2->TransformPointLocalToParent(loc2);
    } else {
        m_loc1 = body1->TransformPointParentToLocal(loc1);
        m_loc2 = body2->TransformPointParentToLocal(loc2);
        m_aloc1 = loc1;
        m_aloc2 = loc2;
    }

    s_0 = (m_aloc1 - m_aloc2).Length();
    s = s_0;
    sd = 0;

    // Resize temporary vectors of generalized body forces
    m_Qforce.setZero(12);
    m_Qunit.setZero(12);

    // If enabled, set up the coupling KRM block
    SetupCouplingKRM();
}

void ChLinearActuator::Update(double time, UpdateFlags update_flags) {
    // If the actuator is attached to bodies
    // - update its length and rate from the body states
    // - apply the generated force to the two bodies
    if (is_attached) {
        m_aloc1 = m_body1->TransformPointLocalToParent(m_loc1);
        m_aloc2 = m_body2->TransformPointLocalToParent(m_loc2);

        auto avel1 = m_body1->PointSpeedLocalToParent(m_loc1);
        auto avel2 = m_body2->PointSpeedLocalToParent(m_loc2);

        ChVector3d dir = (m_aloc1 - m_aloc2).GetNormalized();

        s = (m_aloc1 - m_aloc2).Length();
        sd = Vdot(dir, avel1 - avel2);

        // Generalized force per unit actuator force (also d(sd)/d(body velocities))
        auto atorque1 = Vcross(m_aloc1 - m_body1->GetPos(), dir);   // absolute frame
        auto atorque2 = Vcross(m_aloc2 - m_body2->GetPos(), -dir);  // absolute frame
        m_Qunit.segment(0, 3) = dir.eigen();
        m_Qunit.segment(3, 3) = m_body1->TransformDirectionParentToLocal(atorque1).eigen();
        m_Qunit.segment(6, 3) = (-dir).eigen();
        m_Qunit.segment(9, 3) = m_body2->TransformDirectionParentToLocal(atorque2).eigen();
    }

    // Internal dynamics must use the current connected-body kinematics.
    ChExternalDynamicsODE::Update(time, update_flags);

    if (is_attached) {
        // Actuator force, applied as equal and opposite forces on the two bodies
        m_Qforce = GetActuatorForce(time) * m_Qunit;
    }
}

// -----------------------------------------------------------------------------

ChRotationalActuator::ChRotationalActuator() : ChActuator(), a_0(0), a(0), ad(0), T0(0), m_turns(0), m_last_angle(0) {}

void ChRotationalActuator::SetActuatorInitialAngle(double angle) {
    a_0 = angle;
    a = angle;
    ad = 0;
}

void ChRotationalActuator::SetActuatorAngle(double angle, double ang_vel) {
    // Do nothing if not attached to bodies
    if (is_attached)
        return;

    a = angle;
    ad = ang_vel;
}

void ChRotationalActuator::SetInitialLoad(double initial_load) {
    calculate_consistent_IC = true;
    T0 = initial_load;
}

void ChRotationalActuator::CalcAngle() {
    // Express the RSDA frames in absolute frame
    ChQuaternion<> rot1 = m_body1->GetRot() * m_csys1.rot;
    ChQuaternion<> rot2 = m_body2->GetRot() * m_csys2.rot;

    // Extract unit vectors
    auto f1 = rot1.GetAxisX();
    auto g1 = rot1.GetAxisY();
    auto f2 = rot2.GetAxisX();

    // Calculate sine and cosine of rotation angle
    double s = Vdot(g1, f2);
    double c = Vdot(f1, f2);

    // Get angle (in [-pi , +pi])
    a = std::atan2(s, c);

    // Get angle rate of change
    m_axis = rot1.GetAxisZ();
    ad = Vdot(m_axis, m_body2->GetAngVelParent() - m_body1->GetAngVelParent());
}

void ChRotationalActuator::AdjustAngle() {
    // Check cross at +- pi
    if (m_last_angle > CH_PI_2 && a < 0)
        a += CH_2PI;
    if (m_last_angle < -CH_PI_2 && a > 0)
        a -= CH_2PI;

    // Accumulate full turns
    if (m_last_angle - a > CH_PI)
        m_turns++;
    if (m_last_angle - a < -CH_PI)
        m_turns--;

    if (a < 0) {
        while (m_turns > 0) {
            a += CH_2PI;
            m_turns--;
        }
    }
    if (a > 0) {
        while (m_turns < 0) {
            a -= CH_2PI;
            m_turns++;
        }
    }

    // Update last angle
    m_last_angle = a;
}

void ChRotationalActuator::Connect(std::shared_ptr<ChBody> body1, std::shared_ptr<ChBody> body2, const ChFrame<>& frame) {
    is_attached = true;

    // Call base initialization (external dynamics)
    Initialize();

    // Cache connected bodies
    m_body1 = body1.get();
    m_body2 = body2.get();

    // Coordinate frames on the two bodies
    m_csys1 = m_body1->GetCoordsys().TransformParentToLocal(frame.GetCoordsys());
    m_csys2 = m_body2->GetCoordsys().TransformParentToLocal(frame.GetCoordsys());

    // Calculate initial angle
    CalcAngle();
    a_0 = a;
    m_turns = 0;
    m_last_angle = a;

    // Resize temporary vectors of generalized body forces
    m_Qforce.setZero(12);
    m_Qunit.setZero(12);

    // If enabled, set up the coupling KRM block
    SetupCouplingKRM();
}

void ChRotationalActuator::Update(double time, UpdateFlags update_flags) {
    // If the actuator is attached to bodies
    // - update its angle and rate from the body states
    // - apply the generated torque to the two bodies
    if (is_attached) {
        // Calculate current angle and angle rate
        CalcAngle();
        AdjustAngle();

        // Generalized force per unit actuator torque (also d(ad)/d(body velocities)).
        // The angle measures body 2 relative to body 1, so body 2 receives +axis and body 1 the reaction -axis.
        m_Qunit.segment(0, 3).setZero();
        m_Qunit.segment(3, 3) = m_body1->TransformDirectionParentToLocal(-m_axis).eigen();
        m_Qunit.segment(6, 3).setZero();
        m_Qunit.segment(9, 3) = m_body2->TransformDirectionParentToLocal(m_axis).eigen();
    }

    // Internal dynamics must use the current connected-body kinematics.
    ChExternalDynamicsODE::Update(time, update_flags);

    if (is_attached) {
        // Actuator torque, applied as equal and opposite torques on the two bodies
        m_Qforce = GetActuatorTorque(time) * m_Qunit;
    }
}

}  // end namespace chrono
