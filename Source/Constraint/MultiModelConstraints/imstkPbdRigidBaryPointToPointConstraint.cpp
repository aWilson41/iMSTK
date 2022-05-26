/*=========================================================================

   Library: iMSTK

   Copyright (c) Kitware, Inc. & Center for Modeling, Simulation,
   & Imaging in Medicine, Rensselaer Polytechnic Institute.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0.txt

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

=========================================================================*/

#include "imstkPbdRigidBaryPointToPointConstraint.h"

namespace imstk
{
PbdRigidBaryPointToPointConstraint::PbdRigidBaryPointToPointConstraint(std::shared_ptr<RigidBody> obj1) :
    PbdBaryPointToPointConstraint(),
    RbdConstraint(
        obj1,
        nullptr,
        RbdConstraint::Side::A)
{
}

///
/// \brief compute value and gradient of constraint function and weight it
/// by half to force the motion to the half way point between two bodies
///
/// \param[inout] c constraint value
/// \param[inout] dcdxA constraint gradient for A
/// \param[inout] dcdxB constraint gradient for B
bool
PbdRigidBaryPointToPointConstraint::computeValueAndGradient(PbdState& bodies,
                                                            double& c, std::vector<Vec3d>& dcdx) const
{
    // Compute the middle position between the point on the rigid body and the PBD object
    m_diff = 0.5 * computeInterpolantDifference(bodies);

    c = m_diff.norm();

    if (c < IMSTK_DOUBLE_EPS)
    {
        m_diff = Vec3d::Zero();
        return false;
    }
    m_diff /= c;

    for (size_t i = 0; i < dcdx.size(); i++)
    {
        if (m_bodiesSides[i])
        {
            dcdx[i] = -m_diff * m_weights[i];
        }
        else
        {
            dcdx[i] = m_diff * m_weights[i];
        }
    }

    return true;
}

void
PbdRigidBaryPointToPointConstraint::compute(double dt)
{
    J = Eigen::Matrix<double, 3, 4>::Zero();

    J(0, 0) = -m_diff[0]; J(0, 1) = 0.0;
    J(1, 0) = -m_diff[1]; J(1, 1) = 0.0;
    J(2, 0) = -m_diff[2]; J(2, 1) = 0.0;

    // B stabilization term
    vu = m_diff.norm() * m_beta / dt;
}
} // namespace imstk