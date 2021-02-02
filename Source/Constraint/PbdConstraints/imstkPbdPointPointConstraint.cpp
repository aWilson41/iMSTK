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

#pragma once

#include "imstkPbdPointPointConstraint.h"

namespace imstk
{
void
PbdPointPointConstraint::initConstraint(const Side side,
    Vec3d* pt1, double* invMassA1,
    Vec3d* pt2, double* invMassB1,
    std::shared_ptr<PbdCollisionConstraintConfig> configA, std::shared_ptr<PbdCollisionConstraintConfig> configB)
{
    m_side = side;
    m_normal = *pt2 - *pt1;
    m_penetrationDepth = m_normal.norm();
    m_normal = m_normal.normalized();
    m_bodiesFirst[0].vertex = pt1;
    m_bodiesFirst[0].invMass = invMassA1;
    m_bodiesSecond[0].vertex = pt2;
    m_bodiesSecond[0].invMass = invMassB1;
    m_configA = configA;
    m_configB = configB;
}

bool
PbdPointPointConstraint::computeValueAndGradient(double& c,
                                                 VecDataArray<double, 3>& dcdxA,
                                                 VecDataArray<double, 3>& dcdxB) const
{
    // Current position during solve
    const Vec3d& x0 = *m_bodiesFirst[0].vertex;
    const Vec3d& x1 = *m_bodiesSecond[0].vertex;
    const Vec3d diff = x1 - x0;

    if (m_side == Side::AB)
    {
        // \todo: Not implemented yet
    }
    // Only handle the movement of body A
    else if (m_side == Side::A)
    {
        // Actual penetration depth (thus far in solve)
        c = std::max(std::min(-diff.dot(-m_normal), m_penetrationDepth), 0.0);

        dcdxA[0] = -m_normal;
        dcdxB[0] = Vec3d::Zero();
    }
    else if (m_side == Side::B)
    {
        // Actual penetration depth (thus far in solve)
        c = std::max(std::min(diff.dot(m_normal), m_penetrationDepth), 0.0);

        dcdxA[0] = Vec3d::Zero();
        dcdxB[0] = m_normal;
    }

    return true;
}
} // imstk