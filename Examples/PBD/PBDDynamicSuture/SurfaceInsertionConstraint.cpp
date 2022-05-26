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

#include "SurfaceInsertionConstraint.h"

using namespace imstk;

void
SurfaceInsertionConstraint::initConstraint(
    const Vec3d&         insertionPoint,
    const PbdParticleId& ptB1,
    const PbdParticleId& ptB2,
    const PbdParticleId& ptB3,
    const Vec3d&         contactPt,
    const Vec3d&         barycentricPt,
    double               stiffnessA,
    double               stiffnessB)
{
    m_insertionPoint = insertionPoint;
    m_contactPt      = contactPt;
    m_barycentricPt  = barycentricPt;

    m_particles[0] = { -1, 0 }; // Doesn't matter
    m_particles[1] = ptB1;
    m_particles[2] = ptB2;
    m_particles[3] = ptB3;

    m_stiffness[0] = stiffnessA;
    m_stiffness[1] = stiffnessB;
}

bool
SurfaceInsertionConstraint::computeValueAndGradient(PbdState& bodies,
                                                    double& c, std::vector<Vec3d>& dcdx) const
{
    // Get current position of puncture point
    // Move triangle to match motion of needle

    Vec3d diff = m_contactPt - m_insertionPoint;

    c = diff.norm();

    // If sufficiently close, do not solve constraint
    if (c < 1E-8)
    {
        return false;
    }

    diff.normalize();// gradient dcdx

    // Dont adjust position of needle, force mesh to follow needle
    dcdx[0] = Vec3d::Zero();

    // Weight by berycentric coordinates
    dcdx[1] = diff * m_barycentricPt[0];
    dcdx[2] = diff * m_barycentricPt[1];
    dcdx[3] = diff * m_barycentricPt[2];

    return true;
}
