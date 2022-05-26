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

#include "imstkPbdPointPointConstraint.h"

namespace imstk
{
void
PbdPointPointConstraint::initConstraint(
    const PbdParticleId& ptA, const PbdParticleId& ptB,
    double stiffnessA, double stiffnessB)
{
    m_particles[0] = ptA;
    m_particles[1] = ptB;

    m_stiffness[0] = stiffnessA;
    m_stiffness[1] = stiffnessB;
}

bool
PbdPointPointConstraint::computeValueAndGradient(PbdState& bodies,
                                                 double& c, std::vector<Vec3d>& dcdx) const
{
    // Current position during solve
    const Vec3d& x0 = bodies.getPosition(m_particles[0]);
    const Vec3d& x1 = bodies.getPosition(m_particles[1]);

    const Vec3d diff = x1 - x0;
    c = diff.norm();

    if (c == 0.0)
    {
        return false;
    }

    const Vec3d n = diff / c;

    // A
    dcdx[0] = n;
    // B
    dcdx[1] = -n;

    return true;
}
} // namespace imstk