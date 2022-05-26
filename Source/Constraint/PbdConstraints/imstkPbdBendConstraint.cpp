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

#include "imstkPbdBendConstraint.h"

namespace imstk
{
void
PbdBendConstraint::initConstraint(
    const Vec3d& p0, const Vec3d& p1, const Vec3d& p2,
    const PbdParticleId& pIdx0, const PbdParticleId& pIdx1, const PbdParticleId& pIdx2,
    const double k)
{
    // Instead of using the angle between the segments we can use the distance
    // from the center of the triangle
    const Vec3d& center = (p0 + p1 + p2) / 3.0;

    initConstraint(pIdx0, pIdx1, pIdx2, (p1 - center).norm(), k);
}

void
PbdBendConstraint::initConstraint(
    const PbdParticleId& pIdx0, const PbdParticleId& pIdx1, const PbdParticleId& pIdx2,
    const double restLength,
    const double k)
{
    m_particles[0] = pIdx0;
    m_particles[1] = pIdx1;
    m_particles[2] = pIdx2;

    setStiffness(k);

    m_restLength = restLength;
}

bool
PbdBendConstraint::computeValueAndGradient(PbdState& bodies,
                                           double& c, std::vector<Vec3d>& dcdx) const
{
    const Vec3d& p0 = bodies.getPosition(m_particles[0]);
    const Vec3d& p1 = bodies.getPosition(m_particles[1]);
    const Vec3d& p2 = bodies.getPosition(m_particles[2]);

    // Move towards triangle center
    const Vec3d& center = (p0 + p1 + p2) / 3.0;
    const Vec3d& diff   = p1 - center;
    const double dist   = diff.norm();

    if (dist < m_epsilon)
    {
        return false;
    }

    c = dist - m_restLength;

    dcdx[0] = (-2.0 / dist) * diff;
    dcdx[1] = -2.0 * dcdx[0];
    dcdx[2] = dcdx[0];

    return true;
}
} // namespace imstk