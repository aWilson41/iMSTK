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

#include "imstkPbdAreaConstraint.h"

namespace imstk
{
void
PbdAreaConstraint::initConstraint(
    const Vec3d& p0, const Vec3d& p1, const Vec3d& p2,
    const PbdParticleId& pIdx0, const PbdParticleId& pIdx1, const PbdParticleId& pIdx2,
    const double k)
{
    m_particles[0] = pIdx0;
    m_particles[1] = pIdx1;
    m_particles[2] = pIdx2;
    setStiffness(k);

    m_restArea = 0.5 * (p1 - p0).cross(p2 - p0).norm();
}

bool
PbdAreaConstraint::computeValueAndGradient(PbdState& bodies,
                                           double& c, std::vector<Vec3d>& dcdx) const
{
    const Vec3d& p0 = bodies.getPosition(m_particles[0]);
    const Vec3d& p1 = bodies.getPosition(m_particles[1]);
    const Vec3d& p2 = bodies.getPosition(m_particles[2]);

    const Vec3d e0 = p0 - p1;
    const Vec3d e1 = p1 - p2;
    const Vec3d e2 = p2 - p0;

    Vec3d n = e0.cross(e1);
    c = 0.5 * n.norm();

    if (c < m_epsilon)
    {
        return false;
    }

    n /= 2.0 * c;
    c -= m_restArea;

    dcdx[0] = e1.cross(n);
    dcdx[1] = e2.cross(n);
    dcdx[2] = e0.cross(n);

    return true;
}
} // namespace imstk
