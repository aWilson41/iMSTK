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

#include "imstkPbdVolumeConstraint.h"

namespace imstk
{
void
PbdVolumeConstraint::initConstraint(
    const Vec3d& p0, const Vec3d& p1, const Vec3d& p2, const Vec3d& p3,
    const PbdParticleId& pIdx0, const PbdParticleId& pIdx1,
    const PbdParticleId& pIdx2, const PbdParticleId& pIdx3,
    const double k)
{
    m_particles[0] = pIdx0;
    m_particles[1] = pIdx1;
    m_particles[2] = pIdx2;
    m_particles[3] = pIdx3;

    setStiffness(k);

    m_restVolume = (1.0 / 6.0) * ((p1 - p0).cross(p2 - p0)).dot(p3 - p0);
}

bool
PbdVolumeConstraint::computeValueAndGradient(PbdState& bodies,
                                             double& c, std::vector<Vec3d>& dcdx) const
{
    const Vec3d& x0 = bodies.getPosition(m_particles[0]);
    const Vec3d& x1 = bodies.getPosition(m_particles[1]);
    const Vec3d& x2 = bodies.getPosition(m_particles[2]);
    const Vec3d& x3 = bodies.getPosition(m_particles[3]);

    const double onesixth = 1.0 / 6.0;

    dcdx[0] = onesixth * (x1 - x2).cross(x3 - x1);
    dcdx[1] = onesixth * (x2 - x0).cross(x3 - x0);
    dcdx[2] = onesixth * (x3 - x0).cross(x1 - x0);
    dcdx[3] = onesixth * (x1 - x0).cross(x2 - x0);

    const double volume = dcdx[3].dot(x3 - x0);
    c = volume - m_restVolume;
    return true;
}
} // namespace imstk