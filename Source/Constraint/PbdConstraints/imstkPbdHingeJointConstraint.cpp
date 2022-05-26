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

#include "imstkPbdHingeJointConstraint.h"

namespace imstk
{
void
PbdHingeJointConstraint::initConstraint(
    const PbdParticleId& pIdx0,
    const Vec3d&         hingeAxes,
    const double         k)
{
    m_particles[0] = pIdx0;
    setStiffness(k);

    m_hingeAxes = hingeAxes;
}

bool
PbdHingeJointConstraint::computeValueAndGradient(PbdState& bodies,
                                                 double& c, std::vector<Vec3d>& dcdx) const
{
    // Is this the fastest way to get basis vectors?
    const Vec3d up = bodies.getOrientation(m_particles[0]).toRotationMatrix().col(1);

    // Gives rotation
    Vec3d dir = m_hingeAxes.cross(up);
    dcdx[0] = dir.normalized();
    c       = dir.norm();

    return true;
}
} // namespace imstk