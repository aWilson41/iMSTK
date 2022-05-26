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

#include "imstkPbdPointEdgeConstraint.h"

namespace imstk
{
void
PbdPointEdgeConstraint::initConstraint(
    const PbdParticleId& ptA1,
    const PbdParticleId& ptB1, const PbdParticleId& ptB2,
    double stiffnessA, double stiffnessB)
{
    m_particles[0] = ptA1;
    m_particles[1] = ptB1;
    m_particles[2] = ptB2;

    m_stiffness[0] = stiffnessA;
    m_stiffness[1] = stiffnessB;
}

bool
PbdPointEdgeConstraint::computeValueAndGradient(PbdState& bodies,
                                                double& c, std::vector<Vec3d>& dcdx) const
{
    // Just project x0 onto x3-x2. Get the normal component for distance to line
    const Vec3d& x0 = bodies.getPosition(m_particles[0]);
    const Vec3d& x1 = bodies.getPosition(m_particles[1]);
    const Vec3d& x2 = bodies.getPosition(m_particles[2]);

    const Vec3d  ab     = x2 - x1;
    const double length = ab.norm();
    if (length == 0.0)
    {
        // There is no distance between the edge, can't do anything
        c = 0.0;
        return false;
    }
    const Vec3d dir1 = ab / length;

    // Project onto the line
    const Vec3d  diff = x0 - x1;
    const double p    = dir1.dot(diff);
    if (p < 0.0 || p > length)
    {
        c = 0.0;
        return false;
    }
    // Remove tangent component to get normal
    const Vec3d  diff1 = diff - p * dir1;
    const double l     = diff1.norm();
    if (l == 0.0)
    {
        // The point is on the line
        c = 0.0;
        return false;
    }
    const Vec3d  n = diff1 / l;
    const double u = p / length;

    dcdx[0] = -n;
    dcdx[1] = (1.0 - u) * n;
    dcdx[2] = u * n;

    c = l;

    return true;
}
} // namespace imstk