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
    const BodyVertexId& pIdx0, const BodyVertexId& pIdx1, const BodyVertexId& pIdx2,
    const double k)
{
    m_bodyVertexIds[0] = pIdx0;
    m_bodyVertexIds[1] = pIdx1;
    m_bodyVertexIds[2] = pIdx2;

    this->m_stiffness  = k;
    this->m_compliance = 1.0 / k;

    m_restArea = 0.5 * (p1 - p0).cross(p2 - p0).norm();
}

bool
PbdAreaConstraint::computeValueAndGradient(
    std::vector<PbdBody>& bodies,
    double&               c,
    std::vector<Vec3d>&   dcdx) const
{
    const BodyVertexId& i0 = m_bodyVertexIds[0];
    const BodyVertexId& i1 = m_bodyVertexIds[1];
    const BodyVertexId& i2 = m_bodyVertexIds[2];

    const Vec3d& p0 = (*bodies[i0.first].vertices)[i0.second];
    const Vec3d& p1 = (*bodies[i1.first].vertices)[i1.second];
    const Vec3d& p2 = (*bodies[i2.first].vertices)[i2.second];

    const Vec3d e0 = p0 - p1;
    const Vec3d e1 = p1 - p2;
    const Vec3d e2 = p2 - p0;

    Vec3d n = e0.cross(e1);
    c = 0.5 * n.norm();

    if (c < m_epsilon)
    {
        return false;
    }

    n /= 2 * c;
    c -= m_restArea;

    dcdx[0] = e1.cross(n);
    dcdx[1] = e2.cross(n);
    dcdx[2] = e0.cross(n);

    return true;
}
} // namespace imstk
