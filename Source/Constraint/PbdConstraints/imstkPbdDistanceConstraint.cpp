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

#include "imstkPbdDistanceConstraint.h"

namespace  imstk
{
void
PbdDistanceConstraint::initConstraint(
    const Vec3d& p0, const Vec3d& p1,
    const BodyVertexId& pIdx0, const BodyVertexId& pIdx1,
    const double k)
{
    m_bodyVertexIds[0] = pIdx0;
    m_bodyVertexIds[1] = pIdx1;
    m_stiffness  = k;
    m_compliance = 1.0 / k;

    m_restLength = (p0 - p1).norm();
}

void
PbdDistanceConstraint::initConstraint(
    const double restLength,
    const BodyVertexId& pIdx0, const BodyVertexId& pIdx1,
    const double k)
{
    m_bodyVertexIds[0] = pIdx0;
    m_bodyVertexIds[1] = pIdx1;
    m_stiffness  = k;
    m_compliance = 1.0 / k;

    m_restLength = restLength;
}

bool
PbdDistanceConstraint::computeValueAndGradient(
    std::vector<PbdBody>& bodies,
    double&               c,
    std::vector<Vec3d>&   dcdx) const
{
    const BodyVertexId& i0 = m_bodyVertexIds[0];
    const BodyVertexId& i1 = m_bodyVertexIds[1];

    const Vec3d& p0 = (*bodies[i0.first].vertices)[i0.second];
    const Vec3d& p1 = (*bodies[i1.first].vertices)[i1.second];

    dcdx[0] = p0 - p1;
    const double len = dcdx[0].norm();
    if (len == 0.0)
    {
        return false;
    }
    dcdx[0] /= len;
    dcdx[1]  = -dcdx[0];
    c        = len - m_restLength;

    return true;
}
} // namespace imstk