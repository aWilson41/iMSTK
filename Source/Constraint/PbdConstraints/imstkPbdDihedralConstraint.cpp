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

#include "imstkPbdDihedralConstraint.h"

namespace  imstk
{
void
PbdDihedralConstraint::initConstraint(
    const Vec3d& p0, const Vec3d& p1, const Vec3d& p2, const Vec3d& p3,
    const BodyVertexId& pIdx0, const BodyVertexId& pIdx1,
    const BodyVertexId& pIdx2, const BodyVertexId& pIdx3,
    const double k)
{
    m_bodyVertexIds[0] = pIdx0;
    m_bodyVertexIds[1] = pIdx1;
    m_bodyVertexIds[2] = pIdx2;
    m_bodyVertexIds[3] = pIdx3;

    m_stiffness  = k;
    m_compliance = 1.0 / k;

    const Vec3d n1 = (p2 - p0).cross(p3 - p0).normalized();
    const Vec3d n2 = (p3 - p1).cross(p2 - p1).normalized();

    m_restAngle = atan2(n1.cross(n2).dot(p3 - p2), (p3 - p2).norm() * n1.dot(n2));
}

bool
PbdDihedralConstraint::computeValueAndGradient(
    std::vector<PbdBody>& bodies,
    double&               c,
    std::vector<Vec3d>&   dcdx) const
{
    const BodyVertexId& i0 = m_bodyVertexIds[0];
    const BodyVertexId& i1 = m_bodyVertexIds[1];
    const BodyVertexId& i2 = m_bodyVertexIds[2];
    const BodyVertexId& i3 = m_bodyVertexIds[3];

    const Vec3d& p0 = (*bodies[i0.first].vertices)[i0.second];
    const Vec3d& p1 = (*bodies[i1.first].vertices)[i1.second];
    const Vec3d& p2 = (*bodies[i2.first].vertices)[i2.second];
    const Vec3d& p3 = (*bodies[i3.first].vertices)[i3.second];

    const Vec3d e  = p3 - p2;
    const Vec3d e1 = p3 - p0;
    const Vec3d e2 = p0 - p2;
    const Vec3d e3 = p3 - p1;
    const Vec3d e4 = p1 - p2;
    // To accelerate, all normal (area) vectors and edge length should be precomputed in parallel
    Vec3d        n1 = e1.cross(e);
    Vec3d        n2 = e.cross(e3);
    const double A1 = n1.norm();
    const double A2 = n2.norm();
    n1 /= A1;
    n2 /= A2;

    const double l = e.norm();
    if (l < m_epsilon)
    {
        return false;
    }

    dcdx[0] = -(l / A1) * n1;
    dcdx[1] = -(l / A2) * n2;
    dcdx[2] = (e.dot(e1) / (A1 * l)) * n1 + (e.dot(e3) / (A2 * l)) * n2;
    dcdx[3] = (e.dot(e2) / (A1 * l)) * n1 + (e.dot(e4) / (A2 * l)) * n2;

    c = atan2(n1.cross(n2).dot(e), l * n1.dot(n2)) - m_restAngle;

    return true;
}
} // namespace imstk