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

#include "imstkPbdContactConstraint.h"

namespace imstk
{
void
PbdContactConstraint::projectConstraint(PbdState& bodies,
                                        const double dt, const SolverType& solverType)
{
    if (dt == 0.0)
    {
        return;
    }

    double c      = 0.0;
    bool   update = this->computeValueAndGradient(bodies, c,
        m_n, m_r);
    if (!update)
    {
        return;
    }

    double dcMidc  = 0.0;
    double dlambda = 0.0;
    double alpha   = 0.0;

    for (size_t i = 0; i < m_particles.size(); i++)
    {
        const double invMass = bodies.getInvMass(m_particles[i]);
        if (m_contactTypes[i] == ContactType::RIGID)
        {
            const Mat3d& invInteria = bodies.getInvInertia(m_particles[i]);
            const Vec3d  l = m_r[i].cross(m_n[i]);
            dcMidc += l.transpose() * invInteria * l;
        }
        dcMidc += invMass;
    }

    if (dcMidc < IMSTK_DOUBLE_EPS)
    {
        return;
    }

    switch (solverType)
    {
    case (SolverType::PBD):
        dlambda = -c * m_stiffness / dcMidc;
        break;
    case (SolverType::xPBD):
    default:
        alpha     = m_compliance / (dt * dt);
        dlambda   = -(c + alpha * m_lambda) / (dcMidc + alpha);
        m_lambda += dlambda;
        break;
    }

    for (size_t i = 0; i < m_particles.size(); i++)
    {
        const double invMass = bodies.getInvMass(m_particles[i]);
        const Vec3d  p       = dlambda * m_n[i];
        bodies.getPosition(m_particles[i]) += p * invMass;

        if (m_contactTypes[i] == ContactType::RIGID)
        {
            const Mat3d& invInteria = bodies.getInvInertia(m_particles[i]);
            const Vec3d  l = invInteria * (m_r[i].cross(-p));
            //printf("l: %f, %f, %f\n", l[0], l[1], l[2]);
            Quatd&      q  = bodies.getOrientation(m_particles[0]);
            const Quatd dq = Quatd(0.0, l[0] * 0.5, l[1] * 0.5, l[2] * 0.5) * q;
            q.x() += dq.x();
            q.y() += dq.y();
            q.z() += dq.z();
            q.w() += dq.w();
        }
    }
}
} // namespace imstk