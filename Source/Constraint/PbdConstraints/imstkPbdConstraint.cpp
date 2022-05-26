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

#include "imstkPbdConstraint.h"

namespace imstk
{
void
PbdConstraint::projectConstraint(PbdState& bodies,
                                 const double dt, const SolverType& solverType)
{
    if (dt == 0.0)
    {
        return;
    }

    double c      = 0.0;
    bool   update = this->computeValueAndGradient(bodies, c, m_dcdx);
    if (!update)
    {
        return;
    }

    double dcMidc  = 0.0;
    double dlambda = 0.0;
    double alpha   = 0.0;

    for (size_t i = 0; i < m_particles.size(); i++)
    {
        dcMidc += bodies.getInvMass(m_particles[i]) * m_dcdx[i].squaredNorm();
    }

    if (dcMidc < IMSTK_DOUBLE_EPS)
    {
        return;
    }

    switch (solverType)
    {
    case (SolverType::xPBD):
        alpha     = m_compliance / (dt * dt);
        dlambda   = -(c + alpha * m_lambda) / (dcMidc + alpha);
        m_lambda += dlambda;
        break;
    case (SolverType::PBD):
        dlambda = -c * m_stiffness / dcMidc;
        break;
    default:
        alpha     = m_compliance / (dt * dt);
        dlambda   = -(c + alpha * m_lambda) / (dcMidc + alpha);
        m_lambda += dlambda;
    }

    for (size_t i = 0; i < m_particles.size(); i++)
    {
        const double invMass = bodies.getInvMass(m_particles[i]);
        if (invMass > 0.0)
        {
            bodies.getPosition(m_particles[i]) += invMass * dlambda * m_dcdx[i];
        }
    }
}

void
PbdConstraint::correctVelocity(PbdState& bodies)
{
    const double fricFrac = 1.0 - m_friction;

    for (size_t i = 0; i < m_particles.size(); i++)
    {
        size_t       vertexId = 0;
        size_t       bodyId   = 0;
        const double invMass  = bodies.getInvMass(m_particles[i]);
        if (invMass > 0.0)
        {
            const Vec3d n = m_dcdx[i].normalized();
            Vec3d&      v = bodies.getVelocity(m_particles[i]);

            // Separate velocity into normal and tangent components
            const Vec3d vN = n.dot(v) * n;
            const Vec3d vT = v - vN;

            // Put back together fractionally based on defined restitution and frictional coefficients
            v = vN * m_restitution + vT * fricFrac;
        }
    }
}
} // namespace imstk
