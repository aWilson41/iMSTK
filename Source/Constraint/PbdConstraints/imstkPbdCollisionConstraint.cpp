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

#include "imstkPbdCollisionConstraint.h"

namespace imstk
{
PbdCollisionConstraint::PbdCollisionConstraint(const int numParticlesA, const int numParticlesB) :
    PbdConstraint(numParticlesA + numParticlesB)
{
    m_bodiesSides.resize(numParticlesA + numParticlesB);
    for (int i = 0; i < m_bodiesSides.size(); i++)
    {
        m_bodiesSides[i] = (i >= numParticlesA); // false/0 for A, true/1 for B
    }
}

void
PbdCollisionConstraint::projectConstraint(PbdState& bodies, const double dt, const SolverType& type)
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

    double lambda = 0.0;

    // Sum the mass (so we can weight displacements)
    for (size_t i = 0; i < m_particles.size(); i++)
    {
        lambda += bodies.getInvMass(m_particles[i]) * m_dcdx[i].squaredNorm();
    }

    if (lambda == 0.0)
    {
        return;
    }

    lambda = c / lambda;

    size_t vertexId = 0;
    size_t bodyId   = 0;
    for (size_t i = 0; i < m_particles.size(); i++)
    {
        const double invMass = bodies.getInvMass(m_particles[i]);
        if (invMass > 0.0)
        {
            bodies.getPosition(m_particles[i]) += invMass * lambda *
                                                  m_dcdx[i] * m_stiffness[m_bodiesSides[i]];
        }
    }
}
} // namespace imstk