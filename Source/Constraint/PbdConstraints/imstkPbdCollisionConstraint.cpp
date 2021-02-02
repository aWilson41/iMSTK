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

#include <iostream>

namespace imstk
{
PbdCollisionConstraint::PbdCollisionConstraint(const unsigned int& n1, const unsigned int& n2)
{
    m_bodiesFirst.resize(n1);
    dcdxA.resize(n1);
    m_bodiesSecond.resize(n2);
    dcdxB.resize(n2);
}

void
PbdCollisionConstraint::projectConstraint()
{
    double c;

    const bool update = this->computeValueAndGradient(c, dcdxA, dcdxB);
    if (!update)
    {
        return;
    }

    double lambda = 0.0;

    // Sum the mass (so we can weight displacements)
    for (size_t i = 0; i < m_bodiesFirst.size(); i++)
    {
        lambda += *m_bodiesFirst[i].invMass * dcdxA[i].squaredNorm();
    }

    for (size_t i = 0; i < m_bodiesSecond.size(); i++)
    {
        lambda += *m_bodiesSecond[i].invMass * dcdxB[i].squaredNorm();
    }

    //std::cout << "lambda: " << lambda << std::endl;

    if (lambda == 0.0)
    {
        return;
    }

    lambda = c / lambda;

    for (size_t i = 0; i < m_bodiesFirst.size(); i++)
    {
        if (*m_bodiesFirst[i].invMass > 0.0)
        {
            (*m_bodiesFirst[i].vertex) -= *m_bodiesFirst[i].invMass * lambda * dcdxA[i] * m_configA->m_stiffness;
        }
    }

    for (size_t i = 0; i < m_bodiesSecond.size(); i++)
    {
        if (*m_bodiesSecond[i].invMass > 0.0)
        {
            Vec3d dx = *m_bodiesSecond[i].invMass * lambda * dcdxB[i] * m_configB->m_stiffness;
            //std::cout << "dx " << i << ": " << dx[0] << ", " << dx[1] << ", " << dx[2] << std::endl;
            (*m_bodiesSecond[i].vertex) -= dx;
        }
    }
    //std::cout << std::endl;
}
} // imstk
