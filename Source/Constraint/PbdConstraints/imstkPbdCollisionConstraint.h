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

#pragma once

#include "imstkMath.h"
#include "imstkPbdConstraint.h"
#include "imstkVecDataArray.h"

namespace imstk
{
///
/// \class PbdCollisionConstraint
///
/// \brief The PbdCollisionConstraint implements two sided collision.
/// This allows the usage of differing stiffness for each side which can
/// be useful during solve. Due to differences in definition, collisions
/// do not use XPBD. Only PBD. They are assumed perfectly rigid even though
/// stiffness is modifiable. Given enough iterations in the solve, it will
/// converge to perfectly rigid.
///
/// The collision constraint also provides a correctVelocity function.
/// This may be overriden but by default it will correct velocity along
/// the gradient tangents and normal according to frictional and restitution
/// coefficients.
///
class PbdCollisionConstraint : public PbdConstraint
{
public:
    ~PbdCollisionConstraint() override = default;

public:
    ///
    /// \brief Get/Set stiffness A or B
    ///@{
    double getStiffnessA() const { return m_stiffness[0]; }
    void setStiffnessA(const double stiffnessA) { m_stiffness[0] = stiffnessA; }
    double getStiffnessB() const { return m_stiffness[1]; }
    void setStiffnessB(const double stiffnessB) { m_stiffness[1] = stiffnessB; }
    ///@}

    ///
    /// \brief Performs the actual positional solve
    ///
    void projectConstraint(PbdState& bodies,
                           const double dt, const SolverType& type) override;

protected:
    PbdCollisionConstraint(const int numParticlesA, const int numParticlesB);

    std::vector<bool> m_bodiesSides; ///< Stores 0 or 1 to indicate side of particle
    double m_stiffness[2] = { 1.0, 1.0 };
};
} // namespace imstk
