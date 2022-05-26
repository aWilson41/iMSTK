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

#include "imstkPbdCollisionConstraint.h"

namespace imstk
{
///
/// \class PbdEdgeEdgeConstraint
///
/// \brief Resolves an edge given by two pbd particles to coincide with
/// another edge
///
class PbdEdgeEdgeConstraint : public PbdCollisionConstraint
{
public:
    PbdEdgeEdgeConstraint() : PbdCollisionConstraint(2, 2) { }
    ~PbdEdgeEdgeConstraint() override = default;

public:
    ///
    /// \brief Initialize constraint
    /// \param First point of edge A
    /// \param Second point of edge A
    /// \param First point of edge B
    /// \param Second point of edge B
    /// \param Stiffness of A
    /// \param Stiffness of B
    ///
    void initConstraint(
        const PbdParticleId& ptA1, const PbdParticleId& ptA2,
        const PbdParticleId& ptB1, const PbdParticleId& ptB2,
        double stiffnessA, double stiffnessB);

    ///
    /// \brief Compute value and gradient of constraint function
    /// \param[inout] set of bodies involved in system
    /// \param[inout] c constraint value
    /// \param[inout] dcdx constraint gradient
    ///
    bool computeValueAndGradient(PbdState& bodies,
                                 double& c, std::vector<Vec3d>& dcdx) const override;
};
} // namespace imstk