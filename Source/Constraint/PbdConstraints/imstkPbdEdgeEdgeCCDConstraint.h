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
/// \class PbdEdgeEdgeCCDConstraint
///
/// \brief Pushes an edge "outside" the other edge
///
class PbdEdgeEdgeCCDConstraint : public PbdCollisionConstraint
{
public:
    PbdEdgeEdgeCCDConstraint() : PbdCollisionConstraint(4, 4) { }
    ~PbdEdgeEdgeCCDConstraint() override = default;

public:
    ///
    /// \brief initialize constraint
    /// \return  true if succeeded
    ///
    void initConstraint(
        const PbdParticleId& prevPtA0, const PbdParticleId& prevPtA1,
        const PbdParticleId& prevPtB0, const PbdParticleId& prevPtB1,
        const PbdParticleId& ptA0, const PbdParticleId& ptA1,
        const PbdParticleId& ptB0, const PbdParticleId& ptB1,
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