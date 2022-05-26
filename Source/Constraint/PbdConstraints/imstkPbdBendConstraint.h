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

#include "imstkPbdConstraint.h"

namespace imstk
{
///
/// \class PbdBendConstraint
///
/// \brief Bend constraint between two segments. Maintains angle between
/// the two segments in the initial configuration given
///
class PbdBendConstraint : public PbdConstraint
{
public:
    PbdBendConstraint() : PbdConstraint(3) { }

    /**
        \brief Initialize the constraint
            p0
               \
                \
                p1
                /
               /
            p2
    */
    void initConstraint(
        const Vec3d& initPos0, const Vec3d& initPos1, const Vec3d& initPos2,
        const PbdParticleId& pIdx0, const PbdParticleId& pIdx1, const PbdParticleId& pIdx2,
        const double k);
    void initConstraint(
        const PbdParticleId& pIdx0, const PbdParticleId& pIdx1, const PbdParticleId& pIdx2,
        const double restLength,
        const double k);

    ///
    /// \brief Compute value and gradient of constraint function
    /// \param[inout] set of bodies involved in system
    /// \param[inout] c constraint value
    /// \param[inout] dcdx constraint gradient
    ///
    bool computeValueAndGradient(PbdState& bodies,
                                 double& c, std::vector<Vec3d>& dcdx) const override;

public:
    double m_restLength = 0.; ///< Rest length
};
} // namespace imstk