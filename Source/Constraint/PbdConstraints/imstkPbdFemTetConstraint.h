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

#include "imstkPbdFemConstraint.h"

namespace imstk
{
///
/// \class PbdFemTetConstraint
///
/// \brief The FEMTetConstraint class class for constraint as the elastic energy
/// computed by linear shape functions with tetrahedral mesh.
///
class PbdFemTetConstraint : public PbdFemConstraint
{
public:
    PbdFemTetConstraint(MaterialType mType = MaterialType::StVK) :
        PbdFemConstraint(4, mType) { }

    ///
    /// \brief Initialize the constraint
    ///
    bool initConstraint(
        const Vec3d& p0, const Vec3d& p1, const Vec3d& p2, const Vec3d& p3,
        const PbdParticleId& pIdx0, const PbdParticleId& pIdx1,
        const PbdParticleId& pIdx2, const PbdParticleId& pIdx3,
        std::shared_ptr<PbdFemConstraintConfig> config);

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