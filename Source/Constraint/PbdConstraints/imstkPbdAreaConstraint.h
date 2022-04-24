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
////
/// \class PbdAreaConstraint
///
/// \brief Area constraint for triangular face
///
class PbdAreaConstraint : public PbdConstraint
{
public:
    PbdAreaConstraint() : PbdConstraint(3) { }

    ///
    /// \brief Initializes the area constraint
    ///
    void initConstraint(
        const Vec3d& p0, const Vec3d& p1, const Vec3d& p2,
        const BodyVertexId& pIdx0, const BodyVertexId& pIdx1, const BodyVertexId& pIdx2,
        const double k = 2.5);

    bool computeValueAndGradient(
        std::vector<PbdBody>& bodies,
        double&               c,
        std::vector<Vec3d>&   dcdx) const override;

public:
    double m_restArea = 0.0;  ///< Area at the rest position
};
} // namespace imstk