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
/// \class PbdDistanceConstraint
///
/// \brief Distance constraints between two nodal points
///
class PbdDistanceConstraint : public PbdConstraint
{
public:
    PbdDistanceConstraint() : PbdConstraint(2) { }

    ///
    /// \brief Initializes the distance constraint
    ///@{
    void initConstraint(
        const Vec3d& p0, const Vec3d& p1,
        const BodyVertexId& pIdx0, const BodyVertexId& pIdx1,
        const double k = 1e5);
    void initConstraint(const double restLength,
                        const BodyVertexId& pIdx0, const BodyVertexId& pIdx1,
                        const double k = 1e5);
    ///@}

    bool computeValueAndGradient(
        std::vector<PbdBody>& bodies,
        double&               c,
        std::vector<Vec3d>&   dcdx) const override;

public:
    double m_restLength = 0.0; ///< Rest length between the nodes
};
} // namespace imstk