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
/// \class PbdAngularConstraint
///
/// \brief Applies rotational correction only
///
class PbdAngularConstraint : public PbdConstraint
{
public:
    PbdAngularConstraint(int numParticles) : PbdConstraint(numParticles) { }
    ~PbdAngularConstraint() override = default;

    ///
    /// \brief Update positions by projecting constraints.
    ///
    void projectConstraint(PbdState& bodies, const double dt, const SolverType& type) override;
};
} // namespace imstk