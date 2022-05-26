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
#include "imstkRbdConstraint.h"

using namespace imstk;

///
/// \class ThreadInsertionConstraint
///
/// \brief Constrains an point on a line mesh to a PBD surface mesh
///
class ThreadInsertionConstraint : public PbdCollisionConstraint
{
public:
    ThreadInsertionConstraint() :  PbdCollisionConstraint(2, 3) { }
    ~ThreadInsertionConstraint() override = default;

    void initConstraint(
        const PbdState&      bodies,
        const PbdParticleId& ptA1,
        const PbdParticleId& ptA2,
        const Vec2d&         threadBaryPoint,
        const PbdParticleId& ptB1,
        const PbdParticleId& ptB2,
        const PbdParticleId& ptB3,
        const Vec3d&         triBaryPoint,
        double               stiffnessA,
        double               stiffnessB);

    bool computeValueAndGradient(PbdState& bodies,
                                 double& c, std::vector<Vec3d>& dcdx) const override;

protected:
    Vec2d m_threadBaryPt;
    Vec3d m_triangleBaryPt;
    Vec3d m_triInsertionPoint;
    Vec3d m_threadInsertionPoint;
};