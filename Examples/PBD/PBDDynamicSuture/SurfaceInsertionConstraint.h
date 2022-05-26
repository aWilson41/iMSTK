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

using namespace imstk;

///
/// \class SurfaceInsertionConstraint
///
/// \brief Constrains a barycentric point on a surface mesh to a rigid body arc needle
///
class SurfaceInsertionConstraint : public PbdCollisionConstraint
{
public:
    SurfaceInsertionConstraint() :  PbdCollisionConstraint(1, 3) { }
    ~SurfaceInsertionConstraint() override = default;

    void initConstraint(
        const Vec3d&         insertionPoint,
        const PbdParticleId& ptB1,
        const PbdParticleId& ptB2,
        const PbdParticleId& ptB3,
        const Vec3d&         contactPt,
        const Vec3d&         barycentricPt,
        double               stiffnessA,
        double               stiffnessB);

    bool computeValueAndGradient(PbdState& bodies,
                                 double& c, std::vector<Vec3d>& dcdx) const override;

private:
    Vec3d m_insertionPoint;
    Vec3d m_barycentricPt;
    Vec3d m_contactPt;
};