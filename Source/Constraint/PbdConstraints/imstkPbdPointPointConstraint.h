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
////
/// \class PbdPointPointConstraint
///
/// \brief Pushes points apart from each other (outside direction defined by p2-p1)
///
class PbdPointPointConstraint : public PbdCollisionConstraint
{
public:
    ///
    /// \brief Constructor
    ///
    PbdPointPointConstraint() : PbdCollisionConstraint(1, 1) { }

    ///
    /// \brief Destructor
    ///
    virtual ~PbdPointPointConstraint() override = default;

public:
    ///
    /// \brief Returns the type of the pbd collision constraint
    ///
    Type getType() const { return Type::PointPoint; }

    ///
    /// \brief initialize constraint
    /// \return
    ///
    void initConstraint(const Side side,
        Vec3d* pt1, double* invMassA1,
        Vec3d* pt2, double* invMassB1,
        std::shared_ptr<PbdCollisionConstraintConfig> configA, std::shared_ptr<PbdCollisionConstraintConfig> configB);

    ///
    /// \brief compute value and gradient of constraint function
    ///
    /// \param[inout] c constraint value
    /// \param[inout] dcdx constraint gradient
    ///
    bool computeValueAndGradient(double& c,
                                 VecDataArray<double, 3>& dcdxA,
                                 VecDataArray<double, 3>& dcdxB) const override;

public:
    //double m_penetrationDepth = 0.0;
    Vec3d  m_normal    = Vec3d::Zero();
    double m_penetrationDepth = 0.0;
};
} // imstk