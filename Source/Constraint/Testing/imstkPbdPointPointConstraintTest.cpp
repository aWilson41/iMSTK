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

#include "imstkPbdConstraintTest.h"
#include "imstkPbdPointPointConstraint.h"

#include <gtest/gtest.h>

using namespace imstk;

///
/// \brief Test that two points meet
///
TEST_F(PbdConstraintTest, PointPointConstraint_TestConvergence1)
{
    setNumParticles(2);

    m_invMasses.fill(1.0);

    m_vertices[0] = Vec3d(0.0, 0.0, 0.0);
    m_vertices[1] = Vec3d(0.0, -1.0, 0.0);

    PbdPointPointConstraint constraint;
    m_constraint = &constraint;
    constraint.initConstraint(
        { 0, 0 }, { 0, 1 },
        1.0, 1.0);
    for (int i = 0; i < 3; i++)
    {
        solve(0.01, PbdConstraint::SolverType::PBD);
    }

    ASSERT_EQ(m_vertices[0][1], m_vertices[1][1]);
}