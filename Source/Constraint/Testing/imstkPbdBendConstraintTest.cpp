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

#include "imstkPbdBendConstraint.h"
#include "imstkPbdConstraintTest.h"

#include <gtest/gtest.h>

using namespace imstk;

///
/// \brief Test that two connecting line segments unfold
///
TEST_F(PbdConstraintTest, BendConstraint_TestConvergence1)
{
    setNumParticles(3);

    // Straight line upon initialization
    m_vertices[0] = Vec3d(0.0, 0.0, 0.0);
    m_vertices[1] = Vec3d(0.5, 0.0, 0.0);
    m_vertices[2] = Vec3d(1.0, 0.0, 0.0);

    m_invMasses[0] = 1.0;
    m_invMasses[1] = 0.0; // Center doesn't move
    m_invMasses[2] = 1.0;

    PbdBendConstraint constraint;
    m_constraint = &constraint;
    constraint.initConstraint(
        m_vertices[0], m_vertices[1], m_vertices[2],
        { 0, 0 }, { 0, 1 }, { 0, 2 }, 1e20);

    // Modify it so the line segments look like \/
    m_vertices[0][1] = 0.1;
    m_vertices[2][1] = 0.1;
    for (int i = 0; i < 500; i++)
    {
        solve(0.01, PbdConstraint::SolverType::xPBD);
    }

    // Should resolve back to a flat line
    EXPECT_NEAR(m_vertices[0][1], 0.0, IMSTK_DOUBLE_EPS);
    EXPECT_NEAR(m_vertices[2][1], 0.0, IMSTK_DOUBLE_EPS);
}