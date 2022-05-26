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

#include "imstkPbdSolver.h"
#include "imstkLogger.h"
#include "imstkParallelUtils.h"
#include "imstkPbdCollisionConstraint.h"
#include "imstkPbdConstraintContainer.h"

namespace imstk
{
PbdSolver::PbdSolver() :
    m_constraints(std::make_shared<PbdConstraintContainer>()),
    m_constraintLists(std::make_shared<std::list<std::vector<PbdConstraint*>*>>())
{
}

void
PbdSolver::solve()
{
    const std::vector<std::shared_ptr<PbdConstraint>>&              constraints = m_constraints->getConstraints();
    const std::vector<std::vector<std::shared_ptr<PbdConstraint>>>& partitionedConstraints = m_constraints->getPartitionedConstraints();

    // zero out the Lagrange multiplier
    for (const auto& constraint : constraints)
    {
        constraint->zeroOutLambda();
    }

    for (const auto& constraintPartition : partitionedConstraints)
    {
        ParallelUtils::parallelFor(constraintPartition.size(),
            [&](const size_t idx)
            {
                constraintPartition[idx]->zeroOutLambda();
            });
    }

    unsigned int i = 0;
    while (i++ < m_iterations)
    {
        for (const auto& constraint : constraints)
        {
            constraint->projectConstraint(*m_state, m_dt, m_solverType);
        }

        for (const auto& constraintPartition : partitionedConstraints)
        {
            ParallelUtils::parallelFor(constraintPartition.size(),
                [&](const size_t idx)
                {
                    constraintPartition[idx]->projectConstraint(*m_state, m_dt, m_solverType);
                });
            //// Sequential
            //for (size_t k = 0; k < constraintPartition.size(); k++)
            //{
            //    constraintPartition[k]->projectConstraint(invMasses, m_dt, m_solverType, currPositions);
            //}
        }

        for (auto constraintList : *m_constraintLists)
        {
            const std::vector<PbdConstraint*>& constraints = *constraintList;
            for (size_t j = 0; j < constraints.size(); j++)
            {
                constraints[j]->projectConstraint(*m_state, m_dt, m_solverType);
            }
        }
    }
}
} // namespace imstk