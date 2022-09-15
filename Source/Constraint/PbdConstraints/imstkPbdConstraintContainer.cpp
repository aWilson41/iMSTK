/*
** This file is part of the Interactive Medical Simulation Toolkit (iMSTK)
** iMSTK is distributed under the Apache License, Version 2.0.
** See accompanying NOTICE for details.
*/

#include "imstkPbdConstraintContainer.h"
#include "imstkGraph.h"

#include <unordered_map>

namespace imstk
{
void
PbdConstraintContainer::addConstraint(std::shared_ptr<PbdConstraint> constraint)
{
    m_constraintLock.lock();
    m_constraints.push_back(constraint);
    m_constraintLock.unlock();
}

void
PbdConstraintContainer::removeConstraint(std::shared_ptr<PbdConstraint> constraint)
{
    m_constraintLock.lock();
    iterator i = std::find(m_constraints.begin(), m_constraints.end(), constraint);
    if (i != m_constraints.end())
    {
        m_constraints.erase(i);
    }
    m_constraintLock.unlock();
}

void
PbdConstraintContainer::removeConstraints(std::shared_ptr<std::unordered_set<size_t>> vertices, const int bodyId)
{
    // Remove constraints that contain the given vertices
    auto removeConstraintFunc = [&](std::shared_ptr<PbdConstraint> constraint)
                                {
                                    for (const PbdParticleId& pid : constraint->getParticles())
                                    {
                                        if (pid.first == bodyId && vertices->find(pid.second) != vertices->end())
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };

    m_constraintLock.lock();
    m_constraints.erase(std::remove_if(m_constraints.begin(), m_constraints.end(), removeConstraintFunc),
        m_constraints.end());

    // Also remove partitioned constraints
    for (auto& pc : m_partitionedConstraints)
    {
        pc.erase(std::remove_if(pc.begin(), pc.end(), removeConstraintFunc), pc.end());
    }

    m_constraintLock.unlock();
}

PbdConstraintContainer::iterator
PbdConstraintContainer::eraseConstraint(iterator iter)
{
    m_constraintLock.lock();
    iterator newIter = m_constraints.erase(iter);
    m_constraintLock.unlock();
    return newIter;
}

PbdConstraintContainer::const_iterator
PbdConstraintContainer::eraseConstraint(const_iterator iter)
{
    m_constraintLock.lock();
    const_iterator newIter = m_constraints.erase(iter);
    m_constraintLock.unlock();
    return newIter;
}

void
PbdConstraintContainer::partitionConstraints(const int partitionedThreshold)
{
    // Form the map { vertex : list_of_constraints_involve_vertex }
    std::vector<std::shared_ptr<PbdConstraint>>& allConstraints = m_constraints;

    //std::cout << "---------partitionConstraints: " << allConstraints.size() << std::endl;

    // Determine the maximum body id
    size_t maxBodyId = 0;
    for (size_t constrIdx = 0; constrIdx < allConstraints.size(); ++constrIdx)
    {
        const std::shared_ptr<PbdConstraint>& constr = allConstraints[constrIdx];
        for (const auto& vIds : constr->getParticles())
        {
            maxBodyId = std::max(maxBodyId, static_cast<size_t>(vIds.first));
        }
    }

    // Fully unique hash up until (num particles * maxinum # bodies = max size_t)
    auto getParticleUniqueId = [=](const PbdParticleId& pid)
    {
        return pid.second * maxBodyId + pid.first;
    };

    // Map particle ids -> constraint ids
    std::unordered_map<size_t, std::vector<size_t>> particleIdToConstraintIds;
    for (size_t constrIdx = 0; constrIdx < allConstraints.size(); ++constrIdx)
    {
        const std::shared_ptr<PbdConstraint>& constr = allConstraints[constrIdx];
        for (const PbdParticleId& particleId : constr->getParticles())
        {
            particleIdToConstraintIds[getParticleUniqueId(particleId)].push_back(constrIdx);
        }
    }

    // Add edges to the constraint graph
    // Each edge represent a shared vertex between two constraints
    Graph constraintGraph(allConstraints.size());
    for (const auto& kv : particleIdToConstraintIds)
    {
        // The list of constraints for a particle
        const std::vector<size_t>& constraintIds = kv.second;
        for (size_t i = 0; i < constraintIds.size(); i++)
        {
            for (size_t j = i + 1; j < constraintIds.size(); j++)
            {
                constraintGraph.addEdge(constraintIds[i], constraintIds[j]);
            }
        }
    }
    particleIdToConstraintIds.clear();

    // do graph coloring for the constraint graph
    const Graph::ColorsType coloring = constraintGraph.doColoring(Graph::ColoringMethod::WelshPowell);
    // const auto  coloring = constraintGraph.doColoring(Graph::ColoringMethod::Greedy);
    const std::vector<unsigned short>& partitionIndices = coloring.first;
    const unsigned short numPartitions = coloring.second;
    assert(partitionIndices.size() == allConstraints.size());

    std::vector<std::vector<std::shared_ptr<PbdConstraint>>>& partitionedConstraints = m_partitionedConstraints;
    partitionedConstraints.resize(0);
    partitionedConstraints.resize(static_cast<size_t>(numPartitions));

    for (size_t constraintId = 0; constraintId < partitionIndices.size(); constraintId++)
    {
        const unsigned short partitionIdx = partitionIndices[constraintId];
        partitionedConstraints[partitionIdx].push_back(allConstraints[constraintId]);
    }

    // If a partition has size smaller than the partition threshold, then move its constraints back
    // These constraints will be processed sequentially
    // Because small size partitions yield bad performance upon running in parallel
    allConstraints.resize(0);
    for (const auto& constraints : partitionedConstraints)
    {
        if (constraints.size() < partitionedThreshold)
        {
            for (size_t constraintId = 0; constraintId < constraints.size(); constraintId++)
            {
                allConstraints.push_back(std::move(constraints[constraintId]));
            }
        }
    }

    // Remove all empty partitions
    size_t writeId = 0;
    for (size_t readIdx = 0; readIdx < partitionedConstraints.size(); ++readIdx)
    {
        if (partitionedConstraints[readIdx].size() >= partitionedThreshold)
        {
            if (readIdx != writeId)
            {
                partitionedConstraints[writeId] = std::move(partitionedConstraints[readIdx]);
            }
            writeId++;
        }
    }
    partitionedConstraints.resize(writeId);

    // Print
    /*if (print)
    {
        size_t numConstraints = 0;
        int    idx = 0;
        for (const auto& constraints : partitionedConstraints)
        {
            std::cout << "Partition # " << idx++ << " | # nodes: " << constraints.size() << std::endl;
            numConstraints += constraints.size();
        }
        std::cout << "Sequential processing # nodes: " << allConstraints.size() << std::endl;
        numConstraints += allConstraints.size();
        std::cout << "Total constraints: " << numConstraints << " | Graph size: "
            << constraintGraph.size() << std::endl;
    }*/
}
} // namespace imstk