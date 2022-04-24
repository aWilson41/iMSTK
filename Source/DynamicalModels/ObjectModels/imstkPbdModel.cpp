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

#include "imstkPbdModel.h"
#include "imstkGraph.h"
#include "imstkLineMesh.h"
#include "imstkLogger.h"
#include "imstkParallelUtils.h"
#include "imstkPbdSolver.h"
#include "imstkSurfaceMesh.h"
#include "imstkTaskGraph.h"
#include "imstkTetrahedralMesh.h"
#include "imstkPbdConstraintFunctor.h"

namespace imstk
{
void
PbdModelConfig::computeElasticConstants()
{
    if (std::abs(m_femParams->m_mu) < std::numeric_limits<double>::min()
        && std::abs(m_femParams->m_lambda) < std::numeric_limits<double>::min())
    {
        const double E  = m_femParams->m_YoungModulus;
        const double nu = m_femParams->m_PoissonRatio;
        m_femParams->m_mu     = E / 2.0 / (1.0 + nu);
        m_femParams->m_lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    }
    else
    {
        const double mu     = m_femParams->m_mu;
        const double lambda = m_femParams->m_lambda;
        m_femParams->m_YoungModulus = mu * (3.0 * lambda + 2.0 * mu) / (lambda + mu);
        m_femParams->m_PoissonRatio = lambda / 2.0 / (lambda + mu);
    }
}

void
PbdModelConfig::enableConstraint(ConstraintGenType type, double stiffness, const int bodyId)
{
    auto& funcs = m_functors[type];
    if (type == ConstraintGenType::Distance)
    {
        auto functor = std::make_shared<PbdDistanceConstraintFunctor>();
        funcs.push_back(functor);
        functor->setStiffness(stiffness);
        functor->setBodyIndex(bodyId);
    }
    else if (type == ConstraintGenType::Volume)
    {
        auto functor = std::make_shared<PbdVolumeConstraintFunctor>();
        funcs.push_back(functor);
        functor->setStiffness(stiffness);
        functor->setBodyIndex(bodyId);
    }
    else if (type == ConstraintGenType::Area)
    {
        auto functor = std::make_shared<PbdAreaConstraintFunctor>();
        funcs.push_back(functor);
        functor->setStiffness(stiffness);
        functor->setBodyIndex(bodyId);
    }
    else if (type == ConstraintGenType::Bend)
    {
        auto functor = std::make_shared<PbdBendConstraintFunctor>();
        funcs.push_back(functor);
        functor->setStiffness(stiffness);
        functor->setStride(1);
        functor->setBodyIndex(bodyId);
    }
    else if (type == ConstraintGenType::Dihedral)
    {
        auto functor = std::make_shared<PbdDihedralConstraintFunctor>();
        funcs.push_back(functor);
        functor->setStiffness(stiffness);
        functor->setBodyIndex(bodyId);
    }
    else if (type == ConstraintGenType::ConstantDensity)
    {
        auto functor = std::make_shared<PbdConstantDensityConstraintFunctor>();
        funcs.push_back(functor);
        functor->setStiffness(stiffness);
        functor->setBodyIndex(bodyId);
    }
    else
    {
        LOG(FATAL) << "There exists no standard constraint functor for the ConstraintGenType";
    }
}

void
PbdModelConfig::enableBendConstraint(const double stiffness, const int stride, const bool restLength0, const int bodyId)
{
    auto& funcs = m_functors[ConstraintGenType::Bend];

    // Find the functor with the same stride
    std::shared_ptr<PbdBendConstraintFunctor> foundFunctor = nullptr;
    for (auto functor : funcs)
    {
        auto bendFunctor = std::dynamic_pointer_cast<PbdBendConstraintFunctor>(functor);
        if (bendFunctor->getStride() == stride)
        {
            foundFunctor = bendFunctor;
            break;
        }
    }

    // If one with stride not found, create our own
    if (foundFunctor == nullptr)
    {
        foundFunctor = std::make_shared<PbdBendConstraintFunctor>();
        funcs.push_back(foundFunctor);
    }

    foundFunctor->setRestLength(restLength0 ? 0.0 : -1.0);
    foundFunctor->setBodyIndex(bodyId);
    foundFunctor->setStiffness(stiffness);
    foundFunctor->setStride(stride);
}

void
PbdModelConfig::enableFemConstraint(PbdFemConstraint::MaterialType material, const int bodyId)
{
    auto& funcs = m_functors[ConstraintGenType::FemTet];
    if (funcs.size() == 0)
    {
        funcs.push_back(std::make_shared<PbdFemTetConstraintFunctor>());
    }
    auto functor = std::dynamic_pointer_cast<PbdFemTetConstraintFunctor>(funcs.front());
    functor->setBodyIndex(bodyId);
    functor->setFemConfig(m_femParams);
    functor->setMaterialType(material);
}

PbdModel::PbdModel() : DynamicalModel(DynamicalModelType::PositionBasedDynamics),
    m_config(std::make_shared<PbdModelConfig>())
{
    m_validGeometryTypes = {
        "PointSet",
        "LineMesh",
        "SurfaceMesh",
        "TetrahedralMesh",
        "HexahedralMesh"
    };

    // Setup PBD compute nodes
    m_integrationPositionNode = m_taskGraph->addFunction("PbdModel_IntegratePosition",
        [&]() { integratePosition(); });
    m_solveConstraintsNode = m_taskGraph->addFunction("PbdModel_SolveConstraints",
        [&]() { solveConstraints(); });
    m_updateVelocityNode = m_taskGraph->addFunction("PbdModel_UpdateVelocity",
        [&]() { updateVelocity(); });
}

void
PbdModel::configure(std::shared_ptr<PbdModelConfig> config)
{
    m_config = config;
}

bool
PbdModel::initialize()
{
    getCurrentState()->initialize();
    // Setup the initial state
    getInitialState()->setState(getCurrentState());
    getPreviousState()->setState(getCurrentState());
    // Use ptr reference on initial state to the initial positions in the mesh
    for (auto body : getInitialState()->m_bodies)
    {
        body->vertices = getInitialState()->getBodyGeometry(*body)->getInitialVertexPositions();
    }

    // Initialize constraints
    {
        m_constraints = std::make_shared<PbdConstraintContainer>();

        m_config->computeElasticConstants();

        // Run the functor for every body
        for (const auto& functorVec : m_config->m_functors)
        {
            for (const auto& functor : functorVec.second)
            {
                // If a body functor, run on a specific body
                if (auto bodyFunctor = std::dynamic_pointer_cast<PbdBodyConstraintFunctor>(functor))
                {
                    for (const auto& body : getCurrentState()->m_bodies)
                    {
                        // Only execute if this is that body
                        if (bodyFunctor->m_bodyIndex == body->bodyHandle)
                        {
                            bodyFunctor->setGeometry(getCurrentState()->getBodyGeometry(*body));
                            (*bodyFunctor)(*m_constraints);
                        }
                    }
                }
                // If not a body functor always run it
                else
                {
                    (*functor)(*m_constraints);
                }
            }
        }

        // Partition constraints for parallel computation
        if (m_config->m_doPartitioning)
        {
            m_constraints->partitionConstraints(static_cast<int>(m_partitionThreshold));
        }
        else
        {
            m_constraints->clearPartitions();
        }
    }

    // Setup the default pbd solver if none exists
    if (m_pbdSolver == nullptr)
    {
        m_pbdSolver = std::make_shared<PbdSolver>();
    }

    return true;
}

void
PbdModel::initGraphEdges(std::shared_ptr<TaskNode> source, std::shared_ptr<TaskNode> sink)
{
    // Setup graph connectivity
    m_taskGraph->addEdge(source, m_integrationPositionNode);
    m_taskGraph->addEdge(m_integrationPositionNode, m_solveConstraintsNode);
    m_taskGraph->addEdge(m_solveConstraintsNode, m_updateVelocityNode);
    m_taskGraph->addEdge(m_updateVelocityNode, sink);
}

void
PbdModel::addConstraints(std::shared_ptr<std::unordered_set<size_t>> vertices)
{
    for (const auto& functorVec : m_config->m_functors)
    {
        for (const auto& functor : functorVec.second)
        {
            functor->addConstraints(*m_constraints, vertices);
        }
    }
}

void
PbdModel::integratePosition()
{
    for (const auto& body : getCurrentState()->m_bodies)
    {
        integratePosition(*body);
    }
}

void
PbdModel::integratePosition(const PbdBody& body)
{
    VecDataArray<double, 3>& pos       = *body.vertices;
    VecDataArray<double, 3>& prevPos   = *body.prevVertices;
    VecDataArray<double, 3>& vel       = *body.velocities;
    VecDataArray<double, 3>& accn      = *body.accelerations;
    const DataArray<double>& invMasses = *body.invMasses;

    // Check all the arrays are the same
    const int numVerts = pos.size();
    CHECK(numVerts == prevPos.size()) << "PbdModel data corrupt";
    CHECK(numVerts == vel.size()) << "PbdModel data corrupt";
    CHECK(numVerts == accn.size()) << "PbdModel data corrupt";
    CHECK(numVerts == invMasses.size()) << "PbdModel data corrupt";

    ParallelUtils::parallelFor(numVerts,
        [&](const int i)
        {
            if (std::abs(invMasses[i]) > 0.0)
            {
                vel[i]    += (accn[i] + m_config->m_gravity) * m_config->m_dt;
                accn[i]    = Vec3d::Zero();
                prevPos[i] = pos[i];
                pos[i]    += (1.0 - m_config->m_viscousDampingCoeff) * vel[i] * m_config->m_dt;
            }
        }, numVerts > 50); // Only run parallel when more than 50 pts
}

void
PbdModel::updateVelocity()
{
    for (const auto& body : getCurrentState()->m_bodies)
    {
        updateVelocity(*body);
    }
}

void
PbdModel::updateVelocity(const PbdBody& body)
{
    if (m_config->m_dt > 0.0)
    {
        const VecDataArray<double, 3>& pos       = *body.vertices;
        const VecDataArray<double, 3>& prevPos   = *body.prevVertices;
        VecDataArray<double, 3>&       vel       = *body.velocities;
        const DataArray<double>&       invMasses = *body.invMasses;

        // Check all the arrays are the same
        const int numVerts = pos.size();
        CHECK(numVerts == prevPos.size()) << "PbdModel data corrupt";
        CHECK(numVerts == vel.size()) << "PbdModel data corrupt";
        CHECK(numVerts == invMasses.size()) << "PbdModel data corrupt";

        const double invDt = 1.0 / m_config->m_dt;
        ParallelUtils::parallelFor(numVerts,
            [&](const int i)
            {
                if (std::abs(invMasses[i]) > 0.0)
                {
                    vel[i] = (pos[i] - prevPos[i]) * invDt;
                }
            }, numVerts > 50);
    }
}

void
PbdModel::solveConstraints()
{
    m_pbdSolver->setPbdBodies(getCurrentState()->m_bodies);
    m_pbdSolver->setConstraints(getConstraints());
    m_pbdSolver->setTimeStep(m_config->m_dt);
    m_pbdSolver->setIterations(m_config->m_iterations);
    m_pbdSolver->setSolverType(m_config->m_solverType);
    m_pbdSolver->solve();
}
} // namespace imstk