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
#include "imstkTaskGraph.h"
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
    // Add a dummy pbdy body for virtual particle addition
    addBody();

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
    m_collisionSolveConstraintsNode = m_taskGraph->addFunction("PbdModel_SolveCollisionConstraints",
        [&]() { solveCollisionConstraints(); });
    m_updateVelocityNode = m_taskGraph->addFunction("PbdModel_UpdateVelocity",
        [&]() { updateVelocity(); });
}

void
PbdModel::configure(std::shared_ptr<PbdModelConfig> config)
{
    m_config = config;
}

std::shared_ptr<PbdBody>
PbdModel::addBody()
{
    m_state.m_bodies.push_back(std::make_shared<PbdBody>(m_iterKey));
    m_modified = true;
    m_iterKey++;
    return m_state.m_bodies.back();
}

void
PbdModel::removeBody(std::shared_ptr<PbdBody> body)
{
    auto iter = std::find(m_state.m_bodies.begin(), m_state.m_bodies.end(), body);
    CHECK(iter != m_state.m_bodies.end()) << "removeBody called but could not find PbdyBody in PbdState";
    m_state.m_bodies.erase(iter);
    m_modified = true;
}

PbdParticleId
PbdModel::addVirtualParticle(
    const Vec3d& pos, const Quatd& orientation,
    const double mass, const Mat3d inertia,
    const Vec3d velocity, const Vec3d angularVelocity,
    const Vec3d acceleration, const Vec3d angularAccel)
{
    m_state.m_bodies[0]->vertices->push_back(pos);
    m_state.m_bodies[0]->orientations->push_back(orientation);
    m_state.m_bodies[0]->velocities->push_back(velocity);
    m_state.m_bodies[0]->angularVelocities->push_back(angularVelocity);
    m_state.m_bodies[0]->accelerations->push_back(acceleration);
    m_state.m_bodies[0]->angularAccel->push_back(angularAccel);
    m_state.m_bodies[0]->masses->push_back(mass);
    m_state.m_bodies[0]->inertias->push_back(inertia);
    return { 0, m_state.m_bodies[0]->vertices->size() - 1 };
}

PbdParticleId
PbdModel::addVirtualParticle(
    const Vec3d& pos, const double mass,
    const Vec3d velocity, const Vec3d acceleration)
{
    return addVirtualParticle(pos, Quatd::Identity(),
        mass, Mat3d::Identity(),
        velocity, Vec3d::Zero(),
        acceleration, Vec3d::Zero());
}

void
PbdModel::clearVirtualParticles()
{
    CHECK(m_state.m_bodies.size() != 0 || m_state.m_bodies[0] != nullptr) << "Missing virtual/dummy body";
    resizeBodyParticles(*m_state.m_bodies[0], 0);
}

void
PbdModel::resizeBodyParticles(PbdBody& body, const int particleCount)
{
    body.prevVertices->resize(particleCount);
    body.vertices->resize(particleCount);
    body.velocities->resize(particleCount);
    body.accelerations->resize(particleCount);
    body.masses->resize(particleCount);
    body.invMasses->resize(particleCount);
    if (body.getOriented())
    {
        body.prevOrientations->resize(particleCount);
        body.orientations->resize(particleCount);
        body.angularVelocities->resize(particleCount);
        body.angularAccel->resize(particleCount);
        body.inertias->resize(particleCount);
        body.invInertias->resize(particleCount);
    }
}

bool
PbdModel::initialize()
{
    // Clear the dummy body for virtual particles, init with space for 10 particles
    m_state.m_bodies[0] = std::make_shared<PbdBody>(0);
    m_state.m_bodies[0]->bodyType          = PbdBody::Type::DEFORMABLE_ORIENTED;
    m_state.m_bodies[0]->prevVertices      = std::make_shared<VecDataArray<double, 3>>();
    m_state.m_bodies[0]->vertices          = std::make_shared<VecDataArray<double, 3>>();
    m_state.m_bodies[0]->prevOrientations  = std::make_shared<StdVectorOfQuatd>();
    m_state.m_bodies[0]->orientations      = std::make_shared<StdVectorOfQuatd>();
    m_state.m_bodies[0]->velocities        = std::make_shared<VecDataArray<double, 3>>();
    m_state.m_bodies[0]->angularVelocities = std::make_shared<VecDataArray<double, 3>>();
    m_state.m_bodies[0]->accelerations     = std::make_shared<VecDataArray<double, 3>>();
    m_state.m_bodies[0]->angularAccel      = std::make_shared<VecDataArray<double, 3>>();
    m_state.m_bodies[0]->masses      = std::make_shared<DataArray<double>>();
    m_state.m_bodies[0]->invMasses   = std::make_shared<DataArray<double>>();
    m_state.m_bodies[0]->inertias    = std::make_shared<StdVectorOfMat3d>();
    m_state.m_bodies[0]->invInertias = std::make_shared<StdVectorOfMat3d>();

    // Initialize constraints
    {
        m_constraints = std::make_shared<PbdConstraintContainer>();

        m_config->computeElasticConstants();

        // Run all the functors to generate the constraints
        for (const auto& functorVec : m_config->m_functors)
        {
            for (const auto& functor : functorVec.second)
            {
                (*functor)(*m_constraints);
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
    if (m_pbdCollisionSolver == nullptr)
    {
        m_pbdCollisionSolver = std::make_shared<PbdSolver>();
    }

    return true;
}

void
PbdModel::initGraphEdges(std::shared_ptr<TaskNode> source, std::shared_ptr<TaskNode> sink)
{
    // Setup graph connectivity
    m_taskGraph->addEdge(source, m_integrationPositionNode);
    m_taskGraph->addEdge(m_integrationPositionNode, m_solveConstraintsNode);
    m_taskGraph->addEdge(m_solveConstraintsNode, m_collisionSolveConstraintsNode);
    m_taskGraph->addEdge(m_collisionSolveConstraintsNode, m_updateVelocityNode);
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
    // resize 0 virtual particles (avoids reallocation)
    clearVirtualParticles();

    //for (const auto& body : m_state.m_bodies)
    for (auto bodyIter = std::next(m_state.m_bodies.begin());
         bodyIter != m_state.m_bodies.end(); bodyIter++)
    {
        integratePosition(**bodyIter);
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
    const Vec3d              bodyExternalForce  = body.externalForce;
    const Vec3d              bodyExternalTorque = body.externalTorque;

    // Check all the arrays are the same
    const int numParticles = pos.size();
    CHECK(numParticles == prevPos.size()) << "PbdModel data corrupt";
    CHECK(numParticles == vel.size()) << "PbdModel data corrupt";
    CHECK(numParticles == accn.size()) << "PbdModel data corrupt";
    CHECK(numParticles == invMasses.size()) << "PbdModel data corrupt";

    const double linearVelocityDamp = (1.0 - m_config->m_linearDampingCoeff);
    ParallelUtils::parallelFor(numParticles,
        [&](const int i)
        {
            if (std::abs(invMasses[i]) > 0.0)
            {
                vel[i] += (accn[i] + m_config->m_gravity +
                           bodyExternalForce * invMasses[i]) * m_config->m_dt;
                accn[i]    = Vec3d::Zero();
                prevPos[i] = pos[i];
                pos[i]    += linearVelocityDamp * vel[i] * m_config->m_dt;
            }
        }, numParticles > 50); // Only run parallel when more than 50 pts

    // If using oriented particles update those too
    if (body.getOriented())
    {
        StdVectorOfQuatd&        orientations      = *body.orientations;
        StdVectorOfQuatd&        prevOrientations  = *body.prevOrientations;
        VecDataArray<double, 3>& angularVelocities = *body.angularVelocities;
        VecDataArray<double, 3>& angularAccel      = *body.angularAccel;
        StdVectorOfMat3d&        invInertias       = *body.invInertias;

        // Check all the arrays are the same
        CHECK(numParticles == orientations.size()) << "PbdModel data corrupt";
        CHECK(numParticles == prevOrientations.size()) << "PbdModel data corrupt";
        CHECK(numParticles == angularVelocities.size()) << "PbdModel data corrupt";
        CHECK(numParticles == angularAccel.size()) << "PbdModel data corrupt";
        CHECK(numParticles == invInertias.size()) << "PbdModel data corrupt";

        const double angularVelocityDamp = (1.0 - m_config->m_angularDampingCoeff);
        ParallelUtils::parallelFor(numParticles,
            [&](const int i)
            {
                if (!invInertias[i].isZero())
                {
                    angularVelocities[i] += (angularAccel[i] +
                                             invInertias[i] * bodyExternalTorque) * m_config->m_dt;
                    angularAccel[i]     = Vec3d::Zero();
                    prevOrientations[i] = orientations[i];

                    const Quatd q = Quatd(0.0,
                        angularVelocities[i][0] * angularVelocityDamp * 0.5,
                        angularVelocities[i][1] * angularVelocityDamp * 0.5,
                        angularVelocities[i][2] * angularVelocityDamp * 0.5) * orientations[i];
                    orientations[i].x() += q.x() * m_config->m_dt;
                    orientations[i].y() += q.y() * m_config->m_dt;
                    orientations[i].z() += q.z() * m_config->m_dt;
                    orientations[i].w() += q.w() * m_config->m_dt;
                    orientations[i].normalize();
                }
            }, numParticles > 50); // Only run parallel when more than 50 pts
    }
}

void
PbdModel::updateVelocity()
{
    for (auto bodyIter = std::next(m_state.m_bodies.begin());
         bodyIter != m_state.m_bodies.end(); bodyIter++)
    {
        updateVelocity(**bodyIter);
    }

    // Correctly velocities for friction and restitution
    // Unfortunately the constraint would be clear after a solve
    for (const auto& colConstraintList : m_pbdCollisionSolver->getConstraintLists())
    {
        for (auto& colConstraint : *colConstraintList)
        {
            colConstraint->correctVelocity(m_state);
        }
    }
    m_pbdCollisionSolver->clearConstraintLists();
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
        const int numParticles = pos.size();
        CHECK(numParticles == prevPos.size()) << "PbdModel data corrupt";
        CHECK(numParticles == vel.size()) << "PbdModel data corrupt";
        CHECK(numParticles == invMasses.size()) << "PbdModel data corrupt";

        const double invDt = 1.0 / m_config->m_dt;
        ParallelUtils::parallelFor(numParticles,
            [&](const int i)
            {
                if (std::abs(invMasses[i]) > 0.0)
                {
                    vel[i] = (pos[i] - prevPos[i]) * invDt;
                }
            }, numParticles > 50);

        if (body.getOriented())
        {
            const StdVectorOfQuatd&  orientations      = *body.orientations;
            const StdVectorOfQuatd&  prevOrientations  = *body.prevOrientations;
            VecDataArray<double, 3>& angularVelocities = *body.angularVelocities;
            const StdVectorOfMat3d&  invInertias       = *body.invInertias;

            // Check all the arrays are the same
            CHECK(numParticles == orientations.size()) << "PbdModel data corrupt";
            CHECK(numParticles == prevOrientations.size()) << "PbdModel data corrupt";
            CHECK(numParticles == angularVelocities.size()) << "PbdModel data corrupt";
            CHECK(numParticles == invInertias.size()) << "PbdModel data corrupt";

            ParallelUtils::parallelFor(numParticles,
                [&](const int i)
                {
                    if (!invInertias[i].isZero())
                    {
                        const Quatd dq = orientations[i] * prevOrientations[i].inverse();
                        //const Quatd dq = prevOrientations[i] * orientations[i].inverse();
                        const Vec3d angularVel = 2.0 * Vec3d(dq.x(), dq.y(), dq.z()) * invDt;
                        angularVelocities[i]   = dq.w() >= 0.0 ? angularVel : -angularVel;
                    }
                }, numParticles > 50);
        }
    }
}

void
PbdModel::solveConstraints()
{
    m_pbdSolver->setPbdBodies(&m_state);
    m_pbdSolver->setConstraints(getConstraints());
    m_pbdSolver->setTimeStep(m_config->m_dt);
    m_pbdSolver->setIterations(m_config->m_iterations);
    m_pbdSolver->setSolverType(m_config->m_solverType);
    m_pbdSolver->solve();
}

void
PbdModel::solveCollisionConstraints()
{
    m_pbdCollisionSolver->setPbdBodies(&m_state);
    m_pbdCollisionSolver->setTimeStep(m_config->m_dt);
    m_pbdCollisionSolver->setIterations(m_config->m_collisionIterations);
    m_pbdCollisionSolver->setSolverType(m_config->m_solverType);
    m_pbdCollisionSolver->solve();
}
} // namespace imstk