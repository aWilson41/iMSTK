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

#include "imstkAnalyticalGeometry.h"
#include "imstkCollisionData.h"
#include "imstkPbdPointPointConstraint.h"
#include "imstkPbdPointTriangleConstraint.h"
#include "imstkPbdModel.h"
#include "imstkPbdObject.h"
#include "imstkPBDPickingCH.h"
#include "imstkPbdSolver.h"
#include "imstkPointSet.h"
#include "imstkSurfaceMesh.h"
#include "imstkGeometryMap.h"

namespace imstk
{
PBDPickingCH::PBDPickingCH(const CollisionHandling::Side&       side,
                           const std::shared_ptr<CollisionData> colData,
                           std::shared_ptr<PbdObject>           pbdObj,
                           std::shared_ptr<CollidingObject>     pickObj) :
    CollisionHandling(Type::PBDPicking, side, colData),
    m_pbdObj(pbdObj),
    m_pickObj(pickObj),
    m_pbdCollisionSolver(std::make_shared<PbdCollisionSolver>())
{
    m_isPicking = false;
    m_pickedPtIdxOffset.clear();
}

PBDPickingCH::~PBDPickingCH()
{
    for (const auto ptr : m_PPConstraintPool)
    {
        delete ptr;
    }
}

void
PBDPickingCH::processCollisionData()
{
    CHECK(m_pbdObj != nullptr && m_pickObj != nullptr)
        << "PBDPickingCH::handleCollision error: "
        << "no picking collision handling available the object";

    if (m_isPicking)
    {
        this->updatePickConstraints();
    }
    else
    {
        this->generatePBDConstraints();
        if (m_PBDConstraints.size() == 0)
        {
            return;
        }

        m_pbdCollisionSolver->addCollisionConstraints(&m_PBDConstraints);
    }
}

void
PBDPickingCH::updatePickConstraints()
{
    if (m_pickedPtIdxOffset.size() == 0)
    {
        this->removePickConstraints();
        return;
    }

    std::shared_ptr<PbdModel>                model      = m_pbdObj->getPbdModel();
    std::shared_ptr<AnalyticalGeometry>      pickGeom   = std::dynamic_pointer_cast<AnalyticalGeometry>(m_pickObj->getCollidingGeometry());
    std::shared_ptr<VecDataArray<double, 3>> vertices   = model->getCurrentState()->getPositions();
    VecDataArray<double, 3>&                 vertexData = *vertices;
    for (auto iter = m_pickedPtIdxOffset.begin(); iter != m_pickedPtIdxOffset.end(); iter++)
    {
        auto rot = pickGeom->getRotation();
        vertexData[iter->first] = pickGeom->getPosition() + rot * iter->second;
    }
}

void
PBDPickingCH::addPickConstraints(std::shared_ptr<PbdObject> pbdObj, std::shared_ptr<CollidingObject> pickObj)
{
    if (m_colData->PColData.isEmpty())
    {
        return;
    }

    CHECK(pbdObj != nullptr && pickObj != nullptr)
        << "PBDPickingCH:addPickConstraints error: "
        << "no pdb object or colliding object.";

    std::shared_ptr<PbdModel>           model    = pbdObj->getPbdModel();
    std::shared_ptr<AnalyticalGeometry> pickGeom = std::dynamic_pointer_cast<AnalyticalGeometry>(pickObj->getCollidingGeometry());
    CHECK(pickGeom != nullptr) << "Colliding geometry is analytical geometry ";

    std::shared_ptr<VecDataArray<double, 3>> vertices   = model->getCurrentState()->getPositions();
    VecDataArray<double, 3>&                 vertexData = *vertices;

    ParallelUtils::SpinLock lock;
    ParallelUtils::parallelFor(m_colData->PColData.getSize(),
        [&](const size_t idx) {
            const auto& cd = m_colData->PColData[idx];
            if (m_pickedPtIdxOffset.find(cd.nodeIdx) == m_pickedPtIdxOffset.end() && model->getInvMass(cd.nodeIdx) != 0.0)
            {
                auto pv  = cd.penetrationVector;
                auto rot = pickGeom->getRotation().transpose();
                auto relativePos = rot * (vertexData[cd.nodeIdx] - pv - pickGeom->getPosition());

                lock.lock();
                m_pickedPtIdxOffset[cd.nodeIdx] = relativePos;
                model->setFixedPoint(cd.nodeIdx);
                vertexData[cd.nodeIdx] = pickGeom->getPosition() + rot.transpose() * relativePos;
                lock.unlock();
            }
    });
}

void
PBDPickingCH::removePickConstraints()
{
    std::shared_ptr<PbdModel> model = m_pbdObj->getPbdModel();
    m_isPicking = false;
    for (auto iter = m_pickedPtIdxOffset.begin(); iter != m_pickedPtIdxOffset.end(); ++iter)
    {
        model->setPointUnfixed(iter->first);
    }
    m_pickedPtIdxOffset.clear();
}

void
PBDPickingCH::activatePickConstraints()
{
    if (!m_colData->PColData.isEmpty())
    {
        this->addPickConstraints(m_pbdObj, m_pickObj);
        m_isPicking = true;
    }
}

void
PBDPickingCH::generatePBDConstraints()
{
    const size_t totalPPConstraints = m_colData->PDColData.getSize() + m_colData->PColData.getSize();
    const size_t totalVTConstraints = m_colData->TFVColData.getSize();
    m_PPConstraintPoolSize = 0;
    m_VTConstraintPoolSize = 0;

    if (m_PPConstraintPool.size() < totalPPConstraints)
    {
        for (size_t i = m_PPConstraintPool.size(); i < totalPPConstraints; i++)
        {
            m_PPConstraintPool.push_back(new PbdPointPointConstraint);
        }
    }
    if (m_VTConstraintPool.size() < totalVTConstraints)
    {
        for (size_t i = m_VTConstraintPool.size(); i < totalVTConstraints; i++)
        {
            m_VTConstraintPool.push_back(new PbdPointTriangleConstraint);
        }
    }

    // Fixed vertex constraints
    // The constraints use pointer values
    m_fixedVertexPairs.resize(
        m_colData->PDColData.getSize() +
        m_colData->PColData.getSize() +
        m_colData->TFVColData.getSize());
    int fixedPairIter = 0;

    // Currently handles both PColData and PDColData
    std::shared_ptr<PbdCollisionConstraintConfig> config   = m_pbdObj->getPbdModel()->getConfig()->m_collisionParams;
    std::shared_ptr<PointSet>                     pointSet = std::dynamic_pointer_cast<PointSet>(m_pbdObj->getPhysicsGeometry());

    std::shared_ptr<VecDataArray<double, 3>> verticesPtr = pointSet->getVertexPositions();
    VecDataArray<double, 3>&                 vertices    = *verticesPtr;

    std::shared_ptr<DataArray<double>> invMassesPtr = m_pbdObj->getPbdModel()->getInvMasses();
    DataArray<double>&                 invMasses    = *invMassesPtr;

    std::shared_ptr<GeometryMap> map = m_pbdObj->getPhysicsToCollidingMap();

    // Process Position-Direction Data
    if (m_colData->PDColData.getSize() > 0)
    {
        for (int i = 0; i < m_colData->PDColData.getSize(); i++)
        {
            const PositionDirectionCollisionDataElement& colData = m_colData->PDColData[i];

            const Vec3d penVec = colData.dirAtoB * colData.penetrationDepth;
            m_fixedVertexPairs[fixedPairIter++] = std::pair<Vec3d, double>(colData.posB - penVec, 0.0);

            int v1 = colData.nodeIdx;
            if (map && map->getType() == GeometryMap::Type::OneToOne)
            {
                // Get the ids of the simulated mesh
                v1 = map->getMapIdx(v1);
            }

            m_PPConstraintPool[i]->initConstraint(
                PbdCollisionConstraint::Side::B,
                &m_fixedVertexPairs[fixedPairIter - 1].first, &m_fixedVertexPairs[fixedPairIter - 1].second,
                &vertices[v1], &invMasses[v1],
                nullptr, config);
            m_PPConstraintPoolSize++;
        }
    }

    // Process Penetration Data
    if (m_colData->PColData.getSize() > 0)
    {
        const size_t shift = m_PPConstraintPoolSize;
        for (int i = 0; i < m_colData->PColData.getSize(); i++)
        {
            const PenetrationCollisionDataElement& colData = m_colData->PColData[i];

            int v1 = colData.nodeIdx;
            if (map && map->getType() == GeometryMap::Type::OneToOne)
            {
                // Get the ids of the simulated mesh
                v1 = map->getMapIdx(v1);
            }

            m_fixedVertexPairs[fixedPairIter++] = std::pair<Vec3d, double>(vertices[v1] - colData.penetrationVector, 0.0);

            m_PPConstraintPool[i]->initConstraint(
                PbdCollisionConstraint::Side::B,
                &m_fixedVertexPairs[fixedPairIter - 1].first, &m_fixedVertexPairs[fixedPairIter - 1].second,
                &vertices[v1], &invMasses[v1],
                nullptr, config);
            m_PPConstraintPoolSize++;
        }
    }

    // Process Triangle-FixedVertex Data
    if (m_colData->TFVColData.getSize() > 0)
    {
        // With VT constraints assume obj1 has SurfaceMesh
        std::shared_ptr<VecDataArray<int, 3>> indicesPtr = std::dynamic_pointer_cast<SurfaceMesh>(m_pbdObj->getCollidingGeometry())->getTriangleIndices();
        const VecDataArray<int, 3>& indices = *indicesPtr;

        const size_t shift = m_VTConstraintPoolSize;
        for (int i = 0; i < m_colData->TFVColData.getSize(); i++)
        {
            const TriangleFixedVertexCollisionDataElement& colData = m_colData->TFVColData[i];
            const Vec3i& triVerts = indices[colData.triIdx];

            size_t v1, v2, v3;
            if (map && map->getType() == GeometryMap::Type::OneToOne)
            {
                v1 = map->getMapIdx(triVerts[0]);
                v2 = map->getMapIdx(triVerts[1]);
                v3 = map->getMapIdx(triVerts[2]);
            }
            else
            {
                v1 = triVerts[0];
                v2 = triVerts[1];
                v3 = triVerts[2];
            }

            m_fixedVertexPairs[fixedPairIter++] = std::pair<Vec3d, double>(colData.vertexPt, 0.0);

           m_VTConstraintPool[i + shift]->initConstraint(
                PbdCollisionConstraint::Side::B,
                &m_fixedVertexPairs[fixedPairIter - 1].first, &m_fixedVertexPairs[fixedPairIter - 1].second,
                &vertices[v1], &invMasses[v1], &vertices[v2], &invMasses[v2], &vertices[v3], &invMasses[v3],
                nullptr, config, colData.dir, colData.closestDistance);
            m_VTConstraintPoolSize++;
        }
    }

    // Copy constraints
    m_PBDConstraints.resize(0);
    m_PBDConstraints.reserve(totalPPConstraints + totalVTConstraints);

    for (size_t i = 0; i < totalPPConstraints; i++)
    {
        m_PBDConstraints.push_back(static_cast<PbdCollisionConstraint*>(m_PPConstraintPool[i]));
    }
    for (size_t i = 0; i < totalVTConstraints; i++)
    {
        m_PBDConstraints.push_back(static_cast<PbdCollisionConstraint*>(m_VTConstraintPool[i]));
    }
}
}
