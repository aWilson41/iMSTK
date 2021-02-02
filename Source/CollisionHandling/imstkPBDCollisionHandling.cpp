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

#include "imstkPBDCollisionHandling.h"
#include "imstkCollisionData.h"
#include "imstkGeometryMap.h"
#include "imstkPbdEdgeEdgeConstraint.h"
#include "imstkPbdModel.h"
#include "imstkPbdObject.h"
#include "imstkPbdPointPointConstraint.h"
#include "imstkPbdPointTriangleConstraint.h"
#include "imstkPbdSolver.h"
#include "imstkSurfaceMesh.h"

namespace imstk
{
PBDCollisionHandling::PBDCollisionHandling(const std::shared_ptr<CollisionData> colData,
                                           std::shared_ptr<PbdObject>           object1,
                                           std::shared_ptr<PbdObject>           object2) :
    CollisionHandling(Type::PBD, Side::AB, colData), // Handle both sides
    m_obj1(object1),
    m_obj2(object2),
    m_pbdCollisionSolver(std::make_shared<PbdCollisionSolver>())
{
}

PBDCollisionHandling::PBDCollisionHandling(const std::shared_ptr<CollisionData> colData,
                                           std::shared_ptr<PbdObject>           object1,
                                           std::shared_ptr<CollidingObject>     object2) :
    CollisionHandling(Type::PBD, Side::A, colData), // Only handle side A
    m_obj1(object1),
    m_obj2(object2),
    m_pbdCollisionSolver(std::make_shared<PbdCollisionSolver>())
{
}

PBDCollisionHandling::~PBDCollisionHandling()
{
    for (const auto ptr: m_EEConstraintPool)
    {
        delete ptr;
    }
    for (const auto ptr: m_VTConstraintPool)
    {
        delete ptr;
    }
    for (const auto ptr : m_PPConstraintPool)
    {
        delete ptr;
    }
}

void
PBDCollisionHandling::processCollisionData()
{
    generateConstraints();

    if (m_PBDConstraints.size() != 0)
    {
        m_pbdCollisionSolver->addCollisionConstraints(&m_PBDConstraints);
    }
}

void
PBDCollisionHandling::generateConstraints()
{
    // Get and assume obj1/sideA geometry is at least a pointset with pbd
    std::shared_ptr<PointSet>                obj1PointSet  = std::dynamic_pointer_cast<PointSet>(m_obj1->getPhysicsGeometry());
    std::shared_ptr<VecDataArray<double, 3>> vertices1Ptr  = obj1PointSet->getVertexPositions();
    VecDataArray<double, 3>&                 vertices1     = *vertices1Ptr;
    std::shared_ptr<DataArray<double>>       invMasses1Ptr = std::dynamic_pointer_cast<DataArray<double>>(obj1PointSet->getVertexAttribute("InvMass"));
    DataArray<double>&                       invMasses1    = *invMasses1Ptr;

    // For obj2/sideB we don't neccesarily have a pbd object and model
    std::shared_ptr<PointSet> obj2PointSet = std::dynamic_pointer_cast<PointSet>(m_obj2->getCollidingGeometry());
    if (auto obj2Pbd = std::dynamic_pointer_cast<PbdObject>(m_obj2))
    {
        obj2PointSet = std::dynamic_pointer_cast<PointSet>(obj2Pbd->getPhysicsGeometry());
    }

    // Get the configs, nullptr for config2 if obj2 is not a PbdObject
    std::shared_ptr<PbdCollisionConstraintConfig> config1 = m_obj1->getPbdModel()->getConfig()->m_collisionParams;
    std::shared_ptr<GeometryMap>                  map1    = m_obj1->getPhysicsToCollidingMap();

    std::shared_ptr<PbdCollisionConstraintConfig> config2 = nullptr;
    std::shared_ptr<GeometryMap>                  map2    = nullptr;

    if (auto obj2Pbd = std::dynamic_pointer_cast<PbdObject>(m_obj2))
    {
        config2 = obj2Pbd->getPbdModel()->getConfig()->m_collisionParams;
    }
    if (auto obj2Dynamic = std::dynamic_pointer_cast<DynamicObject>(m_obj2))
    {
        map2 = obj2Dynamic->getPhysicsToCollidingMap();
    }

    const PbdCollisionConstraint::Side pbdSide = (m_side == Side::A) ? PbdCollisionConstraint::Side::A : PbdCollisionConstraint::Side::AB;

    // Upsize pools if needed to maximum size
    const size_t totalEEConstraints = m_colData->EEColData.getSize() + m_colData->EFEColData.getSize();
    const size_t totalPPConstraints = m_colData->PDColData.getSize() + m_colData->PColData.getSize();
    const size_t totalVTConstraints = m_colData->VTColData.getSize() + m_colData->TFVColData.getSize() + m_colData->TVColData.getSize();

    // Keep track of the actual size
    m_EEConstraintPoolSize = m_VTConstraintPoolSize = m_PPConstraintPoolSize = 0;

    if (m_EEConstraintPool.size() < totalEEConstraints)
    {
        // Allocate the new ones as size has increased
        for (size_t i = m_EEConstraintPool.size(); i < totalEEConstraints; i++)
        {
            m_EEConstraintPool.push_back(new PbdEdgeEdgeConstraint);
        }
    }
    if (m_VTConstraintPool.size() < totalVTConstraints)
    {
        for (size_t i = m_VTConstraintPool.size(); i < totalVTConstraints; i++)
        {
            m_VTConstraintPool.push_back(new PbdPointTriangleConstraint);
        }
    }
    if (m_PPConstraintPool.size() < totalPPConstraints)
    {
        for (size_t i = m_PPConstraintPool.size(); i < totalPPConstraints; i++)
        {
            m_PPConstraintPool.push_back(new PbdPointPointConstraint);
        }
    }

    // Fixed vertex constraints
    // The constraints use pointer values
    m_fixedVertexPairs.resize(
        m_colData->EFEColData.getSize() * 2 +
        m_colData->TFVColData.getSize() +
        m_colData->PDColData.getSize() +
        m_colData->PColData.getSize());
    int fixedPairIter = 0;

    // Process Edge-Edge Data
    if (m_colData->EEColData.getSize() > 0)
    {
        // With EE constraints assume obj2 has PointSet
        std::shared_ptr<VecDataArray<double, 3>> vertices2Ptr  = std::dynamic_pointer_cast<PointSet>(m_obj2->getCollidingGeometry())->getVertexPositions();
        VecDataArray<double, 3>&                 vertices2     = *vertices2Ptr;
        std::shared_ptr<DataArray<double>>       invMasses2Ptr = std::dynamic_pointer_cast<DataArray<double>>(obj2PointSet->getVertexAttribute("InvMass"));
        DataArray<double>&                       invMasses2    = *invMasses2Ptr;

        for (int i = 0; i < m_colData->EEColData.getSize(); i++)
        {
            const auto& colData = m_colData->EEColData[i];

            size_t edgeA1, edgeA2;
            if (map1 && map1->getType() == GeometryMap::Type::OneToOne)
            {
                // Get the ids of the simulated mesh
                edgeA1 = map1->getMapIdx(colData.edgeIdA.first);
                edgeA2 = map1->getMapIdx(colData.edgeIdA.second);
            }
            else // Collision was done using simulated mesh
            {
                edgeA1 = colData.edgeIdA.first;
                edgeA2 = colData.edgeIdA.second;
            }

            size_t edgeB1, edgeB2;
            if (map2 && map2->getType() == GeometryMap::Type::OneToOne)
            {
                edgeB1 = map2->getMapIdx(colData.edgeIdB.first);
                edgeB2 = map2->getMapIdx(colData.edgeIdB.second);
            }
            else
            {
                edgeB1 = colData.edgeIdB.first;
                edgeB2 = colData.edgeIdB.second;
            }

            m_EEConstraintPool[i]->initConstraint(
                pbdSide,
                &vertices1[edgeA1], &invMasses1[edgeA1],
                &vertices1[edgeA2], &invMasses1[edgeA2],
                &vertices2[edgeB1], &invMasses2[edgeB1],
                &vertices2[edgeB2], &invMasses2[edgeB2],
                config1, config2);
            m_EEConstraintPoolSize++;
        }
    }

    // Process Vertex-Triangle Data
    if (m_colData->VTColData.getSize() > 0)
    {
        // With VT constraints assume obj2 has SurfaceMesh
        std::shared_ptr<VecDataArray<double, 3>> vertices2Ptr  = std::dynamic_pointer_cast<PointSet>(m_obj2->getCollidingGeometry())->getVertexPositions();
        VecDataArray<double, 3>&                 vertices2     = *vertices2Ptr;
        std::shared_ptr<DataArray<double>>       invMasses2Ptr = std::dynamic_pointer_cast<DataArray<double>>(obj2PointSet->getVertexAttribute("InvMass"));
        DataArray<double>&                       invMasses2    = *invMasses2Ptr;

        // If vertex-triangle data, assume SurfaceMesh
        std::shared_ptr<VecDataArray<int, 3>> indicesPtr = std::static_pointer_cast<SurfaceMesh>(m_obj2->getCollidingGeometry())->getTriangleIndices();
        const VecDataArray<int, 3>&           indices    = *indicesPtr;
        for (int i = 0; i < m_colData->VTColData.getSize(); i++)
        {
            const VertexTriangleCollisionDataElement& colData  = m_colData->VTColData[i];
            const Vec3i&                              triVerts = indices[colData.triIdx];

            size_t v1, v2, v3;
            if (map2 && map2->getType() == GeometryMap::Type::OneToOne)
            {
                v1 = map2->getMapIdx(triVerts[0]);
                v2 = map2->getMapIdx(triVerts[1]);
                v3 = map2->getMapIdx(triVerts[2]);
            }
            else
            {
                v1 = triVerts[0];
                v2 = triVerts[1];
                v3 = triVerts[2];
            }

           /* m_VTConstraintPool[i]->initConstraint(
                pbdSide,
                &vertices1[colData.vertexIdx], &invMasses1[colData.vertexIdx],
                &vertices2[v1], &invMasses2[v1], &vertices2[v2], &invMasses2[v2], &vertices2[v3], &invMasses2[v3],
                config1, config2);*/
            m_VTConstraintPoolSize++;
        }
    }

    // Process Edge-FixedEdge Data
    if (m_colData->EFEColData.getSize() > 0)
    {
        const size_t shift = m_EEConstraintPoolSize;
        for (int i = 0; i < m_colData->EFEColData.getSize(); i++)
        {
            const EdgeFixedEdgeCollisionDataElement& colData = m_colData->EFEColData[i];

            size_t edgeA1, edgeA2;
            if (map1 && map1->getType() == GeometryMap::Type::OneToOne)
            {
                // Get the ids of the simulated mesh
                edgeA1 = map1->getMapIdx(colData.edgeIdA.first);
                edgeA2 = map1->getMapIdx(colData.edgeIdA.second);
            }
            else // Collision was done using simulated mesh
            {
                edgeA1 = colData.edgeIdA.first;
                edgeA2 = colData.edgeIdA.second;
            }

            // Add fixed virtual points
            m_fixedVertexPairs[fixedPairIter++] = std::pair<Vec3d, double>(colData.edgePtsB.first, 0.0);
            m_fixedVertexPairs[fixedPairIter++] = std::pair<Vec3d, double>(colData.edgePtsB.second, 0.0);

            m_EEConstraintPool[i + shift]->initConstraint(
                PbdCollisionConstraint::Side::B,
                &vertices1[edgeA1], &invMasses1[edgeA1],
                &vertices1[edgeA2], &invMasses1[edgeA2],
                &m_fixedVertexPairs[fixedPairIter - 2].first, &m_fixedVertexPairs[fixedPairIter - 2].second,
                &m_fixedVertexPairs[fixedPairIter - 1].first, &m_fixedVertexPairs[fixedPairIter - 1].second,
                config2, config1);
            m_EEConstraintPoolSize++;
        }
    }

    // Process Triangle-FixedVertex Data
    if (m_colData->TFVColData.getSize() > 0)
    {
        // With VT constraints assume obj1 has SurfaceMesh
        std::shared_ptr<VecDataArray<int, 3>> indicesPtr = std::dynamic_pointer_cast<SurfaceMesh>(m_obj1->getCollidingGeometry())->getTriangleIndices();
        const VecDataArray<int, 3>&           indices    = *indicesPtr;

        const size_t shift = m_VTConstraintPoolSize;
        for (int i = 0; i < m_colData->TFVColData.getSize(); i++)
        {
            const TriangleFixedVertexCollisionDataElement& colData  = m_colData->TFVColData[i];
            const Vec3i&                                   triVerts = indices[colData.triIdx];

            size_t v1, v2, v3;
            if (map1 && map1->getType() == GeometryMap::Type::OneToOne)
            {
                v1 = map1->getMapIdx(triVerts[0]);
                v2 = map1->getMapIdx(triVerts[1]);
                v3 = map1->getMapIdx(triVerts[2]);
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
                &vertices1[v1], &invMasses1[v1], &vertices1[v2], &invMasses1[v2], &vertices1[v3], &invMasses1[v3],
                config2, config1, colData.dir, colData.closestDistance);
            m_VTConstraintPoolSize++;
        }
    }

    // Process Position-Direction Data
    if (m_colData->PDColData.getSize() > 0)
    {
        for (int i = 0; i < m_colData->PDColData.getSize(); i++)
        {
            const PositionDirectionCollisionDataElement& colData = m_colData->PDColData[i];

            const Vec3d penVec = colData.dirAtoB * colData.penetrationDepth;
            m_fixedVertexPairs[fixedPairIter++] = std::pair<Vec3d, double>(colData.posB - penVec, 0.0);

            int v1 = colData.nodeIdx;
            if (map1 && map1->getType() == GeometryMap::Type::OneToOne)
            {
                // Get the ids of the simulated mesh
                v1 = map1->getMapIdx(v1);
            }

            m_PPConstraintPool[i]->initConstraint(
                PbdCollisionConstraint::Side::B,
                &m_fixedVertexPairs[fixedPairIter - 1].first, &m_fixedVertexPairs[fixedPairIter - 1].second,
                &vertices1[v1], &invMasses1[v1],
                config2, config1);
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
            if (map1 && map1->getType() == GeometryMap::Type::OneToOne)
            {
                // Get the ids of the simulated mesh
                v1 = map1->getMapIdx(v1);
            }

            m_fixedVertexPairs[fixedPairIter++] = std::pair<Vec3d, double>(vertices1[v1] - colData.penetrationVector, 0.0);

            m_PPConstraintPool[i]->initConstraint(
                PbdCollisionConstraint::Side::B,
                &m_fixedVertexPairs[fixedPairIter - 1].first, &m_fixedVertexPairs[fixedPairIter - 1].second,
                &vertices1[v1], &invMasses1[v1],
                config2, config1);
            m_PPConstraintPoolSize++;
        }
    }

    // Copy constraints
    m_PBDConstraints.resize(0);
    m_PBDConstraints.reserve(totalEEConstraints + totalVTConstraints + totalPPConstraints);

    for (size_t i = 0; i < totalEEConstraints; i++)
    {
        m_PBDConstraints.push_back(static_cast<PbdCollisionConstraint*>(m_EEConstraintPool[i]));
    }
    for (size_t i = 0; i < totalVTConstraints; i++)
    {
        m_PBDConstraints.push_back(static_cast<PbdCollisionConstraint*>(m_VTConstraintPool[i]));
    }
    for (size_t i = 0; i < totalPPConstraints; i++)
    {
        m_PBDConstraints.push_back(static_cast<PbdCollisionConstraint*>(m_PPConstraintPool[i]));
    }
}
}
