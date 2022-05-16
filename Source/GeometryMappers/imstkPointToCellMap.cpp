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

#include "imstkPointToCellMap.h"
#include "imstkLineMesh.h"
#include "imstkLogger.h"
#include "imstkParallelUtils.h"
#include "imstkSurfaceMesh.h"
#include "imstkTetrahedralMesh.h"
#include "imstkVecDataArray.h"
#include "imstkCollisionUtils.h"

namespace imstk
{
//static std::unordered_map<int, std::vector<int>> triCaseToVertIdMap =
//{
//    { 0, { 0 } },
//    { 1, { 1 } },
//    { 2, { 2 } },
//    { 3, { 0, 1 } },
//    { 4, { 1, 2 } },
//    { 5, { 2, 3 } },
//    { 6, { 0, 1, 2 } }
//};
//static std::unordered_map<int, std::vector<int>> lineCaseToVertIdMap =
//{
//    { 0, { 0 } },
//    { 1, { 1 } },
//    { 2, { 0, 1 } }
//};

PointToCellMap::PointToCellMap()
{
    setRequiredInputType<PointSet>(0);
    setRequiredInputType<PointSet>(1);
}

PointToCellMap::PointToCellMap(
    std::shared_ptr<Geometry> parent,
    std::shared_ptr<Geometry> child)
{
    setRequiredInputType<PointSet>(0);
    setRequiredInputType<PointSet>(1);

    CHECK(std::dynamic_pointer_cast<LineMesh>(parent) != nullptr
        || std::dynamic_pointer_cast<SurfaceMesh>(parent) != nullptr
        || std::dynamic_pointer_cast<TetrahedralMesh>(parent) != nullptr) <<
        "PointToCellMap parent must be a LineMesh, SurfaceMesh, or TetrahedralMesh";

    setParentGeometry(parent);
    setChildGeometry(child);
}

void
PointToCellMap::compute()
{
    if (!areInputsValid())
    {
        LOG(WARNING) << "PointToCellMap failed to run, inputs not satisfied";
        return;
    }

    // Clear the map weights
    m_weights.clear();

    auto childMesh = std::dynamic_pointer_cast<PointSet>(getChildGeometry());

    // Setup the bone origin and transform for each cell/element
    if (auto tetMesh = std::dynamic_pointer_cast<TetrahedralMesh>(getParentGeometry()))
    {
        auto parentMesh = std::dynamic_pointer_cast<TetrahedralMesh>(getParentGeometry());
        //parentMesh->computeVertexToCells();

        LOG(WARNING) << "Not implemented yet";
    }
    else if (auto triMesh = std::dynamic_pointer_cast<SurfaceMesh>(getParentGeometry()))
    {
        auto                     parentMesh  = std::dynamic_pointer_cast<SurfaceMesh>(getParentGeometry());
        VecDataArray<double, 3>& parentVerts = *parentMesh->getVertexPositions();
        VecDataArray<int, 3>&    parentCells = *parentMesh->getIndices();

        for (int i = 0; i < parentMesh->getNumCells(); i++)
        {
            const Vec3i& cell = parentCells[i];
        }
    }
    else if (auto lineMesh = std::dynamic_pointer_cast<LineMesh>(getParentGeometry()))
    {
        auto                     parentMesh  = std::dynamic_pointer_cast<LineMesh>(getParentGeometry());
        VecDataArray<double, 3>& parentVerts = *parentMesh->getVertexPositions();
        VecDataArray<int, 2>&    parentCells = *parentMesh->getIndices();

        for (int i = 0; i < parentMesh->getNumCells(); i++)
        {
            const Vec2i& cell = parentCells[i];
            m_origin[i] = parentVerts[cell[0]];

            const Vec3d a    = parentVerts[cell[0]];
            const Vec3d b    = parentVerts[cell[1]];
            Vec3d       diff = (b - a).normalized();

            // Automatically determine a basis for each cell
            Mat3d transform = Mat3d::Identity();
            transform.col(0) = diff.cross(Vec3d(diff[1], -diff[0], -diff[2])).normalized();
            transform.col(1) = diff;
            transform.col(2) = transform.col(0).cross(transform.col(1)).normalized();
            m_transform[i]   = transform;
        }
    }
    else
    {
        LOG(FATAL) << "PointToCellMap failed to run, geometry case not implemented";
    }

    if (auto tetMesh = std::dynamic_pointer_cast<TetrahedralMesh>(getParentGeometry()))
    {
        auto parentMesh = std::dynamic_pointer_cast<TetrahedralMesh>(getParentGeometry());
        //parentMesh->computeVertexToCells();

        LOG(WARNING) << "Not implemented yet";
    }
    else if (auto triMesh = std::dynamic_pointer_cast<SurfaceMesh>(getParentGeometry()))
    {
        auto parentMesh = std::dynamic_pointer_cast<SurfaceMesh>(getParentGeometry());
        parentMesh->computeVertexToCells();

        // For every child vertex, find the nearest element transform with
        ParallelUtils::parallelFor(childMesh->getNumVertices(),
            [&](const int ptId)
            {
                const Vec3d& childPos = childMesh->getVertexPosition(ptId);
                computeTriangleWeights(ptId, childPos, m_weights);
                normalizeWeights(m_weights[ptId]);
            }, false);
    }
    else if (auto lineMesh = std::dynamic_pointer_cast<LineMesh>(getParentGeometry()))
    {
        auto parentMesh = std::dynamic_pointer_cast<LineMesh>(getParentGeometry());
        parentMesh->computeVertexToCells();

        // For every child vertex, find the nearest element transform with
        ParallelUtils::parallelFor(childMesh->getNumVertices(),
            [&](const int ptId)
            {
                const Vec3d& childPos = childMesh->getVertexPosition(ptId);
                computeSegmentWeights(ptId, childPos, m_weights);
                normalizeWeights(m_weights[ptId]);
            }, false);
    }
    else
    {
        LOG(FATAL) << "PointToCellMap failed to run, geometry case not implemented";
    }
}

void
PointToCellMap::requestUpdate()
{
    auto parentPointSet = std::dynamic_pointer_cast<PointSet>(getParentGeometry());
    auto childPointSet  = std::dynamic_pointer_cast<PointSet>(getChildGeometry());

    std::shared_ptr<VecDataArray<double, 3>> parentVerticesPtr = childPointSet->getVertexPositions();
    const VecDataArray<double, 3>&           parentVertices    = *parentVerticesPtr;

    std::shared_ptr<VecDataArray<double, 3>> childVerticesPtr = childPointSet->getVertexPositions();
    VecDataArray<double, 3>&                 childVertices    = *childVerticesPtr;

    ParallelUtils::parallelFor(childVertices.size(),
        [&](const int i)
        {
            const std::vector<VertexWeight>& weights = m_weights[i];
            childVertices[i] = Vec3d::Zero();
            for (int j = 0; j < weights.size(); j++)
            {
                const VertexWeight& weight = weights[j];
                childVertices[i] += parentVertices[weight.first] * weight.second;
            }
            // Account for local displacement from basis
            //childVertices[i] += displacement * rot;
        }, false);

    setOutput(childPointSet);
}

void
PointToCellMap::computeTetWeights(const int posKey, const Vec3d& pos,
                                  std::unordered_map<int, std::vector<VertexWeight>>& vertexWeights)
{
}

void
PointToCellMap::computeTriangleWeights(const int posKey, const Vec3d& pos,
                                       std::unordered_map<int, std::vector<VertexWeight>>& vertexWeights)
{
    auto                                        mesh = std::dynamic_pointer_cast<SurfaceMesh>(getParentGeometry());
    const std::vector<std::unordered_set<int>>& vertexToCells = mesh->getVertexToCells();

    std::shared_ptr<VecDataArray<double, 3>> verticesPtr = mesh->getVertexPositions();
    const VecDataArray<double, 3>&           vertices    = *verticesPtr;
    std::shared_ptr<VecDataArray<int, 3>>    indicesPtr  = mesh->getIndices();
    const VecDataArray<int, 3>&              indices     = *indicesPtr;

    // Find the closest point out of all the cells
    double closestDistanceSqr  = IMSTK_DOUBLE_MAX;
    int    closestCellId       = -1;
    int    closestCellCaseType = -1;
    Vec3d  closestPos = Vec3d::Zero();
    for (int i = 0; i < indices.size(); i++)
    {
        const Vec3d& a = vertices[indices[i][0]];
        const Vec3d& b = vertices[indices[i][1]];
        const Vec3d& c = vertices[indices[i][2]];
        int          caseType = -1;
        const Vec3d  closestPtOnCell =
            CollisionUtils::closestPointOnTriangle(pos, a, b, c, caseType);

        // If closer than the closest
        const double sqrDist = (pos - closestPtOnCell).squaredNorm();
        if (sqrDist < closestDistanceSqr)
        {
            // Record new closest
            closestDistanceSqr  = sqrDist;
            closestCellId       = i;
            closestCellCaseType = caseType;
            closestPos = closestPtOnCell;
        }
    }

    // After determining the closest feature compute weights for the given position

    // Closest to a/0, b/1, c/2
    if (closestCellCaseType == 0 || closestCellCaseType == 1 || closestCellCaseType == 2)
    {
        // If closest point is the vertex itself, project to both neighbor cells and sum their weights
        // normalize later when we have all the weights
        const std::unordered_set<int>& cells = vertexToCells[indices[closestCellId][closestCellCaseType]];
        for (const auto& cellId : cells)
        {
            const int v0 = indices[cellId][0];
            const int v1 = indices[cellId][1];
            const int v2 = indices[cellId][2];

            const Vec3d bcCoord = baryCentric(closestPos, vertices[v0], vertices[v1], vertices[v2]);
            vertexWeights[posKey].push_back({ v0, bcCoord[0] });
            vertexWeights[posKey].push_back({ v1, bcCoord[1] });
            vertexWeights[posKey].push_back({ v2, bcCoord[2] });
        }
    }
    // Closest to edge ab, bc, or ca
    else if (closestCellCaseType == 3 || closestCellCaseType == 4 || closestCellCaseType == 5)
    {
        Vec2i vertexIds = { indices[closestCellId][0], indices[closestCellId][1] }; // ab
        if (closestCellCaseType == 4)
        {
            vertexIds = { indices[closestCellId][1], indices[closestCellId][2] }; // bc
        }
        else if (closestCellCaseType == 5)
        {
            vertexIds = { indices[closestCellId][2], indices[closestCellId][0] }; // ca
        }

        // Find the triangles connected to this edge
        for (int neighborFaceIndex : vertexToCells[indices[closestCellId][0]])
        {
            const Vec3i& cell     = indices[neighborFaceIndex];
            bool         found[2] = { false, false };
            for (int j = 0; j < 3; j++)
            {
                if (cell[j] == vertexIds[0])
                {
                    found[0] = true;
                }
                else if (cell[j] == vertexIds[1])
                {
                    found[1] = true;
                }
            }
            // If it contains both vertices its a face to the edge
            if (found[0] && found[1])
            {
                const int v0 = indices[neighborFaceIndex][0];
                const int v1 = indices[neighborFaceIndex][1];
                const int v2 = indices[neighborFaceIndex][2];

                const Vec3d bcCoord = baryCentric(closestPos, vertices[v0], vertices[v1], vertices[v2]);
                vertexWeights[posKey].push_back({ v0, bcCoord[0] });
                vertexWeights[posKey].push_back({ v1, bcCoord[1] });
                vertexWeights[posKey].push_back({ v2, bcCoord[2] });
            }
        }
    }
    // Closest to face/triangle
    else if (closestCellCaseType == 6)
    {
        const int v0 = indices[closestCellId][0];
        const int v1 = indices[closestCellId][1];
        const int v2 = indices[closestCellId][2];

        const Vec3d bcCoord = baryCentric(closestPos, vertices[v0], vertices[v1], vertices[v2]);
        vertexWeights[posKey].push_back({ v0, bcCoord[0] });
        vertexWeights[posKey].push_back({ v1, bcCoord[1] });
        vertexWeights[posKey].push_back({ v2, bcCoord[2] });
    }
}

void
PointToCellMap::computeSegmentWeights(const int posKey, const Vec3d& pos,
                                      std::unordered_map<int, std::vector<VertexWeight>>& vertexWeights)
{
    auto                                        mesh = std::dynamic_pointer_cast<LineMesh>(getParentGeometry());
    const std::vector<std::unordered_set<int>>& vertexToCells = mesh->getVertexToCells();

    std::shared_ptr<VecDataArray<double, 3>> verticesPtr = mesh->getVertexPositions();
    const VecDataArray<double, 3>&           vertices    = *verticesPtr;
    std::shared_ptr<VecDataArray<int, 2>>    indicesPtr  = mesh->getIndices();
    const VecDataArray<int, 2>&              indices     = *indicesPtr;

    // Find the closest point out of all the cells
    double closestDistanceSqr  = IMSTK_DOUBLE_MAX;
    int    closestCellId       = -1;
    int    closestCellCaseType = -1;
    Vec3d  closestPos = Vec3d::Zero();
    for (int i = 0; i < indices.size(); i++)
    {
        const Vec3d& a = vertices[indices[i][0]];
        const Vec3d& b = vertices[indices[i][1]];
        int          caseType = -1;
        const Vec3d  closestPtOnCell =
            CollisionUtils::closestPointOnSegment(pos, a, b, caseType);

        // If closer than the closest
        const double sqrDist = (pos - closestPtOnCell).squaredNorm();
        if (sqrDist < closestDistanceSqr)
        {
            // Record new closest
            closestDistanceSqr = sqrDist;
            closestCellId      = i;
            closestPos = closestPtOnCell;
            closestCellCaseType = caseType;
        }
    }

    // Closest to edge
    if (closestCellCaseType == 2)
    {
        const int v0 = indices[closestCellId][0];
        const int v1 = indices[closestCellId][1];

        // This is the easiest case
        const Vec2d bcCoord = baryCentric(closestPos, vertices[v0], vertices[v1]);

        vertexWeights[posKey].push_back({ v0, bcCoord[0] });
        vertexWeights[posKey].push_back({ v1, bcCoord[1] });

        // Compute the perp
        Vec3d diff      = vertices[v1] - vertices[v0];
        Mat3d transform = Mat3d::Identity();
        transform.col(0) = diff.cross(Vec3d(diff[1], -diff[0], -diff[2])).normalized();
        transform.col(1) = diff;
        transform.col(2) = transform.col(0).cross(transform.col(1)).normalized();

        const double length0 = transform.col(0).dot(pos - vertices[v0]);
        const double length2 = transform.col(2).dot(pos - vertices[v0]);

        // Apply both of these completely to the sum
        m_orthoWeights[posKey].push_back({ transform.col(0) * length0, 1.0 });
        m_orthoWeights[posKey].push_back({ transform.col(2) * length2, 1.0 });
    }
    // Closest to a/1 or b/0
    else if (closestCellCaseType == 1 || closestCellCaseType == 0)
    {
        // If closest point is the vertex itself, project to both segments and sum their weights
        // normalize later when we have all the weights
        const std::unordered_set<int>& cells = vertexToCells[indices[closestCellId][closestCellCaseType]];
        //const Vec3d closestPt = vertices[indices[closestCellId][closestCellCaseType]];
        for (const auto& cellId : cells)
        {
            const int v0 = indices[cellId][0];
            const int v1 = indices[cellId][1];

            const Vec2d bcCoord = baryCentric(closestPos, vertices[v0], vertices[v1]);
            vertexWeights[posKey].push_back({ v0, bcCoord[0] });
            vertexWeights[posKey].push_back({ v1, bcCoord[1] });

            // Compute the perp
            Vec3d diff      = vertices[v1] - vertices[v0];
            Mat3d transform = Mat3d::Identity();
            transform.col(0) = diff.cross(Vec3d(diff[1], -diff[0], -diff[2])).normalized();
            transform.col(1) = diff.normalized();;
            transform.col(2) = transform.col(0).cross(transform.col(1)).normalized();

            const double length0 = diff.normalized().dot(pos - vertices[v0]);
            const double length2 = transform.col(2).dot(pos - vertices[v0]);

            // Apply both of these completely to the sum
            m_orthoWeights[posKey].push_back({ transform.col(0) * length0, 0.5 });
            m_orthoWeights[posKey].push_back({ transform.col(2) * length2, 0.5 });
        }
    }
}

void
PointToCellMap::normalizeWeights(std::vector<std::pair<int, double>>& vertexWeights)
{
    double sum = 0.0;
    for (int i = 0; i < vertexWeights.size(); i++)
    {
        sum += vertexWeights[i].second;
    }
    sum /= vertexWeights.size();
    for (int i = 0; i < vertexWeights.size(); i++)
    {
        vertexWeights[i].second = vertexWeights[i].second / sum;
    }
}
} // namespace imstk
