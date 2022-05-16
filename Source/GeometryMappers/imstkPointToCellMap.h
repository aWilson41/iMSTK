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

#include "imstkGeometryMap.h"
#include "imstkMath.h"
#include "imstkTypes.h"

namespace imstk
{
template<typename T, int N> class VecDataArray;

///
/// \class PointToCellMap
///
/// \brief PointToCellMap can map a PointSet to a LineMesh, SurfaceMesh,
/// or TetrahedralMesh.
///
/// This effectively does automatic weight computation for linear blend
/// skinning but without any rotations. It can then apply the parent
/// deformation to the child.
///
/// Unlike linear blend skinning we do not have bones with basis'. There
/// are no rotations. So if you map a SurfaceMesh to a LineMesh and rotate
/// the LineMesh, the SurfaceMesh will not rotate.
///
/// Usage is to have a detailed surface mesh deform according to a simulated
/// mesh.
///  - Can be used to deform a detailed thread SurfaceMesh with a simulated LineMesh
///  - Can be used to deform a detailed thick SurfaceMesh with a thin simulated SurfaceMesh
///  - Can be used to deform a detailed SurfaceMesh with a coarse TetrahedralMesh
///  - Can be used to deform a LineMesh (such as vessels) according to a TetrahedralMesh tissue.
///  - WIP: can be used to deform a TetrahedralMesh according to a LineMesh.
///
/// For points coincident (in or on) to a cell (line, triangle, or tet)
/// a barycentric interpolation is used for weights (bc weights).
///
/// For points not coincident (in or on) to a cell. We find the nearest
/// point of all cells and compute bc weights.
///
class PointToCellMap : public GeometryMap
{
public:
    PointToCellMap();
    PointToCellMap(
        std::shared_ptr<Geometry> parent,
        std::shared_ptr<Geometry> child);
    ~PointToCellMap() override = default;

    IMSTK_TYPE_NAME(PointToCellMap)

    ///
    /// \brief Compute the tetra-triangle mesh map
    ///
    void compute() override;

    void computeTetWeights();
    void computeTriangleWeights();
    void computeSegmentWeights();

protected:
    ///
    /// \brief Apply (if active) the tetra-triangle mesh map
    ///
    void requestUpdate() override;

    ///
    /// \brief Find the closest tetrahedron based on the distance to their
    /// centroids for a given point in 3D space
    ///
    void computeTetWeights(const int posKey, const Vec3d& pos,
                           std::unordered_map<int, std::vector<std::pair<int, double>>>& vertexWeights);

    ///
    /// \brief Find the closest triangle based on the distance
    /// \return triangle id and case type (face-0,edge-1,vertex-2)
    ///
    void computeTriangleWeights(const int posKey, const Vec3d& pos,
                                std::unordered_map<int, std::vector<std::pair<int, double>>>& vertexWeights);

    ///
    /// \brief Find the closest  line based on the distance
    /// \return line id and case type (edge-1,vertex-2)
    ///
    void computeSegmentWeights(const int posKey, const Vec3d& pos,
                               std::unordered_map<int, std::vector<std::pair<int, double>>>& vertexWeights);

    ///
    /// \brief Sum then divide by 1
    ///
    void normalizeWeights(std::vector<std::pair<int, double>>& vertexWeights);

protected:
    using VertexWeight      = std::pair<int, double>;   ///< Vertex id that contributes, weight it contributes
    using VertexOrthoWeight = std::pair<Vec3d, double>; ///< Local difference transformed along the ortho, weight it contributes

    // Child vertex to weights map (weight being parent vertex id, bc weight)
    // child vertex = vertex_i * weight_i
    std::unordered_map<int, std::vector<VertexWeight>> m_weights;

    // Orthogonal weights should be summed with normal weights
    // These are locally defined and transformed according to the cell
    std::unordered_map<int, std::vector<VertexOrthoWeight>> m_orthoWeights;

    // Parent cell/bone id -> origin and transform
    std::vector<Vec3d> m_origin;    ///< Bone origin
    std::vector<Mat3d> m_transform; ///< Bone Transform
};
} // namespace imstk
