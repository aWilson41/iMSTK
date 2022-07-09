/*=========================================================================
   Library: iMSTK
   Copyright (c) Kitware, Inc.

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

#include "imstkPointSet.h"
#include "imstkVecDataArray.h"

#include <unordered_set>

#pragma once

namespace imstk
{
///
/// \class AbstractCellMesh
///
/// \brief Provides non templated base for cell based meshes
///
class AbstractCellMesh : public PointSet
{
public:
    AbstractCellMesh() = default;
    ~AbstractCellMesh() override = default;

    ///
    /// \brief Returns true if the geometry is a mesh, else returns false
    ///
    bool isMesh() const override { return true; }

    void clear() override;

    ///
    /// \brief Print the surface mesh
    ///
    void print() const override;

    virtual int getNumCells() const = 0;

    ///
    /// \brief Computes neighboring cells for all vertices
    ///
    virtual void computeVertexToCellMap() { }

    ///
    /// \brief Computes neighboring vertices for all vertices
    ///
    virtual void computeVertexNeighbors() { }

    ///
    /// \brief Get cells as abstract array. Overridden by derived classes to return
    ///     cells as point indices.
    ///
    virtual std::shared_ptr<AbstractDataArray> getIndices() const = 0;

    ///
    /// \brief Returns map of vertices to cells that contain the vertex (reverse linkage)
    ///
    const std::vector<std::unordered_set<int>>& getVertexToCellMap() const { return m_vertexToCells; }

    ///
    /// \brief Returns map of vertices to neighboring vertices
    ///
    const std::vector<std::unordered_set<int>>& getVertexNeighbors() const { return m_vertexToNeighborVertex; }

// Attributes
    ///
    /// \brief Get the cell attributes map
    ///
    const std::unordered_map<std::string, std::shared_ptr<AbstractDataArray>>& getCellAttributes() const { return m_cellAttributes; }

    void setCellAttribute(const std::string& arrayName, std::shared_ptr<AbstractDataArray> arr)
    {
        m_cellAttributes[arrayName] = arr;
    }

    std::shared_ptr<AbstractDataArray> getCellAttribute(const std::string& arrayName) const
    {
        auto it = m_cellAttributes.find(arrayName);
        if (it == m_cellAttributes.end())
        {
            LOG(WARNING) << "No array with such name holds any cell data.";
            return nullptr;
        }
        return it->second;
    }

    ///
    /// \brief Check if a specific data array exists.
    ///
    bool hasCellAttribute(const std::string& arrayName) const
    {
        return (m_cellAttributes.find(arrayName) != m_cellAttributes.end());
    }

    ///
    /// \brief Set the cell attributes map
    ///
    void setCellAttributes(std::unordered_map<std::string, std::shared_ptr<AbstractDataArray>> attributes) { m_cellAttributes = attributes; }

    ///
    /// \brief Get/Set the active scalars
    ///@{
    void setCellScalars(const std::string& arrayName, std::shared_ptr<AbstractDataArray> scalars);
    void setCellScalars(const std::string& arrayName);
    std::string getActiveCellScalars() const { return m_activeCellScalars; }
    std::shared_ptr<AbstractDataArray> getCellScalars() const;
    ///@}

    ///
    /// \brief Get/Set the active normals
    ///@{
    void setCellNormals(const std::string& arrayName, std::shared_ptr<VecDataArray<double, 3>> normals);
    void setCellNormals(const std::string& arrayName);
    std::string getActiveCellNormals() const { return m_activeCellNormals; }
    std::shared_ptr<VecDataArray<double, 3>> getCellNormals() const;
    ///@}

    ///
    /// \brief Get/Set the active tangents
    ///@{
    void setCellTangents(const std::string& arrayName, std::shared_ptr<VecDataArray<double, 3>> tangents);
    void setCellTangents(const std::string& arrayName);
    std::string getActiveCellTangents() const { return m_activeCellTangents; }
    std::shared_ptr<VecDataArray<double, 3>> getCellTangents() const;
///@}

protected:
    void setCellActiveAttribute(std::string& activeAttributeName, std::string attributeName,
                                const int expectedNumComponents, const ScalarTypeId expectedScalarType);

    std::vector<std::unordered_set<int>> m_vertexToCells;          ///< Map of vertices to neighbor cells
    std::vector<std::unordered_set<int>> m_vertexToNeighborVertex; ///< Map of vertice sto neighbor vertices

    ///< Per cell attributes
    std::unordered_map<std::string, std::shared_ptr<AbstractDataArray>> m_cellAttributes;

    std::string m_activeCellNormals  = "";
    std::string m_activeCellTangents = "";
    std::string m_activeCellScalars  = "";
};

///
/// \class CellMesh
///
/// \brief Abstract template base class for all meshes that have homogenous
/// cell types. This class allows templated access to cells.
///
template<int N>
class CellMesh : public AbstractCellMesh
{
public:
    static constexpr int CellVertexCount = N;

    CellMesh() : m_indices(std::make_shared<VecDataArray<int, N>>()) { }
    ~CellMesh() override = default;

    ///
    /// \brief Initializes the rest of the data structures given vertex
    /// positions and connectivity
    ///
    void initialize(std::shared_ptr<VecDataArray<double, 3>> vertices,
                    std::shared_ptr<VecDataArray<int, N>> indices)
    {
        clear();

        PointSet::initialize(vertices);
        setCells(indices);
    }

    void clear() override
    {
        AbstractCellMesh::clear();
        if (m_indices != nullptr)
        {
            m_indices->clear();
        }
    }

    ///
    /// \brief Computes neighboring triangles for all vertices
    ///
    void computeVertexToCellMap() override
    {
        m_vertexToCells.clear();
        m_vertexToCells.resize(m_vertexPositions->size());

        int cellId = 0;
        for (const auto& cell : *m_indices)
        {
            // \todo: Could variadic unfold
            // Subclasses can implement a more efficient version if needbe
            for (int i = 0; i < N; i++)
            {
                m_vertexToCells.at(cell[i]).insert(cellId);
            }
            cellId++;
        }
    }

    ///
    /// \brief Computes neighboring vertices for all vertices
    ///
    void computeVertexNeighbors() override
    {
        m_vertexToNeighborVertex.clear();
        m_vertexToNeighborVertex.resize(m_vertexPositions->size());
        this->computeVertexToCellMap();

        // For every vertex
        const VecDataArray<int, N>& indices = *m_indices;
        for (int vertexId = 0; vertexId < m_vertexToNeighborVertex.size(); vertexId++)
        {
            // For every cell it is connected too
            for (const int cellId : m_vertexToCells.at(vertexId))
            {
                // For every vertex of that cell
                for (int i = 0; i < N; i++)
                {
                    // So long as its not the source vertex (not a neighbor of itself)
                    const int vertexId2 = indices[cellId][i];
                    if (vertexId2 != vertexId)
                    {
                        m_vertexToNeighborVertex.at(vertexId).insert(vertexId2);
                    }
                }
            }
        }
    }

    ///
    /// \brief compute the barycentric weights of a given point in 3D space for a given the cell
    ///
    virtual Eigen::Vector<double, N> computeBarycentricWeights(const int    imstkNotUsed(cellId),
                                                               const Vec3d& imstkNotUsed(pos)) const
    {
        return Eigen::Vector<double, N>::Zero();
    }

    std::shared_ptr<AbstractDataArray> getIndices() const override { return m_indices; }

    ///
    /// \brief Get/Set cell connectivity
    ///@{
    void setCells(std::shared_ptr<VecDataArray<int, N>> indices) { m_indices = indices; }
    std::shared_ptr<VecDataArray<int, N>> getCells() const { return m_indices; }
    ///@}

    ///
    /// \brief Returns the number of cells
    ///
    int getNumCells() const override { return m_indices->size(); }

protected:
    std::shared_ptr<VecDataArray<int, N>> m_indices = nullptr;
};
} // namespace imstk