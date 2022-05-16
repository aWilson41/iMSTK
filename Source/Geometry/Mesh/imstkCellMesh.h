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

#include "imstkPointSet.h"

#include "imstkVecDataArray.h"
#include <unordered_set>

#pragma once

namespace imstk
{
///
/// \class CellMesh
///
/// \brief Abstract template base class for all meshes that have homogenous
/// cell types. This class allows templated access to cells.
///
template<int N>
class CellMesh : public PointSet
{
public:
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
        setIndices(indices);
    }

    void clear()
    {
        PointSet::clear();
        if (m_indices != nullptr)
        {
            m_indices->clear();
        }
        m_vertexToCells.clear();
        m_vertexToCells.clear();
        for (auto i : m_cellAttributes)
        {
            i.second->clear();
        }
    }

    ///
    /// \brief Print the surface mesh
    ///
    void print() const override
    {
        PointSet::print();

        LOG(INFO) << "Number of cells: " << getNumCells();
    }

    ///
    /// \brief Computes neighboring triangles for all vertices
    ///
    void computeVertexToCells()
    {
        m_vertexToCells.clear();
        m_vertexToCells.resize(m_vertexPositions->size());

        int cellId = 0;
        for (const auto& cell : *m_indices)
        {
            // \todo: Could variadic unfold
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
    void computeVertexNeighbors()
    {
        m_vertexToNeighborVertex.clear();
        m_vertexToNeighborVertex.resize(m_vertexPositions->size());
        this->computeVertexToCells();

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
    /// \brief Returns true if the geometry is a mesh, else returns false
    ///
    bool isMesh() const override { return true; }

// Accessors
    ///
    /// \brief Get/Set triangle connectivity
    ///@{
    void setIndices(std::shared_ptr<VecDataArray<int, N>> indices) { m_indices = indices; }
    std::shared_ptr<VecDataArray<int, N>> getIndices() const { return m_indices; }
    ///@}

    ///
    /// \brief Returns the number of cells
    ///
    int getNumCells() const { return m_indices->size(); }

    ///
    /// \brief Returns map of vertices to cells that contain the vertex (reverse linkage)
    ///
    const std::vector<std::unordered_set<int>>& getVertexToCells() { return m_vertexToCells; }

    ///
    /// \brief Returns map of vertices to neighboring vertices
    ///
    const std::vector<std::unordered_set<int>>& getVertexNeighbors() { return m_vertexToNeighborVertex; }

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

protected:
    void setCellActiveAttribute(std::string& activeAttributeName, std::string attributeName,
                                const int expectedNumComponents, const ScalarTypeId expectedScalarType)
    {
        std::shared_ptr<AbstractDataArray> attribute = m_cellAttributes[attributeName];
        if (attribute->getNumberOfComponents() != expectedNumComponents)
        {
            LOG(WARNING) << "Failed to set cell attribute on SurfaceMesh " + getName() + " with "
                         << attribute->getNumberOfComponents() << " components. Expected " <<
                expectedNumComponents << " components.";
            return;
        }
        else if (attribute->getScalarType() != expectedScalarType)
        {
            LOG(INFO) << "Tried to set cell attribute on SurfaceMesh " + getName() + " with scalar type "
                      << static_cast<int>(attribute->getScalarType()) << ". Casting to "
                      << static_cast<int>(expectedScalarType) << " scalar type";
            m_cellAttributes[attributeName] = attribute->cast(expectedScalarType);
        }
        activeAttributeName = attributeName;
    }

    std::shared_ptr<VecDataArray<int, N>> m_indices = nullptr;

    std::vector<std::unordered_set<int>> m_vertexToCells;          ///< Map of vertices to neighbor cells
    std::vector<std::unordered_set<int>> m_vertexToNeighborVertex; ///< Map of vertice sto neighbor vertices

    ///< Per cell attributes
    std::unordered_map<std::string, std::shared_ptr<AbstractDataArray>> m_cellAttributes;
};
} // namespace imstk