/*
** This file is part of the Interactive Medical Simulation Toolkit (iMSTK)
** iMSTK. It is not distributed under the same license.
*/

#pragma once

#include "imstkGeometryAlgorithm.h"
#include "imstkMath.h"

namespace imstk
{
class SurfaceMesh;
class TetrahedralMesh;

///
/// \class Tetrahedralize
///
/// \brief This filter uses TetWild to tetrahedralize meshes
///
class Tetrahedralize : public GeometryAlgorithm
{
public:
    Tetrahedralize();
    ~Tetrahedralize() override = default;

public:
    ///
    /// \brief Required input, port 0
    ///
    void setInputMesh(std::shared_ptr<SurfaceMesh> surfMesh);

    std::shared_ptr<TetrahedralMesh> getOutputTetMesh();

    double getIdealEdgeLength() const { return m_idealEdgeLength; }
    void setIdealEdgeLength(const double idealEdgeLength) { m_idealEdgeLength = idealEdgeLength; }

    int getTargetNumVertices() const { return m_targetNumVertices; }
    void setTargetNumVertices(const int targetNumVertices) { m_targetNumVertices = targetNumVertices; }

    ///
    /// \brief Get/Set whether to use caching. If on, the input surface mesh
    /// file will be searched for in a cache location and if it is the same
    /// as the one we are generating for. It will read it instead of computing
    /// it again.
    ///@{
    const bool getUseCaching() const { return m_useCaching; }
    void setUseCaching(const bool useCaching) { m_useCaching = useCaching; }
///@}

protected:
    void requestUpdate() override;

private:
    double m_idealEdgeLength   = -1.0;
    int    m_targetNumVertices = -1;
    bool   m_useCaching = false;
};
} // namespace imstk