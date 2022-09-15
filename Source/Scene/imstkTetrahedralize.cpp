/*
** This file is part of the Interactive Medical Simulation Toolkit (iMSTK)
** iMSTK. It is not distributed under the same license.
*/

#include "imstkTetrahedralize.h"
#include "imstkLogger.h"
#include "imstkMeshIO.h"
#include "imstkSurfaceMesh.h"
#include "imstkTetrahedralMesh.h"
#include <vtksys/SystemTools.hxx>

namespace imstk
{
Tetrahedralize::Tetrahedralize()
{
    setNumInputPorts(1);
    setRequiredInputType<SurfaceMesh>(0);

    setNumOutputPorts(1);
    setOutput(std::make_shared<TetrahedralMesh>(), 0);
}

void
Tetrahedralize::setInputMesh(std::shared_ptr<SurfaceMesh> mesh)
{
    setInput(mesh, 0);
}

std::shared_ptr<TetrahedralMesh>
Tetrahedralize::getOutputTetMesh()
{
    return std::dynamic_pointer_cast<TetrahedralMesh>(getOutput());
}

static unsigned long
getHashedSizeOf(std::shared_ptr<SurfaceMesh> surfMesh)
{
    unsigned long inputSize = 0;
    inputSize += sizeof(surfMesh);
    std::shared_ptr<VecDataArray<double, 3>> verticesPtr = surfMesh->getVertexPositions();
    std::shared_ptr<VecDataArray<int, 3>>    cellsPtr    = surfMesh->getCells();
    inputSize += sizeof(double) * verticesPtr->size() * 3;
    inputSize += sizeof(int) * cellsPtr->size() * 3;
    return inputSize;
}

void
Tetrahedralize::requestUpdate()
{
    std::shared_ptr<SurfaceMesh>     inputSurfaceMesh = std::dynamic_pointer_cast<SurfaceMesh>(getInput(0));
    std::shared_ptr<TetrahedralMesh> outputTetMesh    = std::dynamic_pointer_cast<TetrahedralMesh>(getOutput(0));

    // Write the input we are about to give it and compare if the inputs differ from last input
    bool readFromExisting = false;
    vtksys::SystemTools::MakeDirectory("tetWild");
    MeshIO::write(inputSurfaceMesh, "tetWild/tetWildStagedInput.stl");
    if (vtksys::SystemTools::FileExists("tetWild/tetWildInput.stl"))
    {
        if (!vtksys::SystemTools::FilesDiffer("tetWild/tetWildInput.stl", "tetWild/tetWildStagedInput.stl"))
        {
            readFromExisting = true;
        }
    }

    // If input is the same and parameters are the same, the output will be the same, there
    // is no need to compute it again. Just read in the existing one.
    if (readFromExisting)
    {
        auto resultsMesh = MeshIO::read<TetrahedralMesh>("tetWild/tetWildInput_.msh");
        *outputTetMesh = *resultsMesh;
        return;
    }
    else
    {
        // Write the input
        MeshIO::write(inputSurfaceMesh, "tetWild/tetWildInput.stl");
        //system("TetWild --help");

        std::stringstream commandStream;
        commandStream << "TetWild --input tetWild/tetWildInput.stl -q"; // q for quiet
        if (m_idealEdgeLength != -1.0)
        {
            commandStream << " --ideal-edge-length " << m_idealEdgeLength;
        }
        if (m_targetNumVertices != -1)
        {
            commandStream << " --targeted-num-v " << m_targetNumVertices;
        }
        system(commandStream.str().c_str());

        // Read the results
        auto resultsMesh = MeshIO::read<TetrahedralMesh>("tetWild/tetWildInput_.msh");

        *outputTetMesh = *resultsMesh;
    }
}
} // namespace imstk