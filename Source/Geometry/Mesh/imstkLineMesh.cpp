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

#include "imstkLineMesh.h"
#include "imstkLogger.h"
#include "imstkVecDataArray.h"

namespace imstk
{
void
LineMesh::clear()
{
    CellMesh::clear();

    m_activeCellScalars = "";
}

void
LineMesh::print() const
{
    CellMesh::print();

    LOG(INFO) << "Active Cell Scalars: " << m_activeCellScalars;
}

void
LineMesh::setCellScalars(const std::string& arrayName, std::shared_ptr<AbstractDataArray> scalars)
{
    m_activeCellScalars = arrayName;
    m_cellAttributes[arrayName] = scalars;
}

void
LineMesh::setCellScalars(const std::string& arrayName)
{
    if (hasCellAttribute(arrayName))
    {
        m_activeCellScalars = arrayName;
    }
}

std::shared_ptr<AbstractDataArray>
LineMesh::getCellScalars() const
{
    if (hasCellAttribute(m_activeCellScalars))
    {
        return m_cellAttributes.at(m_activeCellScalars);
    }
    else
    {
        return nullptr;
    }
}
} // namespace imstk