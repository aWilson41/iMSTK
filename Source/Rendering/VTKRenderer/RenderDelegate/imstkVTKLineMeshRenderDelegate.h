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

#include "imstkVTKRenderDelegate.h"

class vtkDoubleArray;

namespace imstk
{
class LineMesh;

///
/// \class VTKLineMeshRenderDelegate
///
/// \brief
///
class VTKLineMeshRenderDelegate : public VTKRenderDelegate
{
public:
    ///
    /// \brief Constructor
    ///
    explicit VTKLineMeshRenderDelegate(std::shared_ptr<VisualModel> visualModel);

    virtual ~VTKLineMeshRenderDelegate() override = default;

    ///
    /// \brief Update line source based on the line mesh
    ///
    void updateDataSource() override;

protected:

    vtkSmartPointer<vtkDoubleArray> m_mappedVertexArray;    ///> Mapped array of vertices
};
}
