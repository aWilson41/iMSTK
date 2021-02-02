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

#include "imstkVTKConeRenderDelegate.h"
#include "imstkCone.h"
#include "imstkVisualModel.h"

#include <vtkActor.h>
#include <vtkConeSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkTransform.h>

namespace imstk
{
VTKConeRenderDelegate::VTKConeRenderDelegate(std::shared_ptr<VisualModel> visualModel) : VTKPolyDataRenderDelegate(visualModel)
{
    auto geometry = std::static_pointer_cast<Cone>(visualModel->getGeometry());

    m_coneSource = vtkSmartPointer<vtkConeSource>::New();
    m_coneSource->SetHeight(geometry->getHeight());
    m_coneSource->SetAngle(geometry->getAngle() * 180.0 / PI);
    m_coneSource->SetDirection(0.0, 1.0, 0.0);
    m_coneSource->SetCenter(0.0, geometry->getHeight() * 0.5, 0.0);
    m_coneSource->SetResolution(30);

    // Setup mapper
    {
        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputConnection(m_coneSource->GetOutputPort());
        vtkNew<vtkActor> actor;
        actor->SetMapper(mapper);
        actor->SetUserTransform(m_transform);
        m_mapper = mapper;
        m_actor = actor;
    }

    update();
    updateRenderProperties();
}

void
VTKConeRenderDelegate::processEvents()
{
    VTKRenderDelegate::processEvents();

    // Don't use events for primitives, just always update
    auto geometry = std::static_pointer_cast<Cone>(m_visualModel->getGeometry());

    m_coneSource->SetHeight(geometry->getHeight());
    m_coneSource->SetAngle(geometry->getAngle() * 180.0 / PI);

    AffineTransform3d T = AffineTransform3d::Identity();
    T.translate(geometry->getPosition(Geometry::DataType::PostTransform));
    T.rotate(Quatd::FromTwoVectors(UP_VECTOR, geometry->getOrientationAxis(Geometry::DataType::PostTransform)));
    T.scale(1.0);
    T.matrix().transposeInPlace();

    m_transform->SetMatrix(T.data());
}
} // imstk
