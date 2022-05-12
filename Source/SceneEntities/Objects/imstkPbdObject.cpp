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

#include "imstkPbdObject.h"
#include "imstkLogger.h"
#include "imstkPbdModel.h"
#include "imstkPointSet.h"

namespace imstk
{
std::shared_ptr<PbdModel>
PbdObject::getPbdModel()
{
    m_pbdModel = std::dynamic_pointer_cast<PbdModel>(m_dynamicalModel);
    return m_pbdModel;
}

bool
PbdObject::initialize()
{
    m_pbdModel = std::dynamic_pointer_cast<PbdModel>(m_dynamicalModel);
    if (m_pbdModel == nullptr)
    {
        LOG(FATAL) << "PbdObject " << m_name << " was not given a PbdModel. Please PbdObject::setDynamicalModel";
        return false;
    }

    DynamicObject::initialize();

    return true;
}

void
PbdObject::setDynamicalModel(std::shared_ptr<AbstractDynamicalModel> dynaModel)
{
    // todo: If already has another model, should remove the corresponding body?
    m_pbdModel       = std::dynamic_pointer_cast<PbdModel>(dynaModel);
    m_dynamicalModel = dynaModel;

    // If the model already has a pbd body for this PbdObject remove it from
    // that prior model
    if (m_pbdBody != nullptr)
    {
        CHECK(m_pbdModel != nullptr) <<
            "PbdObject has a PbdBody but cannot find associated PbdModel?";
        m_pbdModel->removePbdBody(m_pbdBody);
    }
    m_pbdBody = m_pbdModel->addPbdBody();
    if (m_physicsGeometry != nullptr)
    {
        auto pointSet = std::dynamic_pointer_cast<PointSet>(m_physicsGeometry);
        CHECK(pointSet != nullptr) << "PbdObject " << m_name << " only supports PointSet geometries";
        m_pbdModel->getCurrentState()->setBodyGeometry(*m_pbdBody, pointSet);
    }
}

void
PbdObject::setPhysicsGeometry(std::shared_ptr<Geometry> geometry)
{
    DynamicObject::setPhysicsGeometry(geometry);
    if (m_pbdModel != nullptr)
    {
        auto pointSet = std::dynamic_pointer_cast<PointSet>(m_physicsGeometry);
        CHECK(pointSet != nullptr) << "PbdObject " << m_name << " only supports PointSet geometries";
        m_pbdModel->getCurrentState()->setBodyGeometry(*m_pbdBody, pointSet);
    }
}
} // namespace imstk