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

#include "imstkCone.h"
#include "imstkLogger.h"

namespace imstk
{
void
Cone::print() const
{
    AnalyticalGeometry::print();
    LOG(INFO) << "Height: " << m_height;
    LOG(INFO) << "Angle: " << m_angle;
}

double
Cone::getHeight(DataType type /* = DataType::PostTransform */)
{
    if (type == DataType::PostTransform)
    {
        this->updatePostTransformData();
        return m_heightPostTransform;
    }
    return m_height;
}

void
Cone::applyScaling(const double s)
{
    this->setHeight(m_height * s);
    this->modified();
}

void
Cone::updatePostTransformData() const
{
    if (m_transformApplied)
    {
        return;
    }
    AnalyticalGeometry::updatePostTransformData();
    m_heightPostTransform = m_scaling * m_height;
    m_transformApplied = true;
}
}