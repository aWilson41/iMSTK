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

#include "imstkVirtualCouplingObject.h"

#include <memory>

#include "imstkGeometry.h"
#include "imstkGeometryMap.h"

#include <g3log/g3log.hpp>

namespace imstk {

void
VirtualCouplingObject::updateFromDevice()
{
    Vec3d p;
    Quatd r;

    if (!this->computeTrackingData(p, r))
    {
        LOG(WARNING) << "VirtualCouplingObject::updateFromDevice warning: could not update tracking info.";
        return;
    }

    // Update colliding geometry
    m_collidingGeometry->setPosition(p);
    m_collidingGeometry->setOrientation(r);

    // Update visual geometry
    if(m_collidingToVisualMap)
    {
        m_collidingToVisualMap->apply();
    }
}
}
