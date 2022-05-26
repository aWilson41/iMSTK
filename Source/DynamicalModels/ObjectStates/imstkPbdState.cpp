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

#include "imstkPbdState.h"
#include "imstkPointSet.h"

namespace imstk
{
void
PbdStateDummy::setState(std::shared_ptr<PbdStateDummy> src)
{
    if (m_bodies.size() != src->m_bodies.size())
    {
        m_bodies.resize(src->m_bodies.size());
    }
    for (size_t i = 0; i < src->m_bodies.size(); i++)
    {
        // Try not to reallocate the body as it would lose the current buffer
        // pointers
        if (m_bodies[i] == nullptr)
        {
            m_bodies[i] = std::make_shared<PbdBody>();
        }
        m_bodies[i]->deepCopy(*src->m_bodies[i]);
        m_bodies[i]->vertices->postModified();
    }
}

void
PbdStateDummy::initialize()
{
}
} // namespace imstk