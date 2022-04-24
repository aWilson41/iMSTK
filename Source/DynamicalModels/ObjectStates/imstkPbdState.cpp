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
PbdState::setState(std::shared_ptr<PbdState> src)
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

    bodyIter = src->bodyIter;
    m_bodyGeometries = src->m_bodyGeometries;
    m_modified       = src->m_modified;
}

void
PbdState::initialize()
{
    // Only run initialize if a body has been added/removed
    if (!m_modified)
    {
        return;
    }

    for (auto& body : m_bodies)
    {
        initState(*body);
    }

    m_modified = false;
}

void
PbdState::initState(PbdBody& body)
{
    auto iter = m_bodyGeometries.find(body.bodyHandle);
    CHECK(iter != m_bodyGeometries.end()) << "No geometry assigned to pbd body";
    std::shared_ptr<PointSet> mesh = iter->second;
    body.vertices     = mesh->getVertexPositions();
    body.prevVertices = std::make_shared<VecDataArray<double, 3>>(*body.vertices);

    const int numParticles = body.vertices->size();

    // Initialize Mass+InvMass
    {
        // If the input mesh has per vertex masses, use those
        std::shared_ptr<AbstractDataArray> masses = mesh->getVertexAttribute("Mass");
        if (masses != nullptr && masses->getNumberOfComponents() == 1
            && masses->getScalarType() == IMSTK_DOUBLE && masses->size() == numParticles)
        {
            body.masses = std::dynamic_pointer_cast<DataArray<double>>(masses);
            body.invMasses->resize(body.masses->size());
            for (int i = 0; i < body.masses->size(); i++)
            {
                (*body.invMasses)[i] = ((*body.masses)[i] == 0.0) ? 0.0 : 1.0 / (*body.masses)[i];
            }
        }
        // If not, initialize as uniform and put on mesh
        else
        {
            // Initialize as uniform
            body.masses    = std::make_shared<DataArray<double>>(numParticles);
            body.invMasses = std::make_shared<DataArray<double>>(numParticles);

            body.masses->fill(body.uniformMassValue);
            body.invMasses->fill((body.uniformMassValue != 0.0) ?
                1.0 / body.uniformMassValue : 0.0);

            mesh->setVertexAttribute("Mass", body.masses);
        }
        mesh->setVertexAttribute("InvMass", body.invMasses);
    }

    // Initialize Velocities
    {
        // If the input mesh has per vertex velocities, use those
        std::shared_ptr<AbstractDataArray> velocities = mesh->getVertexAttribute("Velocities");
        if (velocities != nullptr && velocities->getNumberOfComponents() == 3 && velocities->getScalarType() == IMSTK_DOUBLE
            && std::dynamic_pointer_cast<VecDataArray<double, 3>>(velocities)->size() == numParticles)
        {
            body.velocities = std::dynamic_pointer_cast<VecDataArray<double, 3>>(velocities);
        }
        // If not, put existing (0 initialized velocities) on mesh
        else
        {
            body.velocities = std::make_shared<VecDataArray<double, 3>>(numParticles);
            body.velocities->fill(Vec3d::Zero());
            mesh->setVertexAttribute("Velocities", body.velocities);
        }
    }

    // Initialize Accelerations
    {
        // If the input mesh has per vertex accelerations, use those
        std::shared_ptr<AbstractDataArray> accelerations = mesh->getVertexAttribute("Accelerations");
        if (accelerations != nullptr && accelerations->getNumberOfComponents() == 3 && accelerations->getScalarType() == IMSTK_DOUBLE
            && std::dynamic_pointer_cast<VecDataArray<double, 3>>(accelerations)->size() == numParticles)
        {
            body.accelerations = std::dynamic_pointer_cast<VecDataArray<double, 3>>(accelerations);
        }
        // If not, put existing (0 initialized velocities) on mesh
        else
        {
            body.accelerations = std::make_shared<VecDataArray<double, 3>>(numParticles);
            body.accelerations->fill(Vec3d::Zero());
            mesh->setVertexAttribute("Accelerations", body.accelerations);
        }
    }

    // Overwrite some masses for specified fixed points
    body.fixedNodeInvMass = std::unordered_map<int, double>();
    for (auto i : body.fixedNodeIds)
    {
        DataArray<double>& invMasses = *body.invMasses;
        CHECK(i < numParticles && i >= 0) << "Tried to fix particle " << i
                                          << " but there only exist " << numParticles << " particles";
        body.fixedNodeInvMass[i] = invMasses[i];
        invMasses[i] = 0.0;
    }
}
} // namespace imstk