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

#include "imstkPbdConstraint.h"

namespace imstk
{
class PointSet;

///
/// \class PbdState
///
/// \brief State of the body governed by PBD mathematical model
///
/// \todo: Attempt a fixed size state & absbtraction to RigidBodyModel2
///
class PbdState
{
public:
    PbdState() = default;
    virtual ~PbdState() = default;

    ///
    /// \brief Set the state to a given one, copies vector values by value instead of references
    ///
    void setState(std::shared_ptr<PbdState> rhs);

    std::shared_ptr<PbdBody> addBody()
    {
        auto body = std::make_shared<PbdBody>(bodyIter);
        bodyIter++;
        m_bodies.push_back(body);
        m_modified = true;
        return body;
    }

    void removeBody(std::shared_ptr<PbdBody> body)
    {
        auto iter = std::find(m_bodies.begin(), m_bodies.end(), body);
        CHECK(iter != m_bodies.end()) << "removeBody called but could not find PbdyBody in PbdState";
        m_bodies.erase(iter);
    }

    std::shared_ptr<PointSet> getBodyGeometry(const PbdBody& body)
    {
        auto iter = m_bodyGeometries.find(body.bodyHandle);
        if (iter == m_bodyGeometries.end())
        {
            LOG(FATAL) << "Tried to get geometry for body that doesn't exist";
            return nullptr;
        }
        return iter->second;
    }

    void setBodyGeometry(const PbdBody& body, std::shared_ptr<PointSet> geometry)
    {
        m_bodyGeometries[body.bodyHandle] = geometry;
    }

    ///
    /// \brief Initialize the state of the bodies using geometries
    ///
    void initialize();

    ///
    /// \brief Initialize an individual body
    ///
    void initState(PbdBody& body);

public:
    std::vector<std::shared_ptr<PbdBody>> m_bodies;
    int bodyIter = 0; ///< Iterative key for bodies
    ///< The geometries corresponding to each body
    std::unordered_map<int, std::shared_ptr<PointSet>> m_bodyGeometries;
    bool m_modified = true;
};
} // namespace imstk