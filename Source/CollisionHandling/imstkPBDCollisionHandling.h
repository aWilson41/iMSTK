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

#include "imstkCollisionHandling.h"
#include "imstkMath.h"

#include <vector>

namespace imstk
{
class CollidingObject;
class PbdObject;
class PbdCollisionConstraint;
class PbdEdgeEdgeConstraint;
class PbdPointTriangleConstraint;
class PbdPointPointConstraint;
class PbdCollisionSolver;
struct CollisionData;

///
/// \class PBDCollisionHandling
///
/// \brief Implements PBD based collision handling
///
class PBDCollisionHandling : public CollisionHandling
{
public:
    ///
    /// \brief Constructor
    ///
    PBDCollisionHandling(const std::shared_ptr<CollisionData> colData,
                         std::shared_ptr<PbdObject>           pbdObject1,
                         std::shared_ptr<PbdObject>           pbdObject2);

    ///
    /// \brief Constructor
    ///
    PBDCollisionHandling(const std::shared_ptr<CollisionData> colData,
                         std::shared_ptr<PbdObject>           pbdObject1,
                         std::shared_ptr<CollidingObject>     collidingObject2);

    ///
    /// \brief Destructor, clear memory pool
    ///
    virtual ~PBDCollisionHandling() override;

public:
    ///
    /// \brief Compute forces based on collision data
    ///
    void processCollisionData() override;

    std::shared_ptr<PbdCollisionSolver> getCollisionSolver() const { return m_pbdCollisionSolver; }

protected:
    void generateConstraints();

private:
    std::shared_ptr<PbdObject>       m_obj1 = nullptr; ///> PBD object
    std::shared_ptr<CollidingObject> m_obj2 = nullptr;

    std::shared_ptr<PbdCollisionSolver> m_pbdCollisionSolver = nullptr;

    std::vector<PbdCollisionConstraint*> m_PBDConstraints; ///> List of PBD constraints

    // Types locked in because we won't reallocate pools
    std::vector<PbdEdgeEdgeConstraint*> m_EEConstraintPool;
    size_t m_EEConstraintPoolSize = 0;
    std::vector<PbdPointTriangleConstraint*> m_VTConstraintPool;
    size_t m_VTConstraintPoolSize = 0;
    std::vector<PbdPointPointConstraint*> m_PPConstraintPool;
    size_t m_PPConstraintPoolSize = 0;

    std::vector<std::pair<Vec3d, double>> m_fixedVertexPairs;
};
}
