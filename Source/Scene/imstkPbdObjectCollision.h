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

#include "imstkCollisionInteraction.h"
#include "imstkMacros.h"

namespace imstk
{
class PbdObject;

///
/// \class PbdObjectCollision
///
/// \brief This class defines a collision interaction between two PbdObjects
/// or PbdObject & CollidingObject
///
class PbdObjectCollision : public CollisionInteraction
{
public:
    ///
    /// \brief Constructor for PbdObject-PbdObject or PbdObject-CollidingObject collisions
    ///
    PbdObjectCollision(std::shared_ptr<PbdObject> obj1, std::shared_ptr<CollidingObject> obj2,
                       std::string cdType = "ClosedSurfaceMeshToMeshCD");
    ~PbdObjectCollision() override = default;

    IMSTK_TYPE_NAME(PbdObjectCollision)

    void setRestitution(const double restitution);
    double getRestitution() const;

    void setFriction(const double friction);
    double getFriction() const;

    void setCollisionIterations(const int iterations);
    int getCollisionIterations() const;

    ///
    /// \brief Setup connectivity of task graph
    ///
    void initGraphEdges(std::shared_ptr<TaskNode> source, std::shared_ptr<TaskNode> sink) override;
};
} // namespace imstk