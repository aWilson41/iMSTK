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

#include "imstkCollisionPair.h"
#include "imstkCollisionDetection.h"

namespace imstk
{
class PbdObject;

///
/// \class PbdObjectCollisionPair
///
/// \brief This class defines a collision interaction between two PbdObjects.
/// or a PbdObject and CollidingObject.
/// The Pbd Interaction interjects two steps to the model,
/// 1.) The collision detection step, which happens after the "Integrate Position" and
/// "Update Collision Geometry" steps. But before the internal constraints "Solve" step.
/// 2.) The collision constraint solve step, which happens after the internal constraint
/// "Solve" step.
///
class PbdObjectCollisionPair : public CollisionPair
{
public:
    ///
    /// \brief Constructor for two way PbdObject interaction
    ///
    PbdObjectCollisionPair(std::shared_ptr<PbdObject> obj1, std::shared_ptr<PbdObject> obj2,
                           CollisionDetection::Type cdType = CollisionDetection::Type::MeshToMeshBruteForce);

    ///
    /// \brief Constructor for one way PbdObject interaction
    ///
    PbdObjectCollisionPair(std::shared_ptr<PbdObject> obj1, std::shared_ptr<CollidingObject> obj2,
                           CollisionDetection::Type cdType = CollisionDetection::Type::MeshToMeshBruteForce);

    ///
    /// \brief Destructor
    ///
    virtual ~PbdObjectCollisionPair() override = default;

public:
    ///
    /// \brief Applies modifications to TaskGraph
    ///
    void apply() override;

private:
    // Pbd defines two interactions (one at CD and one at solver)
    Inputs  m_solveNodeInputs;
    Outputs m_solveNodeOutputs;
    std::shared_ptr<TaskNode> m_collisionSolveNode = nullptr;
};
}