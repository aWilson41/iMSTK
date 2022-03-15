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
#include "imstkMath.h"

#include <unordered_map>

namespace imstk
{
class CollidingObject;
class PbdCollisionConstraint;
class PbdObject;

///
/// \class PbdObjectPicking
///
/// \brief This class defines a picking interaction between a PbdObject and
/// a CollidingObject with AnalyticalGeometry. It works with 3 modes.
///
///
class PbdObjectPicking : public SceneObject
{
public:
    enum class Mode
    {
        PickVertex, // Grab a vertex (most performant)
        PickPoint,  // Grab a point on an element
        PickElement // Grab the entire element
        // PickNearest // Grab the point on the nearest element
    };

public:
    PbdObjectPicking(std::shared_ptr<PbdObject> obj1, std::shared_ptr<CollidingObject> obj2, std::string cdType);
    ~PbdObjectPicking() override = default;

    const std::string getTypeName() const override { return "PbdObjectPicking"; }

    ///
    /// \brief Set/Get the picking mode. Indicating how to pick the object
    ///@{
    void setPickingMode(Mode pickingMode) { m_pickingMode = pickingMode; }
    Mode getPickingMode() const { return m_pickingMode; }
    ///@}

    ///
    /// \brief Set/Get the stiffness, 0-1 value that alters the step size in
    /// the solver
    ///@{
    void setStiffness(const double stiffness) { m_stiffness = stiffness; }
    double getStiffness() const { return m_stiffness; }
    ///@}

    ///
    /// \brief End picking (picking will end on next update)
    ///
    void endPick();

    ///
    /// \brief Begin picking (picking will start on next update)
    ///
    void beginPick();

    ///
    /// \brief Compute/generate the constraints for picking
    ///
    void addPickConstraints();

    ///
    /// \brief Remove the constraints for picking
    ///
    void removePickConstraints();

    ///
    /// \brief Add constraint between virtual point
    ///
    virtual void addConstraint(
        Vec3d* fixedPt, Vec3d* fixedPtVel,
        Vec3d* vertex2, double invMass2, Vec3d* velocity2);

    std::shared_ptr<TaskNode> getPickingNode() const { return m_pickingNode; }

    void initGraphEdges(std::shared_ptr<TaskNode> source, std::shared_ptr<TaskNode> sink) override;

protected:
    ///
    /// \brief Update picking state
    ///
    void updatePicking();

    ///
    /// \brief Update the constraints used for picking
    ///
    void updateConstraints();

protected:
    std::shared_ptr<TaskNode> m_pickingNode = nullptr;

    std::shared_ptr<PbdObject>       m_objA = nullptr;
    std::shared_ptr<CollidingObject> m_objB = nullptr;
    std::string m_cdType = "";

    Mode m_pickingMode = Mode::PickVertex;

    bool m_isPicking     = false;
    bool m_isPrevPicking = false;

    double m_stiffness = 0.4;

    std::shared_ptr<CollisionDetectionAlgorithm> m_colDetect = nullptr;

    std::unordered_map<size_t, Vec3d> m_pickedPtIdxOffset; ///> Map for picked nodes.

    std::list<std::tuple<int, Vec3d, Vec3d>> m_constraintPts;
    std::vector<std::shared_ptr<PbdCollisionConstraint>> m_constraints; ///> List of PBD constraints
};
} // namespace imstk