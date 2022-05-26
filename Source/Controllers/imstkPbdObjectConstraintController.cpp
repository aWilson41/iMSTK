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

#include "imstkPbdObjectConstraintController.h"
#include "imstkDeviceClient.h"
#include "imstkLogger.h"
#include "imstkPbdObject.h"
#include "imstkPbdDistanceConstraint.h"
#include "imstkPbdModel.h"
#include "imstkPbdSolver.h"

namespace imstk
{
PbdObjectConstraintController::PbdObjectConstraintController(const std::string& name) :
    SceneObjectController(name),
    m_distConstraint(std::make_shared<PbdDistanceConstraint>())
{
}

void
PbdObjectConstraintController::setControlledObject(std::shared_ptr<SceneObject> obj)
{
    m_pbdObject = std::dynamic_pointer_cast<PbdObject>(obj);
    CHECK(m_pbdObject != nullptr) << "Controlled object must be a PbdObject";
    CHECK(m_pbdObject->getPbdBody()->bodyType == PbdBody::Type::RIGID)
        << "PbdObjectController can only operate on pbd rigid bodies";
    SceneObjectController::setControlledObject(obj);
}

void
PbdObjectConstraintController::update(const double dt)
{
    if (!updateTrackingData(dt))
    {
        LOG(WARNING) << "warning: could not update tracking info.";
        return;
    }

    if (m_pbdObject == nullptr)
    {
        return;
    }

    // Implementation partially from otaduy lin's paper eq14
    // "A Modular Haptic Rendering Algorithm for Stable and Transparent 6 - DOF Manipulation"
    if (m_deviceClient->getTrackingEnabled() && m_useSpring)
    {
        const PbdParticleId& vid = m_pbdObject->getPbdModel()->addVirtualParticle(getPosition(), 0.0);
        m_distConstraint->initConstraint(0.0, { m_pbdObject->getPbdBody()->bodyHandle, 0 }, vid, m_linearKs.maxCoeff());
        m_constraints.push_back(m_distConstraint.get());
        m_pbdObject->getPbdModel()->getCollisionSolver()->addConstraints(&m_constraints);
    }
    else
    {
        // Zero out external force/torque
        m_pbdObject->getPbdBody()->externalForce  = Vec3d(0.0, 0.0, 0.0);
        m_pbdObject->getPbdBody()->externalTorque = Vec3d(0.0, 0.0, 0.0);
        // Directly set position/rotation
        (*m_pbdObject->getPbdBody()->vertices)[0]     = getPosition();
        (*m_pbdObject->getPbdBody()->orientations)[0] = getOrientation();
    }

    this->postEvent(Event(PbdObjectConstraintController::modified()));
}

void
PbdObjectConstraintController::initGraphEdges(std::shared_ptr<TaskNode> source,
                                              std::shared_ptr<TaskNode> sink)
{
    SceneObjectController::initGraphEdges(source, sink);

    //auto                         pbdObj1 = std::dynamic_pointer_cast<PbdObject>(m_objA);
    //std::shared_ptr<SceneObject> obj2 = m_objB;

    //std::shared_ptr<TaskNode> chNodeAB = m_collisionHandleANode;

    //// Ensure a complete graph
    //m_taskGraph->addEdge(source, pbdObj1->getTaskGraph()->getSource());
    //m_taskGraph->addEdge(pbdObj1->getTaskGraph()->getSink(), sink);
    //m_taskGraph->addEdge(source, obj2->getTaskGraph()->getSource());
    //m_taskGraph->addEdge(obj2->getTaskGraph()->getSink(), sink);

    //// -------------------------------------------------------------------------
    //// Internal Constraint Solve -> Collision Geometry Update -> Collision Detection ->
    //// PbdHandlerAB -> Collision Constraint Solve -> Update Pbd Velocity -> Correct
    //// Velocities for Collision (restitution+friction) -> Pbd Sink
    //// -------------------------------------------------------------------------
    //m_taskGraph->addEdge(pbdObj1->getPbdModel()->getSolveNode(), m_collisionGeometryUpdateNode);
    //m_taskGraph->addEdge(m_collisionGeometryUpdateNode, m_collisionDetectionNode);
    //m_taskGraph->addEdge(m_collisionDetectionNode, chNodeAB); // A=AB=B
    //m_taskGraph->addEdge(chNodeAB, pbdObj1->getPbdModel()->getCollisionSolveNode());

    //if (std::dynamic_pointer_cast<PbdObject>(obj2) == nullptr)
    //{
    //    m_taskGraph->addEdge(obj2->getUpdateGeometryNode(), m_collisionGeometryUpdateNode);
    //    m_taskGraph->addEdge(m_collisionDetectionNode, obj2->getTaskGraph()->getSink());
    //}
}

void
PbdObjectConstraintController::applyForces()
{
    //if (!m_deviceClient->getButton(0))
    //{
    //    // Apply force back to device
    //    if (m_pbdObject != nullptr && m_useSpring)
    //    {
    //        const Vec3d force = -getDeviceForce();
    //        if (m_forceSmoothening)
    //        {
    //            m_forces.push_back(force);
    //            m_forceSum += force;
    //            if (static_cast<int>(m_forces.size()) > m_smoothingKernelSize)
    //            {
    //                m_forceSum -= m_forces.front();
    //                m_forces.pop_front();
    //            }
    //            const Vec3d avgForce = m_forceSum / m_forces.size();

    //            // Render only the spring force (not the other forces the body has)
    //            m_deviceClient->setForce(avgForce);
    //        }
    //        else
    //        {
    //            // Render only the spring force (not the other forces the body has)
    //            m_deviceClient->setForce(force);
    //        }
    //    }
    //}
    //else
    //{
    //    m_deviceClient->setForce(Vec3d(0.0, 0.0, 0.0));
    //}
}
} // namespace imstk