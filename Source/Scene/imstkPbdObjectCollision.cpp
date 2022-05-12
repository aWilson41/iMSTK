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

#include "imstkPbdObjectCollision.h"
#include "imstkCDObjectFactory.h"
#include "imstkCollisionData.h"
#include "imstkCCDAlgorithm.h"
#include "imstkCollisionDetectionAlgorithm.h"
#include "imstkPbdCollisionHandling.h"
#include "imstkPbdModel.h"
#include "imstkPbdObject.h"
#include "imstkPbdSolver.h"
#include "imstkTaskGraph.h"

namespace imstk
{
PbdObjectCollision::PbdObjectCollision(std::shared_ptr<PbdObject> obj1, std::shared_ptr<CollidingObject> obj2,
                                       std::string cdType) :
    CollisionInteraction("PbdObjectCollision_" + obj1->getName() + "_vs_" + obj2->getName(), obj1, obj2)
{
    std::shared_ptr<PbdModel> pbdModel1 = obj1->getPbdModel();

    // Setup the CD
    std::shared_ptr<CollisionDetectionAlgorithm> cd = CDObjectFactory::makeCollisionDetection(cdType);
    cd->setInput(obj1->getCollidingGeometry(), 0);
    cd->setInput(obj2->getCollidingGeometry(), 1);
    setCollisionDetection(cd);

    // Setup the handler
    std::shared_ptr<PbdCollisionHandling> ch = std::make_shared<PbdCollisionHandling>();
    ch->setInputObjectA(obj1);
    ch->setInputObjectB(obj2);
    ch->setInputCollisionData(cd->getCollisionData());

    setCollisionHandlingAB(ch);

    m_taskGraph->addNode(obj1->getTaskGraph()->getSource());
    m_taskGraph->addNode(obj1->getTaskGraph()->getSink());
    m_taskGraph->addNode(obj2->getTaskGraph()->getSource());
    m_taskGraph->addNode(obj2->getTaskGraph()->getSink());

    m_taskGraph->addNode(pbdModel1->getSolveNode());
    m_taskGraph->addNode(pbdModel1->getCollisionSolveNode());

    if (auto pbdObj2 = std::dynamic_pointer_cast<PbdObject>(obj2))
    {
        CHECK(pbdModel1 == pbdObj2->getPbdModel()) << "PbdObjectCollision may only be used with PbdObjects that share the same PbdModel";
    }
    else
    {
        m_taskGraph->addNode(obj2->getUpdateGeometryNode());
    }
}

void
PbdObjectCollision::setRestitution(const double restitution)
{
    std::dynamic_pointer_cast<PbdCollisionHandling>(getCollisionHandlingA())->setRestitution(restitution);
}

double
PbdObjectCollision::getRestitution() const
{
    return std::dynamic_pointer_cast<PbdCollisionHandling>(getCollisionHandlingA())->getRestitution();
}

void
PbdObjectCollision::setFriction(const double friction)
{
    std::dynamic_pointer_cast<PbdCollisionHandling>(getCollisionHandlingA())->setFriction(friction);
}

double
PbdObjectCollision::getFriction() const
{
    return std::dynamic_pointer_cast<PbdCollisionHandling>(getCollisionHandlingA())->getFriction();
}

void
PbdObjectCollision::initGraphEdges(std::shared_ptr<TaskNode> source, std::shared_ptr<TaskNode> sink)
{
    CollisionInteraction::initGraphEdges(source, sink);

    auto                         pbdObj1 = std::dynamic_pointer_cast<PbdObject>(m_objA);
    std::shared_ptr<SceneObject> obj2    = m_objB;

    std::shared_ptr<TaskNode> chNodeAB = m_collisionHandleANode;

    // Ensure a complete graph
    m_taskGraph->addEdge(source, pbdObj1->getTaskGraph()->getSource());
    m_taskGraph->addEdge(pbdObj1->getTaskGraph()->getSink(), sink);
    m_taskGraph->addEdge(source, obj2->getTaskGraph()->getSource());
    m_taskGraph->addEdge(obj2->getTaskGraph()->getSink(), sink);

    // -------------------------------------------------------------------------
    // Internal Constraint Solve -> Collision Geometry Update -> Collision Detection ->
    // PbdHandlerAB -> Collision Constraint Solve -> Update Pbd Velocity -> Correct
    // Velocities for Collision (restitution+friction) -> Pbd Sink
    // -------------------------------------------------------------------------
    m_taskGraph->addEdge(pbdObj1->getPbdModel()->getSolveNode(), m_collisionGeometryUpdateNode);
    m_taskGraph->addEdge(m_collisionGeometryUpdateNode, m_collisionDetectionNode);
    m_taskGraph->addEdge(m_collisionDetectionNode, chNodeAB); // A=AB=B
    m_taskGraph->addEdge(chNodeAB, pbdObj1->getPbdModel()->getCollisionSolveNode());

    if (std::dynamic_pointer_cast<PbdObject>(obj2) == nullptr)
    {
        m_taskGraph->addEdge(obj2->getUpdateGeometryNode(), m_collisionGeometryUpdateNode);
        m_taskGraph->addEdge(m_collisionDetectionNode, obj2->getTaskGraph()->getSink());
    }
}
} // namespace imstk