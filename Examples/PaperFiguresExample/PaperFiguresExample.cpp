/*
** This file is part of the Interactive Medical Simulation Toolkit (iMSTK)
** iMSTK is distributed under the Apache License, Version 2.0.
** See accompanying NOTICE for details.
*/

#include "imstkCamera.h"
#include "imstkCapsule.h"
#include "imstkDeviceManager.h"
#include "imstkDeviceManagerFactory.h"
#include "imstkDirectionalLight.h"
#include "imstkGeometryUtilities.h"
#include "imstkIsometricMap.h"
#include "imstkMeshIO.h"
#include "imstkOpenVRDeviceClient.h"
#include "imstkPbdContactConstraint.h"
#include "imstkPbdModel.h"
#include "imstkPbdModelConfig.h"
#include "imstkPbdObject.h"
#include "imstkPbdObjectCollision.h"
#include "imstkPbdObjectController.h"
#include "imstkPbdObjectGrasping.h"
#include "imstkPlane.h"
#include "imstkPortHoleInteraction.h"
#include "imstkRenderMaterial.h"
#include "imstkScene.h"
#include "imstkSceneControlText.h"
#include "imstkSceneManager.h"
#include "imstkSimulationManager.h"
#include "imstkSimulationUtils.h"
#include "imstkSphere.h"
#include "imstkVisualModel.h"
#include "imstkVRCameraControl.h"
#include "imstkVTKOpenVRViewer.h"
#include "imstkVTKViewer.h"
#include "imstkDummyClient.h"
#include "imstkRigidObject2.h"
#include "imstkRigidBodyModel2.h"
#include "imstkRbdConstraint.h"
#include "imstkRigidObjectController.h"

#include "imstkCollisionData.h"
#include "imstkCollisionDetectionAlgorithm.h"
#include "imstkPbdObjectCollision.h"
#include "imstkRigidObjectCollision.h"

using namespace imstk;

const double ks = 10000.0;

std::shared_ptr<PbdObject>
makeSpherePbdObj(const std::string& name)
{
    auto sphereObj = std::make_shared<PbdObject>(name);

    auto sphere = std::make_shared<Sphere>(Vec3d(0.0, 0.0, 0.0), 0.02);

    auto model = std::make_shared<PbdModel>();
    model->getConfig()->m_gravity = Vec3d::Zero();
    model->getConfig()->m_dt = 0.001;
    model->getConfig()->m_doPartitioning = false;
    model->getConfig()->m_linearDampingCoeff = 0.0;
    model->getConfig()->m_angularDampingCoeff = 0.0;
    model->getConfig()->m_iterations = 40;

    sphereObj->setDynamicalModel(model);
    sphereObj->setPhysicsGeometry(sphere);
    sphereObj->setCollidingGeometry(sphere);
    sphereObj->setVisualGeometry(sphere);

    sphereObj->getPbdBody()->setRigid(
        Vec3d(1.0, 0.0, 0.0),
        10.0);

    auto dummyClient = std::make_shared<DummyClient>();
    auto controller = sphereObj->addComponent<PbdObjectController>();
    controller->setDevice(dummyClient);
    controller->setControlledObject(sphereObj);
    controller->setLinearKs(ks);
    controller->setUseForceSmoothening(false);

    return sphereObj;
}

std::shared_ptr<RigidObject2>
makeSphereRbdObj(const std::string& name)
{
    auto sphereObj = std::make_shared<RigidObject2>(name);

    auto sphere = std::make_shared<Sphere>(Vec3d(0.0, 0.0, 0.0), 0.02);

    auto model = std::make_shared<RigidBodyModel2>();
    model->getConfig()->m_gravity = Vec3d::Zero();
    model->getConfig()->m_dt = 0.001;
    model->getConfig()->m_angularVelocityDamping = 1.0;
    model->getConfig()->m_velocityDamping = 1.0;
    model->getConfig()->m_maxNumIterations = 40;

    sphereObj->setDynamicalModel(model);
    sphereObj->setPhysicsGeometry(sphere);
    sphereObj->setCollidingGeometry(sphere);
    sphereObj->setVisualGeometry(sphere);

    sphereObj->getRigidBody()->m_initPos = Vec3d(1.0, 0.0, 0.0);
    sphereObj->getRigidBody()->m_mass = 10.0;

    auto dummyClient = std::make_shared<DummyClient>();
    auto controller = sphereObj->addComponent<RigidObjectController>();
    controller->setControlledObject(sphereObj);
    controller->setDevice(dummyClient);
    controller->setLinearKs(ks);
    controller->setUseForceSmoothening(false);

    return sphereObj;
}

std::shared_ptr<CollidingObject>
makePlaneObj(const std::string& name)
{
    auto plane = std::make_shared<Plane>(Vec3d(0.0, 0.0, 0.0), Vec3d(1.0, 0.0, 0.0));

    auto obj = std::make_shared<CollidingObject>(name);
    obj->setCollidingGeometry(plane);
    obj->setVisualGeometry(plane);

    return obj;
}

int
main()
{
    // Write log to stdout and file
    Logger::startLogger();

    auto scene = std::make_shared<Scene>("PaperFigures");

    std::shared_ptr<PbdObject> pbdObj = makeSpherePbdObj("sphere0");
    scene->addSceneObject(pbdObj);
    std::shared_ptr<RigidObject2> rbdObj = makeSphereRbdObj("sphere1");
    scene->addSceneObject(rbdObj);
    std::shared_ptr<CollidingObject> planeObj = makePlaneObj("plane");
    scene->addSceneObject(planeObj);

    auto pbdCollision = std::make_shared<PbdObjectCollision>(pbdObj, planeObj);
    pbdCollision->setRigidBodyCompliance(0.000001);
    scene->addInteraction(pbdCollision);
    auto rbdCollision = std::make_shared<RigidObjectCollision>(rbdObj, planeObj);
    rbdCollision->setBaumgarteStabilization(1.0);
    scene->addInteraction(rbdCollision);

    // Setup a viewer to render in its own thread
    auto viewer = std::make_shared<VTKViewer>();
    viewer->setActiveScene(scene);
    viewer->setDebugAxesLength(0.1, 0.1, 0.1);

    // Setup a scene manager to advance the scene
    auto sceneManager = std::make_shared<SceneManager>();
    sceneManager->setActiveScene(scene);

    std::ofstream filePbd;
    std::ofstream fileRbd;
    fileRbd.open("C:/Users/Andx_/Desktop/CMBBE_Results/SpringAnalysis/rbdForces.dat");
    filePbd.open("C:/Users/Andx_/Desktop/CMBBE_Results/SpringAnalysis/pbdForces.dat");

    bool firstContactPbd = false;
    bool firstContactRbd = false;
    connect<Event>(sceneManager, &SceneManager::postUpdate,
        [&](Event*)
        {
            const double sceneTime = scene->getSceneTime();
            if (!firstContactPbd && pbdCollision->getCollisionDetection()->getCollisionData()->elementsA.size() > 0)
            {
                printf("Pbd Contact at t: %f\n", sceneTime);
                firstContactPbd = true;
            }
            if (!firstContactRbd && rbdCollision->getCollisionDetection()->getCollisionData()->elementsA.size() > 0)
            {
                printf("Rbd Contact at t: %f\n", sceneTime);
                firstContactRbd = true;
            }
            fileRbd << sceneTime << " " <<
                rbdObj->getComponent<RigidObjectController>()->getDeviceForce().norm() << std::endl;
            filePbd << sceneTime << " " <<
                pbdObj->getComponent<PbdObjectController>()->getDeviceForce().norm() << std::endl;
        });

    auto driver = std::make_shared<SimulationManager>();
    driver->addModule(viewer);
    driver->addModule(sceneManager);
    driver->setDesiredDt(0.001);
    driver->setUseRemainderTimeDivide(false);

    // Add default mouse and keyboard controls to the viewer
    std::shared_ptr<Entity> mouseAndKeyControls =
        SimulationUtils::createDefaultSceneControl(driver);
    scene->addSceneObject(mouseAndKeyControls);

    driver->start();

    filePbd.close();
    fileRbd.close();

    return 0;
}
