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

#include "imstkCamera.h"
#include "imstkCapsule.h"
#include "imstkDirectionalLight.h"
#include "imstkGeometryUtilities.h"
#include "imstkHapticDeviceClient.h"
#include "imstkHapticDeviceManager.h"
#include "imstkImageData.h"
#include "imstkKeyboardDeviceClient.h"
#include "imstkKeyboardSceneControl.h"
#include "imstkLaparoscopicToolController.h"
#include "imstkMeshIO.h"
#include "imstkMouseDeviceClient.h"
#include "imstkMouseSceneControl.h"
#include "imstkNew.h"
#include "imstkPbdModel.h"
#include "imstkPbdObject.h"
#include "imstkPbdObjectCollision.h"
#include "imstkPbdObjectGrasping.h"
#include "imstkRenderMaterial.h"
#include "imstkScene.h"
#include "imstkSceneManager.h"
#include "imstkSimulationManager.h"
#include "imstkSurfaceMesh.h"
#include "imstkVisualModel.h"
#include "imstkVTKViewer.h"
#include "imstkSphere.h"
#include "imstkVertexLabelVisualModel.h"

#include "imstkPbdObjectConstraintController.h"

using namespace imstk;

std::shared_ptr<PbdModel> pbdModel = nullptr;

///
/// \brief Creates tissue object
///
static std::shared_ptr<PbdObject>
makeTissueObj(const std::string& name,
              const double       width,
              const double       height,
              const int          rowCount,
              const int          colCount)
{
    // Setup the Geometry
    std::shared_ptr<SurfaceMesh> mesh =
        GeometryUtils::toTriangleGrid(Vec3d::Zero(),
            Vec2d(width, height), Vec2i(rowCount, colCount),
            Quatd::Identity(), 2.0);

    // Setup the VisualModel
    imstkNew<RenderMaterial> material;
    material->setBackFaceCulling(false);
    material->setDisplayMode(RenderMaterial::DisplayMode::Surface);
    material->setShadingModel(RenderMaterial::ShadingModel::PBR);
    auto diffuseTex = MeshIO::read<ImageData>(iMSTK_DATA_ROOT "/textures/fleshDiffuse.jpg");
    material->addTexture(std::make_shared<Texture>(diffuseTex, Texture::Type::Diffuse));
    auto normalTex = MeshIO::read<ImageData>(iMSTK_DATA_ROOT "/textures/fleshNormal.jpg");
    material->addTexture(std::make_shared<Texture>(normalTex, Texture::Type::Normal));
    auto ormTex = MeshIO::read<ImageData>(iMSTK_DATA_ROOT "/textures/fleshORM.jpg");
    material->addTexture(std::make_shared<Texture>(ormTex, Texture::Type::ORM));

    // Setup the Object
    imstkNew<PbdObject> tissueObj(name);
    tissueObj->setVisualGeometry(mesh);
    tissueObj->getVisualModel(0)->setRenderMaterial(material);
    tissueObj->setPhysicsGeometry(mesh);
    tissueObj->setCollidingGeometry(mesh);
    tissueObj->setDynamicalModel(pbdModel);
    for (int x = 0; x < rowCount; x++)
    {
        for (int y = 0; y < colCount; y++)
        {
            if (x == 0 || y == 0 || x == rowCount - 1 || y == colCount - 1)
            {
                tissueObj->getPbdBody()->fixedNodeIds.push_back(x * colCount + y);
            }
        }
    }
    tissueObj->getPbdBody()->uniformMassValue = 0.1;

    pbdModel->getConfig()->enableConstraint(PbdModelConfig::ConstraintGenType::Distance, 1000.0,
        tissueObj->getPbdBody()->bodyHandle);
    pbdModel->getConfig()->enableConstraint(PbdModelConfig::ConstraintGenType::Dihedral, 1.0,
        tissueObj->getPbdBody()->bodyHandle);

    /*imstkNew<VertexLabelVisualModel> model;
    model->setGeometry(mesh);
    tissueObj->addVisualModel(model);*/

    return tissueObj;
}

static std::shared_ptr<PbdObject>
makeRbdObj(std::string name, Vec3d pos)
{
    // Setup the tool to press the tissue
    /* auto toolGeom = std::make_shared<LineMesh>();
     VecDataArray<double, 3> vertices(2);
     vertices[0] = Vec3d(0.0, 0.0, 0.0);
     vertices[1] = Vec3d(0.0, 2.0, 0.0);
     VecDataArray<int, 2> indices(1);
     indices[0] = Vec2i(0, 1);
     toolGeom->initialize(
         std::make_shared<VecDataArray<double, 3>>(vertices),
         std::make_shared<VecDataArray<int, 2>>(indices));
 #ifndef iMSTK_USE_OPENHAPTICS
     toolGeometry->translate(Vec3d(0.5, 2.0, 0.5));
 #endif*/
    auto toolGeom = std::make_shared<Sphere>(Vec3d(0.0, 0.0, 0.0), 0.02);
    //auto toolGeom = std::make_shared<Capsule>(Vec3d(0.0, 0.0, 0.0), 1.0, 2.0);

    auto toolObj = std::make_shared<PbdObject>(name);
    toolObj->setDynamicalModel(pbdModel);
    toolObj->setVisualGeometry(toolGeom);
    toolObj->setPhysicsGeometry(toolGeom);
    toolObj->setCollidingGeometry(toolGeom);
    toolObj->getVisualModel(0)->getRenderMaterial()->setIsDynamicMesh(false);
    toolObj->getVisualModel(0)->getRenderMaterial()->setRecomputeVertexNormals(false);
    toolObj->getVisualModel(0)->getRenderMaterial()->setBackFaceCulling(false);

    toolObj->getPbdBody()->uniformMassValue = 0.05;
    toolObj->getPbdBody()->bodyType    = PbdBody::Type::RIGID;
    toolObj->getPbdBody()->initPosTest = pos;

    /*auto labelModel = std::make_shared<VertexLabelVisualModel>();
    labelModel->setGeometry(toolGeom);
    labelModel->setRenderMaterial(toolObj->getVisualModel(0)->getRenderMaterial());
    toolObj->addVisualModel(labelModel);*/

    return toolObj;
}

///
/// \brief This example demonstrates Pbd grasping. PbdObjectGrasping allows
/// us to hold onto parts of a tissue or other pbd deformable with a tool
///
int
main()
{
    // Setup logger (write to file and stdout)
    Logger::startLogger();

    // Scene
    imstkNew<Scene> scene("PbdThinTissueGraspingExample");
    scene->getActiveCamera()->setPosition(0.001, 0.07, 0.25);
    scene->getActiveCamera()->setFocalPoint(0.0, 0.0, 0.0);
    scene->getActiveCamera()->setViewUp(0.0, 0.96, -0.28);

    // Setup the Parameters
    imstkNew<PbdModelConfig> pbdParams;
    //pbdParams->m_gravity = Vec3d(0.0, -8.0, 0.0);
    pbdParams->m_gravity    = Vec3d(0.0, 0.0, 0.0);
    pbdParams->m_dt         = 0.005;
    pbdParams->m_iterations = 8;
    pbdParams->m_collisionIterations = 10;
    pbdParams->m_angularDampingCoeff = 0.1;
    pbdParams->m_linearDampingCoeff  = 0.1;
    pbdParams->m_contactStiffness    = 0.1;

    // Setup the Model
    pbdModel = std::make_shared<PbdModel>();
    //pbdModel->setModelGeometry(mesh);
    pbdModel->configure(pbdParams);

    std::shared_ptr<PbdObject> rbdObj = makeRbdObj("whateve", Vec3d(0.0, 0.1, 0.0));
    scene->addSceneObject(rbdObj);

    /*imstkNew<Capsule> geomShaft;
    geomShaft->setLength(1.0);
    geomShaft->setRadius(0.005);
    geomShaft->setOrientation(Quatd(Rotd(PI_2, Vec3d(1.0, 0.0, 0.0))));
    geomShaft->setTranslation(Vec3d(0.0, 0.0, 0.5));
    imstkNew<CollidingObject> objShaft("ShaftObject");
    objShaft->setVisualGeometry(MeshIO::read<SurfaceMesh>(iMSTK_DATA_ROOT "/laptool/pivot.obj"));
    objShaft->setCollidingGeometry(geomShaft);
    scene->addSceneObject(objShaft);

    imstkNew<Capsule> geomUpperJaw;
    geomUpperJaw->setLength(0.05);
    geomUpperJaw->setTranslation(Vec3d(0.0, 0.0013, -0.016));
    geomUpperJaw->setRadius(0.004);
    geomUpperJaw->setOrientation(Quatd(Rotd(PI_2, Vec3d(1.0, 0.0, 0.0))));
    imstkNew<CollidingObject> objUpperJaw("UpperJawObject");
    objUpperJaw->setVisualGeometry(MeshIO::read<SurfaceMesh>(iMSTK_DATA_ROOT "/laptool/upper.obj"));
    objUpperJaw->setCollidingGeometry(geomUpperJaw);
    scene->addSceneObject(objUpperJaw);

    imstkNew<Capsule> geomLowerJaw;
    geomLowerJaw->setLength(0.05);
    geomLowerJaw->setTranslation(Vec3d(0.0, -0.0013, -0.016));
    geomLowerJaw->setRadius(0.004);
    geomLowerJaw->setOrientation(Quatd(Rotd(PI_2, Vec3d(1.0, 0.0, 0.0))));
    imstkNew<CollidingObject> objLowerJaw("LowerJawObject");
    objLowerJaw->setVisualGeometry(MeshIO::read<SurfaceMesh>(iMSTK_DATA_ROOT "/laptool/lower.obj"));
    objLowerJaw->setCollidingGeometry(geomLowerJaw);
    scene->addSceneObject(objLowerJaw);

    imstkNew<Capsule> pickGeom;
    pickGeom->setLength(0.05);
    pickGeom->setTranslation(Vec3d(0.0, 0.0, -0.016));
    pickGeom->setRadius(0.006);
    pickGeom->setOrientation(Quatd(Rotd(PI_2, Vec3d(1.0, 0.0, 0.0))));*/

    // 300mm x 300mm patch of tissue
    std::shared_ptr<PbdObject> tissueObj = makeTissueObj("Tissue", 0.1, 0.1, 16, 16);
    scene->addSceneObject(tissueObj);

    imstkNew<HapticDeviceManager> deviceManager;
    deviceManager->setSleepDelay(1.0);
    std::shared_ptr<HapticDeviceClient> client = deviceManager->makeDeviceClient();

    // Create and add virtual coupling object controller in the scene
    /*imstkNew<LaparoscopicToolController> controller;
    controller->setParts(objShaft, objUpperJaw, objLowerJaw, pickGeom);
    controller->setDevice(client);
    controller->setJawAngleChange(1.0);
    controller->setTranslationScaling(0.001);
    scene->addController(controller);*/

    auto rbdGhost    = std::make_shared<SceneObject>("ghost");
    auto ghostSphere = std::make_shared<Sphere>(Vec3d(0.0, 0.0, 0.0), 0.02);
    rbdGhost->setVisualGeometry(ghostSphere);
    rbdGhost->getVisualModel(0)->getRenderMaterial()->setColor(Color::Red);
    rbdGhost->getVisualModel(0)->getRenderMaterial()->setOpacity(0.5);
    scene->addSceneObject(rbdGhost);

    // Create and add virtual coupling object controller in the scene
    /*imstkNew<LaparoscopicToolController> controller(objShaft, objUpperJaw, objLowerJaw, pickGeom, client);
    controller->setJawAngleChange(1.0);
    controller->setTranslationScaling(0.001);
    scene->addController(controller);*/

    // Add collision for both jaws of the tool
    /*auto upperJawCollision = std::make_shared<PbdObjectCollision>(tissueObj, objUpperJaw, "SurfaceMeshToCapsuleCD");
    auto lowerJawCollision = std::make_shared<PbdObjectCollision>(tissueObj, objLowerJaw, "SurfaceMeshToCapsuleCD");
    scene->addInteraction(upperJawCollision);
    scene->addInteraction(lowerJawCollision);*/

    // Add picking interaction for both jaws of the tool
    auto jawPicking = std::make_shared<PbdObjectGrasping>(tissueObj);
    scene->addInteraction(jawPicking);

    // Add a collision interaction between the tools
    //scene->addInteraction(std::make_shared<PbdObjectCollision>(rbdObj, tissueObj, "SurfaceMeshToSphereCD"));

    // Light
    imstkNew<DirectionalLight> light;
    light->setFocalPoint(Vec3d(.0, -1.0, -1.0));
    light->setIntensity(1.0);
    scene->addLight("light", light);

    // Run the simulation
    {
        // Setup a viewer to render
        imstkNew<VTKViewer> viewer;
        viewer->setActiveScene(scene);
        viewer->setDebugAxesLength(0.01, 0.01, 0.01);

        // Setup a scene manager to advance the scene
        imstkNew<SceneManager> sceneManager;
        sceneManager->setActiveScene(scene);
        sceneManager->pause(); // Start simulation paused

        imstkNew<SimulationManager> driver;
        driver->addModule(deviceManager);
        driver->addModule(viewer);
        driver->addModule(sceneManager);
        driver->setDesiredDt(0.005);

        // Add mouse and keyboard controls to the viewer
        {
            auto mouseControl = std::make_shared<MouseSceneControl>();
            mouseControl->setDevice(viewer->getMouseDevice());
            mouseControl->setSceneManager(sceneManager);
            viewer->addControl(mouseControl);

            auto keyControl = std::make_shared<KeyboardSceneControl>();
            keyControl->setDevice(viewer->getKeyboardDevice());
            keyControl->setSceneManager(sceneManager);
            keyControl->setModuleDriver(driver);
            viewer->addControl(keyControl);
        }

        connect<Event>(sceneManager, &SceneManager::postUpdate,
            [&](Event*)
            {
                // Simulate the cloth in real time
                tissueObj->getPbdModel()->getConfig()->m_dt = sceneManager->getDt();
            });
        //connect<Event>(sceneManager, &SceneManager::postUpdate,
        //    [&](Event*)
        //    {
        //        //ghostMat->setOpacity(std::min(1.0, controller->getDeviceForce().norm() / 15.0));

        //        // Also apply controller transform to ghost geometry
        //        ghostSphere->setTranslation(controller->getPosition());
        //        ghostSphere->setRotation(controller->getOrientation());
        //        ghostSphere->updatePostTransformData();
        //        ghostSphere->postModified();
        //    });

        //connect<Event>(controller, &LaparoscopicToolController::JawClosed,
        //    [&](Event*)
        //    {
        //        LOG(INFO) << "Jaw Closed!";

        //       /* upperJawCollision->setEnabled(false);
        //        lowerJawCollision->setEnabled(false);*/
        //        jawPicking->beginCellGrasp(pickGeom, "SurfaceMeshToCapsuleCD");
        //    });
        //connect<Event>(controller, &LaparoscopicToolController::JawOpened,
        //    [&](Event*)
        //    {
        //        LOG(INFO) << "Jaw Opened!";

        //        /*upperJawCollision->setEnabled(true);
        //        lowerJawCollision->setEnabled(true);*/
        //        jawPicking->endGrasp();
        //    });

        driver->start();
    }

    return 0;
}
