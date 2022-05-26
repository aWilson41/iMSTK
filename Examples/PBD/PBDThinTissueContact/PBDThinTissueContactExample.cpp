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
#include "imstkDirectionalLight.h"
#include "imstkGeometryUtilities.h"
#include "imstkImageData.h"
#include "imstkKeyboardDeviceClient.h"
#include "imstkKeyboardSceneControl.h"
#include "imstkMeshIO.h"
#include "imstkMouseDeviceClient.h"
#include "imstkMouseSceneControl.h"
#include "imstkNew.h"
#include "imstkPbdConstraintFunctor.h"
#include "imstkPbdModel.h"
#include "imstkPbdObject.h"
#include "imstkPbdObjectCollision.h"
#include "imstkRenderMaterial.h"
#include "imstkScene.h"
#include "imstkSceneManager.h"
#include "imstkSimulationManager.h"
#include "imstkVertexLabelVisualModel.h"
#include "imstkVisualModel.h"
#include "imstkVTKViewer.h"

#include "imstkSphere.h"
#include "imstkCapsule.h"

#ifdef iMSTK_USE_OPENHAPTICS
#include "imstkHapticDeviceClient.h"
#include "imstkHapticDeviceManager.h"
#include "imstkPbdObjectController.h"
#else
#include "imstkMouseDeviceClient.h"
#endif

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
    std::shared_ptr<SurfaceMesh> clothMesh =
        GeometryUtils::toTriangleGrid(Vec3d::Zero(),
            Vec2d(width, height), Vec2i(rowCount, colCount));

    // Setup the VisualModel
    imstkNew<RenderMaterial> material;
    material->setBackFaceCulling(false);
    material->setDisplayMode(RenderMaterial::DisplayMode::WireframeSurface);
    material->setShadingModel(RenderMaterial::ShadingModel::PBR);
    auto diffuseTex = MeshIO::read<ImageData>(iMSTK_DATA_ROOT "/textures/fleshDiffuse.jpg");
    material->addTexture(std::make_shared<Texture>(diffuseTex, Texture::Type::Diffuse));
    auto normalTex = MeshIO::read<ImageData>(iMSTK_DATA_ROOT "/textures/fleshNormal.jpg");
    material->addTexture(std::make_shared<Texture>(normalTex, Texture::Type::Normal));
    auto ormTex = MeshIO::read<ImageData>(iMSTK_DATA_ROOT "/textures/fleshORM.jpg");
    material->addTexture(std::make_shared<Texture>(ormTex, Texture::Type::ORM));

    imstkNew<VisualModel> visualModel;
    visualModel->setGeometry(clothMesh);
    visualModel->setRenderMaterial(material);

    imstkNew<VertexLabelVisualModel> labelModel;
    labelModel->setGeometry(clothMesh);
    labelModel->setRenderMaterial(material);

    // Setup the Object
    imstkNew<PbdObject> pbdObject(name);
    pbdObject->addVisualModel(visualModel);
    pbdObject->addVisualModel(labelModel);
    pbdObject->setPhysicsGeometry(clothMesh);
    pbdObject->setCollidingGeometry(clothMesh);
    pbdObject->setDynamicalModel(pbdModel);

    pbdObject->getPbdBody()->uniformMassValue = width * height / (rowCount * colCount);
    for (int x = 0; x < rowCount; x++)
    {
        for (int y = 0; y < colCount; y++)
        {
            if (x == 0 || y == 0 || x == rowCount - 1 || y == colCount - 1)
            {
                pbdObject->getPbdBody()->fixedNodeIds.push_back(x * colCount + y);
            }
        }
    }

    return pbdObject;
}

static std::shared_ptr<PbdObject>
makeToolObj(std::string name)
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
    //auto toolGeom = std::make_shared<Sphere>(Vec3d(0.0, 0.0, 0.0), 1.0);
    auto toolGeom = std::make_shared<Capsule>(Vec3d(0.0, 0.0, 0.0), 1.0, 2.0);

    auto toolObj = std::make_shared<PbdObject>(name);
    toolObj->setDynamicalModel(pbdModel);
    toolObj->setVisualGeometry(toolGeom);
    toolObj->setPhysicsGeometry(toolGeom);
    toolObj->setCollidingGeometry(toolGeom);
    toolObj->getVisualModel(0)->getRenderMaterial()->setDisplayMode(RenderMaterial::DisplayMode::Wireframe);
    toolObj->getVisualModel(0)->getRenderMaterial()->setLineWidth(5.0);
    toolObj->getVisualModel(0)->getRenderMaterial()->setRecomputeVertexNormals(false);
    toolObj->getVisualModel(0)->getRenderMaterial()->setBackFaceCulling(false);

    toolObj->getPbdBody()->uniformMassValue = 1.0;
    toolObj->getPbdBody()->bodyType    = PbdBody::Type::RIGID;
    toolObj->getPbdBody()->initPosTest = Vec3d(1.0, 5.0, 0.0);

    /*auto labelModel = std::make_shared<VertexLabelVisualModel>();
    labelModel->setGeometry(toolGeom);
    labelModel->setRenderMaterial(toolObj->getVisualModel(0)->getRenderMaterial());
    toolObj->addVisualModel(labelModel);*/

    return toolObj;
}

///
/// \brief This example demonstrates collision interaction with a 2d pbd
/// simulated tissue/membrane/cloth
///
int
main()
{
    // Setup logger (write to file and stdout)
    Logger::startLogger();

    // Setup the scene
    imstkNew<Scene> scene("PBDThinTissueContact");
    scene->getActiveCamera()->setPosition(0.12, 4.51, 16.51);
    scene->getActiveCamera()->setFocalPoint(0.0, 0.0, 0.0);
    scene->getActiveCamera()->setViewUp(0.0, 0.96, -0.28);

    // Setup the model used for pbd objects
    // Setup the Parameters
    auto pbdParams = std::make_shared<PbdModelConfig>();
    pbdParams->enableConstraint(PbdModelConfig::ConstraintGenType::Distance, 5000.0);
    pbdParams->enableConstraint(PbdModelConfig::ConstraintGenType::Dihedral, 5000.0);
    pbdParams->m_gravity    = Vec3d(0.0, -70.0, 0.0); // Slightly larger gravity to compensate viscosity
    pbdParams->m_dt         = 0.005;
    pbdParams->m_iterations = 2;
    pbdParams->m_linearDampingCoeff  = 0.05;
    pbdParams->m_angularDampingCoeff = 0.05;

    // Setup the Model
    pbdModel = std::make_shared<PbdModel>();
    pbdModel->configure(pbdParams);

    // Setup a tissue
    std::shared_ptr<PbdObject> tissueObj = makeTissueObj("Tissue", 10.0, 10.0, 5, 5);
    scene->addSceneObject(tissueObj);

    // setup the tool
    std::shared_ptr<PbdObject> toolObj = makeToolObj("Tool");
    scene->addSceneObject(toolObj);

    // Add a collision interaction between the tools
    scene->addInteraction(std::make_shared<PbdObjectCollision>(tissueObj, toolObj, "SurfaceMeshToCapsuleCD"));

    // Light
    imstkNew<DirectionalLight> light;
    light->setFocalPoint(Vec3d(5.0, -8.0, -5.0));
    light->setIntensity(1.0);
    scene->addLight("Light", light);

    // Run the simulation
    {
        // Setup a viewer to render
        imstkNew<VTKViewer> viewer;
        viewer->setVtkLoggerMode(VTKViewer::VTKLoggerMode::MUTE);
        viewer->setActiveScene(scene);

        // Setup a scene manager to advance the scene
        imstkNew<SceneManager> sceneManager;
        sceneManager->setActiveScene(scene);
        sceneManager->pause(); // Start simulation paused

        imstkNew<SimulationManager> driver;
        driver->addModule(viewer);
        driver->addModule(sceneManager);
#ifdef iMSTK_USE_OPENHAPTICS
        //imstkNew<HapticDeviceManager> hapticManager;
        //hapticManager->setSleepDelay(1.0); // Delay for 1ms (haptics thread is limited to max 1000hz)
        //std::shared_ptr<HapticDeviceClient> hapticDeviceClient = hapticManager->makeDeviceClient();
        //driver->addModule(hapticManager);
#endif
        driver->setDesiredDt(0.005);

#ifdef iMSTK_USE_OPENHAPTICS
        //auto controller = std::make_shared<PbdObjectController>(toolObj, hapticDeviceClient);
        //controller->setTranslationScaling(0.1);
        //controller->setForceScaling(0.0);
        //controller->setLinearKs(10.0);
        //controller->setAngularKs(10.0);
        //// Damping doesn't work well here. The force is applied at the start of pbd
        //// Because velocities are ulimately computed after the fact from positions
        //controller->setUseCritDamping(true);
        //scene->addController(controller);
#else
        connect<Event>(sceneManager, &SceneManager::preUpdate, [&](Event*)
            {
                const Vec2d mousePos = viewer->getMouseDevice()->getPos();
                const Vec3d worldPos = Vec3d(mousePos[0] - 0.5, mousePos[1] - 0.5, 0.0) * 10.0;

                toolGeometry->setTranslation(worldPos);
                toolGeometry->postModified();
            });
#endif

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

        driver->start();
    }

    return 0;
}