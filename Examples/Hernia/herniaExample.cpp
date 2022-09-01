/*
** This file is part of the Interactive Medical Simulation Toolkit (iMSTK)
** iMSTK is distributed under the Apache License, Version 2.0.
** See accompanying NOTICE for details.
*/

#include "imstkCamera.h"
#include "imstkCapsule.h"
#include "imstkDirectionalLight.h"
#include "imstkGeometryUtilities.h"
#include "imstkKeyboardDeviceClient.h"
#include "imstkKeyboardSceneControl.h"
#include "imstkLineMesh.h"
#include "imstkLogger.h"
#include "imstkMeshIO.h"
#include "imstkMouseDeviceClient.h"
#include "imstkMouseSceneControl.h"
#include "imstkPbdContactConstraint.h"
#include "imstkPbdHingeJointConstraint.h"
#include "imstkPbdModel.h"
#include "imstkPbdModelConfig.h"
#include "imstkPbdObject.h"
#include "imstkPbdObjectCollision.h"
#include "imstkPbdObjectController.h"
#include "imstkPlane.h"
#include "imstkRenderMaterial.h"
#include "imstkScene.h"
#include "imstkSceneManager.h"
#include "imstkSimulationManager.h"
#include "imstkSphere.h"
#include "imstkVisualModel.h"
#include "imstkVTKViewer.h"
#include "imstkCompositeImplicitGeometry.h"
#include "imstkImplicitGeometryToImageData.h"
#include "imstkSurfaceMeshFlyingEdges.h"
#include "imstkSurfaceMeshSubdivide.h"
#include "imstkOrientedBox.h"
#include "imstkPbdObjectController.h"
#include "imstkDummyClient.h"
#include "imstkDeviceManager.h"
#include "imstkDeviceManagerFactory.h"
#include "imstkPbdObjectGrasping.h"
#include "imstkPointwiseMap.h"

using namespace imstk;

#define HERNIA_RESOURCE_DIR "E:/Hernia/vlahhs/resources/"
///
/// \brief Creates tissue object
/// \param name
/// \param physical dimension of tissue
/// \param dimensions of tetrahedral grid used for tissue
/// \param center of tissue block
///
static std::shared_ptr<PbdObject>
makeTissueObj(const std::string&        name,
              std::shared_ptr<PbdModel> model)
{
    // Setup the Geometry
    std::shared_ptr<TetrahedralMesh> tissueMesh = MeshIO::read<TetrahedralMesh>(HERNIA_RESOURCE_DIR "/FundoTest1/stomach_.vtk");
    std::shared_ptr<SurfaceMesh>     surfMesh   = tissueMesh->extractSurfaceMesh();

    // Setup the Parameters
    model->getConfig()->enableConstraint(PbdModelConfig::ConstraintGenType::Distance, 1000.0);
    model->getConfig()->enableConstraint(PbdModelConfig::ConstraintGenType::Dihedral, 0.1);

    /*model->getConfig()->m_femParams->m_YoungModulus = 50000.0;
    model->getConfig()->m_femParams->m_PoissonRatio = 0.48;
    model->getConfig()->enableFemConstraint(PbdFemConstraint::MaterialType::NeoHookean);*/

    // Setup the material
    auto material = std::make_shared<RenderMaterial>();
    material->setDisplayMode(RenderMaterial::DisplayMode::WireframeSurface);
    material->setBackFaceCulling(false);
    material->setOpacity(0.5);

    // Setup the Object
    auto tissueObj = std::make_shared<PbdObject>(name);
    tissueObj->setVisualGeometry(surfMesh);
    tissueObj->getVisualModel(0)->setRenderMaterial(material);
    tissueObj->setPhysicsGeometry(surfMesh);
    tissueObj->setCollidingGeometry(surfMesh);
    //tissueObj->setPhysicsToCollidingMap(std::make_shared<PointwiseMap>(tissueMesh, surfMesh));
    tissueObj->setDynamicalModel(model);
    tissueObj->getPbdBody()->uniformMassValue = 0.1;
    //model->getConfig()->setBodyDamping(tissueObj->getPbdBody()->bodyHandle, 0.1);

    return tissueObj;
}

///
/// \brief Creates a series of pbd rigid body capsules to stick inside the esophagus
///
void
makeCapsules(std::string                              name,
             std::vector<std::shared_ptr<PbdObject>>& capsuleObjs,
             std::shared_ptr<PbdModel>                model,
             const double                             capsuleRadius)
{
    // Read in a medial line
    auto                                     maLineMesh    = MeshIO::read<LineMesh>(HERNIA_RESOURCE_DIR "/FundoTest1/capsuleRefMesh.obj");
    std::shared_ptr<VecDataArray<double, 3>> maVerticesPtr = maLineMesh->getVertexPositions();
    VecDataArray<double, 3>&                 maVertices    = *maVerticesPtr;
    std::shared_ptr<VecDataArray<int, 2>>    maIndicesPtr  = maLineMesh->getCells();
    VecDataArray<int, 2>&                    maIndices     = *maIndicesPtr;

    // For every line segment in the line mesh create a rigid body capsule
    capsuleObjs.clear();
    int    minZIndex = -1;
    double minZ      = IMSTK_FLOAT_MAX;
    int    maxZIndex = -1;
    double maxZ      = IMSTK_FLOAT_MIN;
    for (int i = 0; i < maIndices.size(); i++)
    {
        const Vec2i& cell   = maIndices[i];
        const Vec3d  p      = maVertices[cell[0]];
        const Vec3d  q      = maVertices[cell[1]];
        const Vec3d  center = (q + p) * 0.5;
        const Vec3d  diff   = q - p;
        const Vec3d  dir    = diff.normalized();

        if (center[2] < minZ)
        {
            minZ      = center[2];
            minZIndex = i;
        }
        if (center[2] > maxZ)
        {
            maxZ      = center[2];
            maxZIndex = i;
        }

        //auto capsule = std::make_shared<Capsule>(Vec3d::Zero(), capsuleRadius, diff.norm());
        auto capsule = std::make_shared<Sphere>(Vec3d::Zero(), capsuleRadius);

        auto capsuleObj = std::make_shared<PbdObject>(name);
        capsuleObj->setDynamicalModel(model);
        capsuleObj->setVisualGeometry(capsule);
        capsuleObj->setPhysicsGeometry(capsule);
        capsuleObj->setCollidingGeometry(capsule);
        capsuleObj->getVisualModel(0)->getRenderMaterial()->setColor(Color::Red);
        capsuleObj->getVisualModel(0)->getRenderMaterial()->setBackFaceCulling(false);
        //capsuleObj->getVisualModel(0)->getRenderMaterial()->setOpacity(0.5);
        //
        // Setup body
        const Quatd orientation = Quatd::FromTwoVectors(Vec3d(0.0, 1.0, 0.0), dir);
        capsuleObj->getPbdBody()->setRigid(center, // Position
            10.0,                                  // Mass
            orientation,                           // Orientation
            Mat3d::Identity() * 0.005);            // Inertia

        capsuleObjs.push_back(capsuleObj);
    }

    // Fix every other capsule via joint point
    //stomachObj->getPbdModel()->getConfig()->addPbdConstraintFunctor(
    //    [=](PbdConstraintContainer& container)
    //    {
    //        const PbdParticleId firstPt =
    //            stomachObj->getPbdModel()->addVirtualParticle(maVertices[maIndices[0][0]], 0.0, Vec3d::Zero(), true);
    //        const PbdParticleId lastPt =
    //            stomachObj->getPbdModel()->addVirtualParticle(maVertices[maIndices[maIndices.size() - 1][1]], 0.0, Vec3d::Zero(), true);

    //        // Fix the top and bottom (min max z axes) capsules
    //        {
    //            const Vec2i& cell1 = maIndices[0];
    //            auto distConstraintStart = std::make_shared<PbdBodyToBodyDistanceConstraint>();
    //            distConstraintStart->initConstraint(
    //                stomachObj->getPbdModel()->getBodies(),
    //                { capsuleObjs[0]->getPbdBody()->bodyHandle, 0 },
    //                maVertices[cell1[0]],
    //                firstPt, 0.0, 0.0001);
    //            container.addConstraint(distConstraintStart);

    //            const Vec2i& cell2 = maIndices[maIndices.size() - 1];
    //            auto distConstraintEnd = std::make_shared<PbdBodyToBodyDistanceConstraint>();
    //            distConstraintEnd->initConstraint(
    //                stomachObj->getPbdModel()->getBodies(),
    //                { capsuleObjs[capsuleObjs.size() - 1]->getPbdBody()->bodyHandle, 0 },
    //                maVertices[cell2[1]],
    //                lastPt, 0.0, 0.0001);
    //            container.addConstraint(distConstraintEnd);
    //        }

    //        for (int i = 0; i < maIndices.size() - 1; i++)
    //        {
    //            const Vec2i& cell1 = maIndices[i];
    //            const Vec3d p1 = maVertices[cell1[0]];
    //            const Vec3d q1 = maVertices[cell1[1]];
    //            const Vec3d center1 = (q1 + p1) * 0.5;
    //            const Vec3d diff1 = q1 - p1;
    //            const Vec3d dir1 = diff1.normalized();

    //            const Vec2i& cell2 = maIndices[i + 1];
    //            const Vec3d p2 = maVertices[cell2[0]];
    //            const Vec3d q2 = maVertices[cell2[1]];
    //            const Vec3d center2 = (q2 + p2) * 0.5;
    //            const Vec3d diff2 = q2 - p2;
    //            const Vec3d dir2 = diff2.normalized();

    //            // Assuming capsule lines are ordered a1 and p1 should be equivalent (in case they aren't average)
    //            const Vec3d jointPt = (q1 + p2) * 0.5;

    //            auto bodyDistToDistConstraint = std::make_shared<PbdBodyToBodyDistanceConstraint>();
    //            bodyDistToDistConstraint->initConstraint(
    //                stomachObj->getPbdModel()->getBodies(),
    //                { capsuleObjs[i]->getPbdBody()->bodyHandle, 0 }, jointPt,
    //                { capsuleObjs[i + 1]->getPbdBody()->bodyHandle, 0 }, jointPt,
    //                0.0, 0.0001);
    //            container.addConstraint(bodyDistToDistConstraint);
    //        }
    //    });
}

int
main()
{
    // Write log to stdout and file
    Logger::startLogger();

    // Setup a scene
    auto scene = std::make_shared<Scene>("Hernia");
    scene->getActiveCamera()->setFocalPoint(0.0158845, -0.00876985, -1.21677);
    scene->getActiveCamera()->setPosition(0.0784281, 0.333918, -1.33171);
    scene->getActiveCamera()->setViewUp(-0.00193505, -0.32117, -0.947019);
    scene->getConfig()->debugCamBoundingBox = false;
    *scene->getCamera("debug") = *scene->getActiveCamera();

    auto pbdModel  = std::make_shared<PbdModel>();
    auto pbdConfig = std::make_shared<PbdModelConfig>();
    pbdConfig->m_gravity = Vec3d(0.0, 0.0, 0.0);
    //pbdConfig->m_gravity = Vec3d(0.0, -9.8, 0.0);
    //pbdConfig->m_gravity = Vec3d(9.8, 0.0, 0.0);
    pbdConfig->m_dt = 0.001;
    pbdConfig->m_iterations = 5;
    pbdConfig->m_linearDampingCoeff  = 0.03;
    pbdConfig->m_angularDampingCoeff = 0.01;
    pbdConfig->m_doPartitioning      = false;
    pbdModel->configure(pbdConfig);

    // 10.0 to 0.1 mass

    auto tissueObj = makeTissueObj("Tissue", pbdModel);
    scene->addSceneObject(tissueObj);

    auto lapTool2 = std::make_shared<PbdObject>("lapTool2");
    {
        const double capsuleLength = 0.3;
        auto         toolGeom      = std::make_shared<Capsule>(Vec3d(0.0, 0.0, 0.0),
            0.002, capsuleLength, Quatd::FromTwoVectors(Vec3d(0.0, 1.0, 0.0), Vec3d(0.0, 0.0, 1.0)));

        lapTool2->setDynamicalModel(pbdModel);
        lapTool2->setPhysicsGeometry(toolGeom);
        lapTool2->setCollidingGeometry(toolGeom);
        lapTool2->setVisualGeometry(toolGeom);

        std::shared_ptr<RenderMaterial> material = lapTool2->getVisualModel(0)->getRenderMaterial();
        material->setIsDynamicMesh(false);
        material->setMetalness(1.0);
        material->setRoughness(0.2);
        material->setShadingModel(RenderMaterial::ShadingModel::PBR);

        lapTool2->getPbdBody()->setRigid(
            Vec3d(0.0, 0.0, capsuleLength * 0.5),
            0.1,
            Quatd::Identity(), Mat3d::Identity() * 10000.0);

        auto controller = lapTool2->addComponent<PbdObjectController>();
        controller->setControlledObject(lapTool2);
        controller->setTranslationOffset((*lapTool2->getPbdBody()->vertices)[0] + Vec3d(0.0, 0.0, -1.3));
        controller->setLinearKs(100.0);
        controller->setAngularKs(100000000.0);
        controller->setForceScaling(0.0);
        controller->setSmoothingKernelSize(15);
        controller->setUseForceSmoothening(true);
    }
    scene->addSceneObject(lapTool2);

    std::vector<std::shared_ptr<PbdObject>> capsuleObjs;
    makeCapsules("CapsuleObjs", capsuleObjs, pbdModel, 0.005);
    for (int i = 0; i < capsuleObjs.size(); i++)
    {
        scene->addSceneObject(capsuleObjs[i]);
        auto capsCollision = std::make_shared<PbdObjectCollision>(
            tissueObj, capsuleObjs[i]);
        /* capsCollision->setFriction(0.0);
         capsCollision->setRestitution(0.0);*/
        capsCollision->setRigidBodyCompliance(0.00001);
        scene->addInteraction(capsCollision);
    }

    // Add picking interaction for both jaws of the tool
    auto grasping = std::make_shared<PbdObjectGrasping>(tissueObj);
    grasping->setStiffness(0.01);
    scene->addInteraction(grasping);

    // Light
    auto light = std::make_shared<DirectionalLight>();
    light->setFocalPoint(Vec3d(5.0, -8.0, -5.0));
    light->setIntensity(1.0);
    scene->addLight("Light", light);

    // Run the simulation
    {
        // Setup a viewer to render
        auto viewer = std::make_shared<VTKViewer>();
        viewer->setVtkLoggerMode(VTKViewer::VTKLoggerMode::MUTE);
        viewer->setActiveScene(scene);
        viewer->setDebugAxesLength(0.01, 0.01, 0.01);

        // Setup a scene manager to advance the scene
        auto sceneManager = std::make_shared<SceneManager>();
        sceneManager->setActiveScene(scene);
        sceneManager->pause(); // Start simulation paused

        auto driver = std::make_shared<SimulationManager>();
        driver->addModule(viewer);
        driver->addModule(sceneManager);
        driver->setDesiredDt(0.003);

        auto hapticManager = DeviceManagerFactory::makeDeviceManager();
        driver->addModule(hapticManager);
        std::shared_ptr<DeviceClient> hapticDevice = hapticManager->makeDeviceClient();
        auto                          controller   = lapTool2->getComponent<PbdObjectController>();
        controller->setDevice(hapticDevice);

        // Add mouse and keyboard controls to the viewer
        {
            auto mouseControl = std::make_shared<MouseSceneControl>();
            mouseControl->setDevice(viewer->getMouseDevice());
            mouseControl->setSceneManager(sceneManager);
            scene->addControl(mouseControl);
            connect<KeyEvent>(viewer->getKeyboardDevice(), &KeyboardDeviceClient::keyPress, [&](KeyEvent* e)
                {
                    if (e->m_key == 'g')
                    {
                        if (pbdModel->getConfig()->m_gravity[0] == 0.0)
                        {
                            pbdModel->getConfig()->m_gravity = Vec3d(0.0, -1.0, 0.0);
                        }
                        else
                        {
                            pbdModel->getConfig()->m_gravity = Vec3d::Zero();
                        }
                    }
                    else if (e->m_key == 'u')
                    {
                        scene->advance(sceneManager->getDt());
                        viewer->update();
                    }
                });
            connect<ButtonEvent>(hapticDevice, &DeviceClient::buttonStateChanged,
                [&](ButtonEvent* e)
                {
                    if (e->m_button == 0 && e->m_buttonState == BUTTON_PRESSED)
                    {
                        LOG(INFO) << "Grasp!";
                        grasping->beginVertexGrasp(std::dynamic_pointer_cast<AnalyticalGeometry>(lapTool2->getCollidingGeometry()));
                        /*grasping->beginCellGrasp(
                            std::dynamic_pointer_cast<AnalyticalGeometry>(lapTool2->getCollidingGeometry()),
                            "SurfaceMeshToCapsuleCD");*/
                    }
                });
            connect<ButtonEvent>(hapticDevice, &DeviceClient::buttonStateChanged,
                [&](ButtonEvent* e)
                {
                    if (e->m_button == 0 && e->m_buttonState == BUTTON_RELEASED)
                    {
                        LOG(INFO) << "Release!";
                        grasping->endGrasp();
                    }
                });

            auto keyControl = std::make_shared<KeyboardSceneControl>();
            keyControl->setDevice(viewer->getKeyboardDevice());
            keyControl->setSceneManager(sceneManager);
            keyControl->setModuleDriver(driver);
            scene->addControl(keyControl);
        }

        driver->start();
    }
}