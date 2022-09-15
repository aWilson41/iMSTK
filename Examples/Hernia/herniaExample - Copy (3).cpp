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
#include "imstkSimulationUtils.h"
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

double tissueCapsuleCompliance = 0.0000001;

double capsuleMass = 10.0;
double capsuleInertiaScale = 0.05;
double capsuleDistConstraintStiffness = 10000.0;
double capsuleHingeConstraintCompliance = 0.01;

double tissueParticleMass = 0.1;
double tissueDistStiffness = 10000.0;
double tissueDihedralStiffness = 1.0;

double globalLinearDamping = 0.05;
double globalAngularDamping = 0.01;

int iterations = 4;
int collisionIterations = 10;
double simTimestep = 0.0003;
double executionTimestep = 0.003;

// We also have tool and grasp stiffness

// Stability:
// XPBD appears to handle multi contact poorly (not as well as pbd). It's hard for me to
// fully understand. In pbd if you have multi contact you can drop stiffness and you will
// get smaller steps towards the solution. It's great not to step immediately to the solution
// particularly when multiple constraints are involved. But xpbd is parameterized by dt such
// that it is consistent in time. Say the solution of x(t), x(1.4)=5 and x(1.5)=10. With dt=1
// Every step we do some 10 iterations. There is one step here at dt=1 so we do 10. Then what
// prevents us from jumping to 10 immediately from 5 on the first iteration. This is bad. I'd
// almost want something that parameterizes over the iteration count sucn that 10 is reached at
// the final iteration. In xpdd though this can be thought of exactly like taking smaller
// timesteps though. This is the motivation of the substeps PBD paper from mueller. Why not
// always just use 1 iteration and lower dt.
// 
// The problem, hard to quantify, is likely the other computation in the scene. PBD does
// reprojection, which is a nifty trick to reach better solutions with less computation.
// Other parts of the scene (namely collision detection and rendering) are not being run
// with extra iterations in the solver. We already do this with rendering and substeps.
// We can take say 1000 physics steps for every 1 render. A similar thing might be required
// with collision detection. Running it at a different rate to allow/improve the usage of single
// iterations in xpbd.

//#define USE_TET

///
/// \brief Compute the a vector of points indices that are coincident with the on the parent
///
static std::vector<int>
computeFixedPtsViaMap(std::shared_ptr<PointSet> parent,
    std::shared_ptr<PointSet> child,
    const double tolerance = 0.00001)
{
    std::vector<int> fixedPts;

    auto map = std::make_shared<PointwiseMap>();
    map->setParentGeometry(parent);
    map->setChildGeometry(child);
    map->setTolerance(tolerance);
    map->compute();
    fixedPts.reserve(child->getNumVertices());
    for (int i = 0; i < child->getNumVertices(); i++)
    {
        fixedPts.push_back(map->getParentVertexId(i));
    }
    return fixedPts;
}

///
/// \brief Creates tissue object
/// \param name
/// \param physical dimension of tissue
/// \param dimensions of tetrahedral grid used for tissue
/// \param center of tissue block
///
static std::shared_ptr<PbdObject>
makeTissueObj(const std::string& name,
    std::shared_ptr<PbdModel> model)
{
    // Setup the Geometry
    std::shared_ptr<TetrahedralMesh> tissueMesh = MeshIO::read<TetrahedralMesh>(HERNIA_RESOURCE_DIR "/FundoTest1/stomach_.vtk");
    std::shared_ptr<SurfaceMesh>     surfMesh = tissueMesh->extractSurfaceMesh();

    // Setup the material
    auto material = std::make_shared<RenderMaterial>();
    material->setDisplayMode(RenderMaterial::DisplayMode::WireframeSurface);
    material->setBackFaceCulling(false);
    material->setOpacity(0.5);

    // Setup the Object
    auto tissueObj = std::make_shared<PbdObject>(name);
    tissueObj->setVisualGeometry(surfMesh);
    tissueObj->setCollidingGeometry(surfMesh);
    tissueObj->getVisualModel(0)->setRenderMaterial(material);
#ifdef USE_TET
    tissueObj->setPhysicsGeometry(tissueMesh);
    tissueObj->setPhysicsToCollidingMap(std::make_shared<PointwiseMap>(tissueMesh, surfMesh));
#else
    tissueObj->setPhysicsGeometry(surfMesh);
#endif

    tissueObj->setDynamicalModel(model);
    tissueObj->getPbdBody()->uniformMassValue = tissueParticleMass;

    auto fixedVerts = MeshIO::read<PointSet>(HERNIA_RESOURCE_DIR"/FundoTest1/fixedVerts.obj");
    PointwiseMap fixedMapper;
    fixedMapper.setParentGeometry(tissueObj->getPhysicsGeometry());
    fixedMapper.setChildGeometry(fixedVerts);
    fixedMapper.compute();
    const std::unordered_map<int, int>& fixedMap = fixedMapper.getMap();
    for (auto i : fixedMap)
    {
        tissueObj->getPbdBody()->fixedNodeIds.push_back(i.second);
    }

#ifdef USE_TET
    model->getConfig()->m_femParams->m_YoungModulus = 10000.0;
    model->getConfig()->m_femParams->m_PoissonRatio = 0.48;
    model->getConfig()->enableFemConstraint(PbdFemConstraint::MaterialType::NeoHookean,
        tissueObj->getPbdBody()->bodyHandle);
#else
    // Setup the Parameters
    model->getConfig()->enableConstraint(PbdModelConfig::ConstraintGenType::Distance, tissueDistStiffness,
        tissueObj->getPbdBody()->bodyHandle);
    model->getConfig()->enableConstraint(PbdModelConfig::ConstraintGenType::Dihedral, tissueDihedralStiffness,
        tissueObj->getPbdBody()->bodyHandle);
#endif

    return tissueObj;
}

///
/// \brief Creates and PbdObject that links/connects to the tissueObj provided
/// via computation of overlapping points
/// 
static std::shared_ptr<PbdObject>
makeLigamentObj(std::string name,
    std::string ligamentGeomFileName,
    std::string fixedVertexFileName,
    const double stiffness,
    std::shared_ptr<PbdObject> tissueObj)
{
    auto ligamentObj = std::make_shared<PbdObject>(name);

    auto surfMesh = MeshIO::read<SurfaceMesh>(ligamentGeomFileName);

    std::shared_ptr<PbdModel> model = tissueObj->getPbdModel();
    ligamentObj->setVisualGeometry(surfMesh);
    ligamentObj->getVisualModel(0)->getRenderMaterial()->setColor(
        Color(183.0 / 255.0, 86.0 / 255.0, 49.0 / 255.0));
    ligamentObj->getVisualModel(0)->getRenderMaterial()->setBackFaceCulling(false);
    ligamentObj->getVisualModel(0)->getRenderMaterial()->setOpacity(0.5);
    ligamentObj->getVisualModel(0)->getRenderMaterial()->setDisplayMode(RenderMaterial::DisplayMode::WireframeSurface);
    ligamentObj->setPhysicsGeometry(surfMesh);
    ligamentObj->setCollidingGeometry(surfMesh);
    ligamentObj->setDynamicalModel(model);
    ligamentObj->getPbdBody()->uniformMassValue = tissueParticleMass;
    auto fixedPtMesh = MeshIO::read<LineMesh>(fixedVertexFileName);
    ligamentObj->getPbdBody()->fixedNodeIds = computeFixedPtsViaMap(surfMesh, fixedPtMesh);

    model->getConfig()->enableConstraint(PbdModelConfig::ConstraintGenType::Distance, 1000.0,
        ligamentObj->getPbdBody()->bodyHandle);
    model->getConfig()->addPbdConstraintFunctor([=](PbdConstraintContainer& container)
        {
            // Setup the map between the tissueObj and current object
            auto tissueToLigamentMap = std::make_shared<PointwiseMap>();
            tissueToLigamentMap->setParentGeometry(tissueObj->getPhysicsGeometry());
            tissueToLigamentMap->setChildGeometry(ligamentObj->getPhysicsGeometry());
            tissueToLigamentMap->setTolerance(0.001);
            std::unordered_map<int, int> tissueMap;
            tissueToLigamentMap->computeMap(tissueMap);

            for (auto i : tissueMap)
            {
                const PbdParticleId id1 = { tissueObj->getPbdBody()->bodyHandle, i.second };
                const PbdParticleId id2 = { ligamentObj->getPbdBody()->bodyHandle, i.first };

                auto constraint = std::make_shared<PbdDistanceConstraint>();
                constraint->initConstraint(
                    tissueObj->getPbdModel()->getBodies().getPosition(id1),
                    tissueObj->getPbdModel()->getBodies().getPosition(id2),
                    id1, id2, 10000.0);
                container.addConstraint(constraint);
            }
        });

    return ligamentObj;
}

///
/// \brief Creates a series of pbd rigid body capsules to stick inside the esophagus
///
void
makeCapsules(const std::string& name,
    std::vector<std::shared_ptr<PbdObject>>& capsuleObjs,
    std::shared_ptr<PbdModel> model,
    const double capsuleRadius)
{
    // Read in a medial line
    auto                                     maLineMesh =
        MeshIO::read<LineMesh>(HERNIA_RESOURCE_DIR "/FundoTest1/capsuleRefMesh.obj");
    std::shared_ptr<VecDataArray<double, 3>> maVerticesPtr = maLineMesh->getVertexPositions();
    VecDataArray<double, 3>& maVertices = *maVerticesPtr;
    std::shared_ptr<VecDataArray<int, 2>>    maIndicesPtr = maLineMesh->getCells();
    VecDataArray<int, 2>& maIndices = *maIndicesPtr;

    // For every line segment in the line mesh create a rigid body capsule
    capsuleObjs.clear();
    int    minZIndex = -1;
    double minZ = IMSTK_FLOAT_MAX;
    int    maxZIndex = -1;
    double maxZ = IMSTK_FLOAT_MIN;
    for (int i = 0; i < maIndices.size(); i++)
    {
        const Vec2i& cell = maIndices[i];
        const Vec3d  p = maVertices[cell[0]];
        const Vec3d  q = maVertices[cell[1]];
        const Vec3d  center = (q + p) * 0.5;
        const Vec3d  diff = q - p;
        const Vec3d  dir = diff.normalized();

        if (center[2] < minZ)
        {
            minZ = center[2];
            minZIndex = i;
        }
        if (center[2] > maxZ)
        {
            maxZ = center[2];
            maxZIndex = i;
        }

        auto capsule = std::make_shared<Capsule>(Vec3d::Zero(), capsuleRadius, diff.norm());
        //auto capsule = std::make_shared<Sphere>(Vec3d::Zero(), capsuleRadius);

        auto capsuleObj = std::make_shared<PbdObject>(name);
        capsuleObj->setDynamicalModel(model);
        capsuleObj->setVisualGeometry(capsule);
        capsuleObj->setPhysicsGeometry(capsule);
        capsuleObj->setCollidingGeometry(capsule);
        capsuleObj->getVisualModel(0)->getRenderMaterial()->setColor(Color::Red);
        capsuleObj->getVisualModel(0)->getRenderMaterial()->setBackFaceCulling(false);
        //capsuleObj->getVisualModel(0)->getRenderMaterial()->setOpacity(0.5);

        // Setup body
        const Quatd orientation = Quatd::FromTwoVectors(Vec3d(0.0, 1.0, 0.0), dir);
        capsuleObj->getPbdBody()->setRigid(center, // Position
            capsuleMass,                                  // Mass
            orientation,                           // Orientation
            Mat3d::Identity() * capsuleInertiaScale);            // Inertia

        capsuleObjs.push_back(capsuleObj);
    }

    // This one just connects every single capsule to be fixed in place in positin
    // and orientation
    model->getConfig()->addPbdConstraintFunctor(
        [=](PbdConstraintContainer& container)
        {
            // For every capsule
            for (int i = 0; i < maIndices.size(); i++)
            {
                auto capsule = std::dynamic_pointer_cast<Capsule>(capsuleObjs[i]->getVisualGeometry());
                const Vec3d center = (*capsuleObjs[i]->getPbdBody()->vertices)[0];
                const Quatd orientation = (*capsuleObjs[i]->getPbdBody()->orientations)[0];
                const Vec3d dir = orientation._transformVector(Vec3d(0.0, 1.0, 0.0));
                const PbdParticleId fixedParticle =
                    model->addVirtualParticle(
                        center, // Vertex
                        0.0, // Mass
                        Vec3d::Zero(),
                        true); // Persist

                auto distConstraint = std::make_shared<PbdDistanceConstraint>();
#ifdef USE_TET
                distConstraint->initConstraint(
                    0.0,
                    { capsuleObjs[i]->getPbdBody()->bodyHandle, 0 },
                    fixedParticle, 500000.0);
#else
                distConstraint->initConstraint(
                    0.0,
                    { capsuleObjs[i]->getPbdBody()->bodyHandle, 0 },
                    fixedParticle, capsuleDistConstraintStiffness);
#endif
                container.addConstraint(distConstraint);

                auto hingeConstraint = std::make_shared<PbdHingeJointConstraint>();
                hingeConstraint->initConstraint({ capsuleObjs[i]->getPbdBody()->bodyHandle, 0 },
                    dir, capsuleHingeConstraintCompliance);
                container.addConstraint(hingeConstraint);


            }
        });
}

std::shared_ptr<PbdObject>
makeLapToolObj(const std::string& name,
    std::shared_ptr<PbdModel> model)
{
    auto lapTool2 = std::make_shared<PbdObject>(name);

    const double capsuleLength = 0.3;
    auto         toolGeom = std::make_shared<Capsule>(Vec3d(0.0, 0.0, 0.0),
        0.002, capsuleLength, Quatd::FromTwoVectors(Vec3d(0.0, 1.0, 0.0), Vec3d(0.0, 0.0, 1.0)));

    lapTool2->setDynamicalModel(model);
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
    controller->setLinearKs(1000.0);
    controller->setAngularKs(300000000.0);
    controller->setForceScaling(0.0);
    controller->setSmoothingKernelSize(15);
    controller->setUseForceSmoothening(true);

    return lapTool2;
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

    auto pbdModel = std::make_shared<PbdModel>();
    std::shared_ptr<PbdModelConfig> pbdConfig = pbdModel->getConfig();
    pbdConfig->m_gravity = Vec3d(0.0, 0.0, 0.0);
    pbdConfig->m_gravity = Vec3d(0.0, -9.8, 0.0);
    //pbdConfig->m_gravity = Vec3d(9.8, 0.0, 0.0);
    pbdConfig->m_dt = simTimestep;
    pbdConfig->m_iterations = iterations;
    pbdConfig->m_collisionIterations = collisionIterations;
    pbdConfig->m_linearDampingCoeff = globalLinearDamping;
    pbdConfig->m_angularDampingCoeff = globalAngularDamping;
    pbdConfig->m_doPartitioning = false;

    std::shared_ptr<PbdObject> tissueObj = makeTissueObj("Tissue", pbdModel);
    scene->addSceneObject(tissueObj);

    std::shared_ptr<PbdObject> lapTool2 = makeLapToolObj("lapTool2", pbdModel);
    scene->addSceneObject(lapTool2);

    std::vector<std::shared_ptr<PbdObject>> capsuleObjs;
    makeCapsules("CapsuleObjs", capsuleObjs, pbdModel, 0.006);
    for (int i = 0; i < capsuleObjs.size(); i++)
    {
        scene->addSceneObject(capsuleObjs[i]);
        auto capsCollision = std::make_shared<PbdObjectCollision>(
            tissueObj, capsuleObjs[i]);
        capsCollision->setRestitution(0.0);
        capsCollision->setRigidBodyCompliance(tissueCapsuleCompliance);
        scene->addInteraction(capsCollision);
    }

    std::shared_ptr<PbdObject> lesserOmentumObj = makeLigamentObj("lesser_omentum",
        HERNIA_RESOURCE_DIR "/FundoTest1/lesser_omentum_sim.obj",
        HERNIA_RESOURCE_DIR "/FundoTest1/lesser_omentum_sim_fixed.obj",
        1000.0,
        tissueObj);
    scene->addSceneObject(lesserOmentumObj);

    std::shared_ptr<PbdObject> greaterOmentumObj =
        makeLigamentObj("greater_omentum",
            HERNIA_RESOURCE_DIR "/FundoTest4/greater_omentum_sim.obj",
            HERNIA_RESOURCE_DIR "/FundoTest4/greater_omentum_sim_fixed.obj",
            1000.0,
            tissueObj);
    scene->addSceneObject(greaterOmentumObj);

    std::shared_ptr<PbdObject> gastrosplenicLigamentObj =
        makeLigamentObj("gastrosplenic_ligament",
            HERNIA_RESOURCE_DIR "/FundoTest1/gastrosplenic_ligament_sim.obj",
            HERNIA_RESOURCE_DIR "/FundoTest1/gastrosplenic_ligament_sim_fixed.obj",
            1000.0,
            tissueObj);
    scene->addSceneObject(gastrosplenicLigamentObj);

    // Add picking interaction for both jaws of the tool
    auto grasping = std::make_shared<PbdObjectGrasping>(tissueObj);
    //grasping->setGeometryToPick(tissueObj->getPhysicsGeometry(), nullptr);
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
        driver->setDesiredDt(executionTimestep);

        auto hapticManager = DeviceManagerFactory::makeDeviceManager();
        driver->addModule(hapticManager);
        std::shared_ptr<DeviceClient> hapticDevice = hapticManager->makeDeviceClient();
        auto                          controller = lapTool2->getComponent<PbdObjectController>();
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
                            std::dynamic_pointer_cast<AnalyticalGeometry>(lapTool2->getCollidingGeometry()));*/
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

            // Add default mouse and keyboard controls to the viewer
            std::shared_ptr<Entity> mouseAndKeyControls =
                SimulationUtils::createDefaultSceneControlEntity(driver);
            scene->addSceneObject(mouseAndKeyControls);
        }

        driver->start();
    }
}