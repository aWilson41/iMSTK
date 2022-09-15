/*
** This file is part of the Interactive Medical Simulation Toolkit (iMSTK)
** iMSTK is distributed under the Apache License, Version 2.0.
** See accompanying NOTICE for details.
*/

#include "CutHelp.h"
#include "imstkCamera.h"
#include "imstkDeviceManager.h"
#include "imstkDeviceManagerFactory.h"
#include "imstkDirectionalLight.h"
#include "imstkGeometryUtilities.h"
#include "imstkIsometricMap.h"
#include "imstkKeyboardDeviceClient.h"
#include "imstkKeyboardSceneControl.h"
#include "imstkMeshIO.h"
#include "imstkMouseDeviceClient.h"
#include "imstkMouseSceneControl.h"
#include "imstkPbdModel.h"
#include "imstkPbdModelConfig.h"
#include "imstkPbdObject.h"
#include "imstkPbdObjectCollision.h"
#include "imstkPbdObjectController.h"
#include "imstkPlane.h"
#include "imstkPointwiseMap.h"
#include "imstkRenderMaterial.h"
#include "imstkScene.h"
#include "imstkSceneManager.h"
#include "imstkSimulationManager.h"
#include "imstkSimulationUtils.h"
#include "imstkVisualModel.h"
#include "imstkVTKViewer.h"

using namespace imstk;

static void
addDummyVertexPointSet(std::shared_ptr<PointSet> pointSet)
{
    // Add a dummy vertex to the vertices
    std::shared_ptr<VecDataArray<double, 3>> verticesPtr = pointSet->getVertexPositions();
    VecDataArray<double, 3>&                 vertices    = *verticesPtr;
    vertices.resize(vertices.size() + 1);
    for (int i = vertices.size() - 1; i >= 1; i--)
    {
        vertices[i] = vertices[i - 1];
    }
    vertices[0] = Vec3d(0.0, 0.0, 0.0);
    pointSet->setInitialVertexPositions(std::make_shared<VecDataArray<double, 3>>(*verticesPtr));
}

static void
addDummyVertex(std::shared_ptr<SurfaceMesh> surfMesh)
{
    addDummyVertexPointSet(surfMesh);

    // Then shift all indices by 1
    std::shared_ptr<VecDataArray<int, 3>> indicesPtr = surfMesh->getCells();
    VecDataArray<int, 3>&                 indices    = *indicesPtr;
    for (int i = 0; i < indices.size(); i++)
    {
        indices[i][0]++;
        indices[i][1]++;
        indices[i][2]++;
    }
}

static void
addDummyVertex(std::shared_ptr<TetrahedralMesh> tetMesh)
{
    addDummyVertexPointSet(tetMesh);

    // Then shift all indices by 1
    std::shared_ptr<VecDataArray<int, 4>> tissueIndicesPtr = tetMesh->getCells();
    VecDataArray<int, 4>&                 tissueIndices    = *tissueIndicesPtr;
    for (int i = 0; i < tissueIndices.size(); i++)
    {
        tissueIndices[i][0]++;
        tissueIndices[i][1]++;
        tissueIndices[i][2]++;
        tissueIndices[i][3]++;
    }
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
              const Vec3d& size, const Vec3i& dim, const Vec3d& center,
              std::shared_ptr<PbdModel> model)
{
    // Setup the Geometry
    std::shared_ptr<TetrahedralMesh> tissueMesh = GeometryUtils::toTetGrid(center, size, dim);
    std::shared_ptr<SurfaceMesh>     surfMesh   = tissueMesh->extractSurfaceMesh();

    addDummyVertex(tissueMesh);
    //addDummyVertex(surfMesh);

    // Add a mask of ints to denote how many elements are referencing this vertex
    auto referenceCountPtr = std::make_shared<DataArray<int>>(tissueMesh->getNumVertices());
    referenceCountPtr->fill(0);
    tissueMesh->setVertexAttribute("ReferenceCount", referenceCountPtr);

    // Use FEMTet constraints
    model->getConfig()->m_femParams->m_YoungModulus = 10000.0;
    model->getConfig()->m_femParams->m_PoissonRatio = 0.49;
    model->getConfig()->enableFemConstraint(PbdFemConstraint::MaterialType::StVK);

    // Setup the material
    auto material = std::make_shared<RenderMaterial>();
    material->setBackFaceCulling(false);
    material->setDisplayMode(RenderMaterial::DisplayMode::WireframeSurface);
    material->setShadingModel(RenderMaterial::ShadingModel::PBR);

    // Setup the Object
    auto tissueObj = std::make_shared<PbdObject>(name);
    tissueObj->setPhysicsGeometry(tissueMesh);
    tissueObj->setVisualGeometry(tissueMesh);
    tissueObj->setCollidingGeometry(surfMesh);
    auto map = std::make_shared<PointwiseMap>(tissueMesh, surfMesh);
    tissueObj->setPhysicsToCollidingMap(map);
    tissueObj->getVisualModel(0)->setRenderMaterial(material);
    tissueObj->setDynamicalModel(model);
    tissueObj->getPbdBody()->uniformMassValue = 0.1;
    // Fix the borders
    for (int z = 0; z < dim[2]; z++)
    {
        for (int y = 0; y < dim[1]; y++)
        {
            for (int x = 0; x < dim[0]; x++)
            {
                if (/*x == 0 ||*/ z == 0 || /*x == dim[0] - 1 ||*/ z == dim[2] - 1)
                {
                    tissueObj->getPbdBody()->fixedNodeIds.push_back(x + dim[0] * (y + dim[1] * z) + 1); // +1 for dummy vertex
                }
            }
        }
    }
    tissueObj->getPbdBody()->fixedNodeIds.push_back(0); // Fix dummy vertex

    return tissueObj;
}

static std::shared_ptr<PbdObject>
makeToolObj(std::shared_ptr<PbdModel> model)
{
    auto plane = std::make_shared<Plane>(Vec3d(0.0, 0.0, 0.0), Vec3d(1.0, 0.0, 0.0));
    plane->setWidth(0.005);
    std::shared_ptr<SurfaceMesh> toolGeom = GeometryUtils::toSurfaceMesh(plane);

    auto toolMesh = MeshIO::read<SurfaceMesh>(
        iMSTK_DATA_ROOT "/Surgical Instruments/Scalpel/Scalpel_Hull_Subdivided_Shifted.stl");
    toolMesh->rotate(Vec3d(1.0, 0.0, 0.0), -PI_2, Geometry::TransformType::ApplyToData);
    toolMesh->scale(0.01, Geometry::TransformType::ApplyToData);

    auto toolObj = std::make_shared<PbdObject>("Tool");
    toolObj->setVisualGeometry(toolMesh);
    auto visualModel = std::make_shared<VisualModel>();
    visualModel->setGeometry(plane);
    toolObj->addVisualModel(visualModel);
    toolObj->setCollidingGeometry(toolGeom);
    toolObj->setPhysicsGeometry(plane);
    toolObj->setPhysicsToCollidingMap(std::make_shared<IsometricMap>(plane, toolGeom));
    toolObj->setPhysicsToVisualMap(std::make_shared<IsometricMap>(plane, toolMesh));
    toolObj->setDynamicalModel(model);
    toolObj->getVisualModel(0)->getRenderMaterial()->setColor(Color::Blue);
    toolObj->getVisualModel(0)->getRenderMaterial()->setDisplayMode(RenderMaterial::DisplayMode::WireframeSurface);
    toolObj->getVisualModel(0)->getRenderMaterial()->setBackFaceCulling(false);
    toolObj->getVisualModel(0)->getRenderMaterial()->setLineWidth(1.0);

    toolObj->getPbdBody()->setRigid(
        Vec3d(0.0, 0.8, 0.0),         // Position
        0.2,                          // Mass
        Quatd::Identity(),            // Orientation
        Mat3d::Identity() * 10000.0); // Inertia

    auto controller = toolObj->addComponent<PbdObjectController>();
    controller->setControlledObject(toolObj);
    controller->setLinearKs(1000.0);
    controller->setLinearKd(50.0);
    controller->setAngularKs(10000000.0);
    controller->setAngularKd(500000.0);
    controller->setForceScaling(0.001);

    return toolObj;
}

///
/// \brief This example demonstrates tetrahedral removal of a pbd simulated mesh
/// using a haptic device. Hold the button the device whilst moving it over elements
/// to remove
///
int
main()
{
    // Setup logger (write to file and stdout)
    Logger::startLogger();

    // Setup the scene
    auto scene =  std::make_shared<Scene>("PbdTissueCut");
    scene->getActiveCamera()->setPosition(0.0, 0.1, 0.1);
    scene->getActiveCamera()->setFocalPoint(0.0, 0.0, 0.0);
    scene->getActiveCamera()->setViewUp(-0.01, 0.48, -0.88);

    // Setup the Model/System
    auto pbdModel = std::make_shared<PbdModel>();
    pbdModel->getConfig()->m_doPartitioning = false;
    pbdModel->getConfig()->m_gravity    = Vec3d(0.0, -0.2, 0.0);
    pbdModel->getConfig()->m_dt         = 0.005;
    pbdModel->getConfig()->m_iterations = 2;

    // Setup a tissue
    std::shared_ptr<PbdObject> tissueObj = makeTissueObj("Tissue",
        Vec3d(0.1, 0.03, 0.1), Vec3i(15, 5, 15), Vec3d(0.0, 0.0, 0.0),
        pbdModel);
    scene->addSceneObject(tissueObj);

    std::shared_ptr<PbdObject> toolObj = makeToolObj(pbdModel);
    scene->addSceneObject(toolObj);

    /*auto interaction = std::make_shared<PbdObjectCollision>(
        toolObj, tissueObj, "ClosedSurfaceMeshToMeshCD");
    scene->addInteraction(interaction);*/

    // Light
    auto light = std::make_shared<DirectionalLight>();
    light->setFocalPoint(Vec3d(5.0, -8.0, -5.0));
    light->setIntensity(1.0);
    scene->addLight("Light", light);

    // Run the simulation
    {
        // Setup a viewer to render
        auto viewer = std::make_shared<VTKViewer>();
        viewer->setActiveScene(scene);
        viewer->setVtkLoggerMode(VTKViewer::VTKLoggerMode::MUTE);
        viewer->setDebugAxesLength(0.1, 0.1, 0.1);

        // Setup a scene manager to advance the scene
        auto sceneManager = std::make_shared<SceneManager>();
        sceneManager->setActiveScene(scene);
        sceneManager->pause(); // Start simulation paused

        auto driver = std::make_shared<SimulationManager>();
        driver->addModule(viewer);
        driver->addModule(sceneManager);
        driver->setDesiredDt(0.005);

        // Setup default haptics manager
        std::shared_ptr<DeviceManager> hapticManager = DeviceManagerFactory::makeDeviceManager();
        std::shared_ptr<DeviceClient>  deviceClient  = hapticManager->makeDeviceClient();
        driver->addModule(hapticManager);

        auto controller = toolObj->getComponent<PbdObjectController>();
        controller->setDevice(deviceClient);

        connect<Event>(sceneManager, &SceneManager::preUpdate, [&](Event*)
            {
                // Keep the tool moving in real time
                pbdModel->getConfig()->m_dt = sceneManager->getDt();
            });

        connect<Event>(sceneManager, &SceneManager::postUpdate, [&](Event*)
            {
                if (deviceClient->getButton(0))
                {
                    auto tissueMesh = std::dynamic_pointer_cast<TetrahedralMesh>(tissueObj->getPhysicsGeometry());
                    auto toolGeom   = std::dynamic_pointer_cast<SurfaceMesh>(toolObj->getCollidingGeometry());

                    // Default config of the tool is pointing downwards on y
                    const Mat3d rot     = toolGeom->getRotation();
                    const Vec3d forward = (rot * Vec3d(0.0, 0.0, 1.0)).normalized();
                    const Vec3d left    = (rot * Vec3d(1.0, 0.0, 0.0)).normalized();
                    const Vec3d n       = (rot * Vec3d(0.0, 1.0, 0.0)).normalized();

                    const Vec3d planePos    = toolGeom->getTranslation();
                    const Vec3d planeNormal = n;
                    auto planeGeom = std::dynamic_pointer_cast<Plane>(toolObj->getPhysicsGeometry());
                    const double planeWidth     = planeGeom->getWidth() * 1.1; // Slightly larger than collision geometry
                    const double planeHalfWidth = planeWidth * 0.5;

                    std::shared_ptr<VecDataArray<double, 3>> tissueVerticesPtr = tissueMesh->getVertexPositions();
                    std::shared_ptr<VecDataArray<int, 4>> tissueIndicesPtr     = tissueMesh->getCells();
                    VecDataArray<double, 3>& tissueVertices = *tissueVerticesPtr;
                    VecDataArray<int, 4>& tissueIndices     = *tissueIndicesPtr;

                    // Compute which tets should be removed
                    std::unordered_set<int> removedTets;
                    for (int i = 0; i < tissueIndices.size(); i++)
                    {
                        Vec4i& tet = tissueIndices[i];
                        std::array<Vec3d, 4> tetVerts;
                        tetVerts[0] = tissueVertices[tet[0]];
                        tetVerts[1] = tissueVertices[tet[1]];
                        tetVerts[2] = tissueVertices[tet[2]];
                        tetVerts[3] = tissueVertices[tet[3]];

                        if (splitTest(tetVerts, planePos, left, planeHalfWidth, forward, planeHalfWidth, n))
                        {
                            // Remove the tet being split
                            removedTets.insert(i);
                        }
                    }

                    // Deal with diffs
                    std::shared_ptr<PbdConstraintContainer> constraintsPtr = tissueObj->getPbdModel()->getConstraints();
                    const std::vector<std::shared_ptr<PbdConstraint>>& constraints = constraintsPtr->getConstraints();

                    // First process all removed tets by removing the constraints and setting the element to the dummy vertex
                    for (auto i : removedTets)
                    {
                        Vec4i& tet = tissueIndices[i];

                        // Find and remove the associated constraints
                        for (auto j = constraints.begin(); j != constraints.end(); j++)
                        {
                            const std::vector<PbdParticleId>& vertexIds = (*j)->getParticles();
                            bool isSameTet = true;
                            for (int k = 0; k < 4; k++)
                            {
                                if (vertexIds[k].second != tet[k])
                                {
                                    isSameTet = false;
                                    break;
                                }
                            }
                            if (isSameTet)
                            {
                                constraintsPtr->eraseConstraint(j);
                                break;
                            }
                        }

                        // Set removed tet to dummy vertex
                        tet = Vec4i(0, 0, 0, 0);
                    }

                    if (removedTets.size() > 0)
                    {
                        //// Update collision geometry by re-extracting the entire mesh
                        //auto map = std::dynamic_pointer_cast<PointwiseMap>(tissueObj->getPhysicsToCollidingMap());
                        //std::shared_ptr<SurfaceMesh> colMesh = tissueMesh->extractSurfaceMesh();
                        //map->setChildGeometry(colMesh);
                        //map->compute();
                        //map->update();
                        //colMesh->computeVertexNormals();
                        //tissueObj->setCollidingGeometry(colMesh);

                        tissueIndicesPtr->postModified();
                        tissueMesh->postModified();
                    }
                }
        });

        // Add default mouse and keyboard controls to the viewer
        std::shared_ptr<Entity> mouseAndKeyControls =
            SimulationUtils::createDefaultSceneControl(driver);
        scene->addSceneObject(mouseAndKeyControls);

        driver->start();
    }

    return 0;
}