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
#include "imstkKeyboardSceneControl.h"
#include "imstkLight.h"
#include "imstkLogger.h"
#include "imstkMouseSceneControl.h"
#include "imstkNew.h"
#include "imstkPbdModel.h"
#include "imstkPbdObject.h"
#include "imstkRenderMaterial.h"
#include "imstkScene.h"
#include "imstkSceneManager.h"
#include "imstkSurfaceMesh.h"
#include "imstkVisualModel.h"
#include "imstkVTKViewer.h"
#include "imstkCone.h"
#include "imstkCollisionPair.h"
#include "imstkImplicitGeometryToSurfaceMeshCD.h"
#include "imstkObjectInteractionFactory.h"
#include "imstkCollisionGraph.h"
#include "imstkPlane.h"
#include "imstkSphere.h"
#include "imstkKeyboardDeviceClient.h"
#include "imstkDebugRenderGeometry.h"
#include "imstkCollisionData.h"
#include "imstkPBDCollisionHandling.h"
#include "imstkPbdSolver.h"
#include "imstkTaskGraph.h"

#define DEBUG

using namespace imstk;

// Parameters to play with
const double width  = 5.0;
const double height = 5.0;
const int    nRows  = 5;
const int    nCols  = 5;

///
/// \brief Creates cloth geometry
///
static std::shared_ptr<SurfaceMesh>
makeClothGeometry(const double width,
                  const double depth,
                  const int    nRows,
                  const int    nCols)
{
    imstkNew<SurfaceMesh> clothMesh;

    imstkNew<VecDataArray<double, 3>> verticesPtr(nRows * nCols);
    VecDataArray<double, 3>& vertices = *verticesPtr.get();
    imstkNew<VecDataArray<double, 3>> velocitiesPtr(nRows * nCols);
    VecDataArray<double, 3>& velocities = *velocitiesPtr.get();
    const double                      dy = width / static_cast<double>(nCols - 1);
    const double                      dx = height / static_cast<double>(nRows - 1);
    for (int i = 0; i < nRows; ++i)
    {
        for (int j = 0; j < nCols; j++)
        {
            vertices[i * nCols + j] = Vec3d(dx * static_cast<double>(i), 2.0, dy * static_cast<double>(j)) - Vec3d(width * 0.5, 0.0, height * 0.5);
            velocities[i * nCols + j] = Vec3d(0.0, -5.0, 0.0);
        }
    }

    // Add connectivity data
    imstkNew<VecDataArray<int, 3>> indicesPtr;
    VecDataArray<int, 3>& indices = *indicesPtr.get();
    for (int i = 0; i < nRows - 1; ++i)
    {
        for (int j = 0; j < nCols - 1; j++)
        {
            const int index1 = i * nCols + j;
            const int index2 = index1 + nCols;
            const int index3 = index1 + 1;
            const int index4 = index2 + 1;

            // Interleave [/][\]
            if (i % 2 ^ j % 2)
            {
                indices.push_back(Vec3i(index1, index2, index3));
                indices.push_back(Vec3i(index4, index3, index2));
            }
            else
            {
                indices.push_back(Vec3i(index2, index4, index1));
                indices.push_back(Vec3i(index4, index3, index1));
            }
        }
    }

    clothMesh->initialize(verticesPtr, indicesPtr);
    clothMesh->setVertexAttribute("Velocities", velocitiesPtr);

    return clothMesh;
}

///
/// \brief Creates cloth object
///
static std::shared_ptr<PbdObject>
makeClothObj(const std::string& name,
             const double       width,
             const double       height,
             const int          nRows,
             const int          nCols)
{
    imstkNew<PbdObject> clothObj(name);

    // Setup the Geometry
    std::shared_ptr<SurfaceMesh> clothMesh = makeClothGeometry(width, height, nRows, nCols);
    //clothMesh->rotate(Vec3d(0.0, 0.0, 1.0), 1.0, Geometry::TransformType::ApplyToData);
    //clothMesh->translate(0.0, 2.0, 0.0, Geometry::TransformType::ApplyToData);

    // Setup the Parameters
    imstkNew<PBDModelConfig> pbdParams;
    pbdParams->enableConstraint(PbdConstraint::Type::Distance, 1.0e2);
    //pbdParams->enableConstraint(PbdConstraint::Type::Dihedral, 1.0e1);
    pbdParams->m_uniformMassValue = width * height / ((double)nRows * (double)nCols);
    pbdParams->m_gravity    = Vec3d(0.0, -9.8, 0.0);
    pbdParams->m_defaultDt  = 0.0005;
    pbdParams->m_iterations = 5;
    pbdParams->m_collisionParams->m_proximity = 0.0;
    pbdParams->m_collisionParams->m_stiffness = 1.0;

    // Setup the Model
    imstkNew<PbdModel> pbdModel;
    pbdModel->setModelGeometry(clothMesh);
    pbdModel->configure(pbdParams);

    // Setup the VisualModel
    imstkNew<RenderMaterial> material;
    material->setBackFaceCulling(false);
    material->setDisplayMode(RenderMaterial::DisplayMode::WireframeSurface);
    material->setEdgeColor(Color::Red);

    imstkNew<VisualModel> visualModel(clothMesh);
    visualModel->setRenderMaterial(material);

    // Setup the Object
    clothObj->addVisualModel(visualModel);
    clothObj->setPhysicsGeometry(clothMesh);
    clothObj->setCollidingGeometry(clothMesh);
    clothObj->setDynamicalModel(pbdModel);

    return clothObj;
}

static std::shared_ptr<PbdObject>
makeSinglePbdTriangle(const std::string& name)
{
    const double dt = 0.0005;

    imstkNew<PbdObject> clothObj(name);

    // Setup the geometry
    std::shared_ptr<SurfaceMesh> clothMesh = std::make_shared<SurfaceMesh>();
    std::shared_ptr<VecDataArray<double, 3>> verticesPtr = std::make_shared<VecDataArray<double, 3>>(3);
    VecDataArray<double, 3>& vertices = *verticesPtr.get();
    vertices[0] = Vec3d(0.0, 2.0, 1.078);
    vertices[1] = Vec3d(-1.2, 2.0, -1.0);
    vertices[2] = Vec3d(1.2, 2.0, -1.0);
    std::shared_ptr<VecDataArray<int, 3>> indicesPtr = std::make_shared<VecDataArray<int, 3>>(1);
    VecDataArray<int, 3>& indices = *indicesPtr.get();
    indices[0] = Vec3i(0, 1, 2);

    clothMesh->initialize(verticesPtr, indicesPtr);

    std::shared_ptr<VecDataArray<double, 3>> velocitiesPtr = std::make_shared<VecDataArray<double, 3>>(3);
    VecDataArray<double, 3>& velocities = *velocitiesPtr;
    velocities[0] = Vec3d(0.0, -0.05, 0.0) / dt;
    velocities[1] = Vec3d(0.0, -0.05, 0.0) / dt;
    velocities[2] = Vec3d(0.0, -0.05, 0.0) / dt; // Displace into the cone 0.2m
    //clothMesh->setVertexAttribute("Velocities", velocitiesPtr);

    // Setup the Parameters
    imstkNew<PBDModelConfig> pbdParams;
    pbdParams->enableConstraint(PbdConstraint::Type::Distance, 1.0e2);
    //pbdParams->enableConstraint(PbdConstraint::Type::Dihedral, 1.0e1);
    pbdParams->m_uniformMassValue = 1.0;
    pbdParams->m_defaultDt  = dt;
    pbdParams->m_iterations = 10;
    pbdParams->m_collisionParams->m_proximity = 0.01;
    pbdParams->m_collisionParams->m_stiffness = 0.01;

    // Setup the Model
    imstkNew<PbdModel> pbdModel;
    pbdModel->setModelGeometry(clothMesh);
    pbdModel->configure(pbdParams);

    // Setup the VisualModel
    imstkNew<RenderMaterial> material;
    material->setBackFaceCulling(false);
    material->setDisplayMode(RenderMaterial::DisplayMode::WireframeSurface);
    material->setEdgeColor(Color::Red);

    imstkNew<VisualModel> visualModel(clothMesh);
    visualModel->setRenderMaterial(material);

    // Setup the Object
    clothObj->addVisualModel(visualModel);
    clothObj->setPhysicsGeometry(clothMesh);
    clothObj->setCollidingGeometry(clothMesh);
    clothObj->setDynamicalModel(pbdModel);

    return clothObj;
}

static std::shared_ptr<CollidingObject>
makeConeObj()
{
    imstkNew<CollidingObject> obj("cone");

    imstkNew<Cone> cone(Vec3d(0.1, -2.0, 0.0), 4.0);
    obj->setCollidingGeometry(cone);
    obj->setVisualGeometry(cone);
    obj->getVisualModel(0)->getRenderMaterial()->setDisplayMode(RenderMaterial::DisplayMode::WireframeSurface);
    obj->getVisualModel(0)->getRenderMaterial()->setOpacity(0.5);

    return obj;
}

///
/// \brief This example demonstrates the cloth simulation
/// using Position based dynamics
///
int
main()
{
    // Write log to stdout and file
    Logger::startLogger();

    // Setup a scene
    imstkNew<Scene> scene("PBDClothCollision");

    //std::shared_ptr<PbdObject> clothObj = makeClothObj("Cloth", width, height, nRows, nCols);
    std::shared_ptr<PbdObject> clothObj = makeSinglePbdTriangle("Cloth");
    scene->addSceneObject(clothObj);

    std::shared_ptr<CollidingObject> coneObj = makeConeObj();
    scene->addSceneObject(coneObj);

    // Collision with cone
    std::shared_ptr<CollisionPair> collisionPair = std::dynamic_pointer_cast<CollisionPair>(makeObjectInteractionPair(clothObj, coneObj,
        InteractionType::PbdObjToCollidingObjCollision, CollisionDetection::Type::SurfaceMeshToImplicit));
    std::shared_ptr<ImplicitGeometryToSurfaceMeshCD> cd = std::dynamic_pointer_cast<ImplicitGeometryToSurfaceMeshCD>(collisionPair->getCollisionDetection());
    cd->setEpsilon(0.01);
    /* std::shared_ptr<ImplicitGeometryToPointSetCD> cd = std::dynamic_pointer_cast<ImplicitGeometryToPointSetCD>(collisionPair->getCollisionDetection());*/
    scene->getCollisionGraph()->addInteraction(collisionPair);


    //imstkNew<Plane> floorGeom(Vec3d(0.0, 0.0, 0.0), Vec3d(0.0, 1.0, 0.0));
    //floorGeom->setWidth(100.0);
    //imstkNew<Sphere> floorGeom(Vec3d(0.0, -1.0, 0.0), 2.0);

    //imstkNew<CollidingObject> floorObj("Floor");
    //floorObj->setCollidingGeometry(floorGeom);
    //floorObj->setVisualGeometry(floorGeom);
    //floorObj->getVisualModel(0)->getRenderMaterial()->setDisplayMode(RenderMaterial::DisplayMode::WireframeSurface);
    //floorObj->getVisualModel(0)->getRenderMaterial()->setOpacity(0.5);

    //scene->addSceneObject(floorObj);

    //// Collision with plane
    //scene->getCollisionGraph()->addInteraction(makeObjectInteractionPair(clothObj, floorObj,
    //    InteractionType::PbdObjToCollidingObjCollision, CollisionDetection::Type::SurfaceMeshToImplicit));

    // Light (white)
    imstkNew<DirectionalLight> whiteLight("whiteLight");
    whiteLight->setFocalPoint(Vec3d(5.0, -8.0, -5.0));
    whiteLight->setIntensity(1.0);
    scene->addLight(whiteLight);

    // Adjust camera
    scene->getActiveCamera()->setFocalPoint(0.0, 0.0, 0.0);
    scene->getActiveCamera()->setPosition(0.0, -5.0, 25.0);

    imstkNew<DebugRenderPoints> debugPts("Debug Points");
    imstkNew<VisualModel>       debugPtsVsModel(debugPts.get());
    debugPtsVsModel->getRenderMaterial()->setDisplayMode(RenderMaterial::DisplayMode::WireframeSurface);
    debugPtsVsModel->getRenderMaterial()->setColor(Color::Yellow);
    debugPtsVsModel->getRenderMaterial()->setPointSize(20.0);
    debugPtsVsModel->getRenderMaterial()->setRenderPointsAsSpheres(true);
    scene->addDebugVisualModel(debugPtsVsModel);

    imstkNew<DebugRenderLines> debugLines("Debug Lines");
    imstkNew<VisualModel>      debugLinesVsModel(debugLines.get());
    debugLinesVsModel->getRenderMaterial()->setDisplayMode(RenderMaterial::DisplayMode::Wireframe);
    debugLinesVsModel->getRenderMaterial()->setColor(Color::Green);
    debugLinesVsModel->getRenderMaterial()->setLineWidth(4.0);
    debugLinesVsModel->getRenderMaterial()->setBackFaceCulling(false);
    scene->addDebugVisualModel(debugLinesVsModel);

#ifdef DEBUG
    // Run the simulation
    {
        // Setup a viewer to render in its own thread
        imstkNew<VTKViewer> viewer("Viewer");
        viewer->setActiveScene(scene);

        // Setup a scene manager to advance the scene in its own thread
        imstkNew<SceneManager> sceneManager("Scene Manager");
        sceneManager->setActiveScene(scene);
        //viewer->addChildThread(sceneManager); // SceneManager will start/stop with viewer
        

        //bool test = false;

        connect<Event>(scene, EventType::Configure, [&](Event*)
            {
                scene->getTaskGraph()->insertAfter(cd->getTaskNode(), std::make_shared<TaskNode>([&]()
                    {
                        VecDataArray<double, 3>& vertices = *std::dynamic_pointer_cast<PointSet>(clothObj->getPhysicsGeometry())->getVertexPositions();
                        VecDataArray<int, 3>& indices = *std::dynamic_pointer_cast<SurfaceMesh>(clothObj->getPhysicsGeometry())->getTriangleIndices();

                        //debugPts->clear();
                        //debugLines->clear();
                        //for (int i = 0; i < cd->getCollisionData()->TFVColData.getSize(); i++)
                        //{
                        //    //const Vec3d& p1 = vertices[cd->getCollisionData()->PDColData[i].nodeIdx];
                        //    //const Vec3d p2 = p1 + cd->getCollisionData()->PDColData[i].penetrationDepth * cd->getCollisionData()->PDColData[i].dirAtoB;
                        //    const Vec3d v = cd->getCollisionData()->TFVColData[i].vertexPt;

                        //    debugPts->appendVertex(v);

                        //    //debugLines->appendVertex(p1);
                        //    //debugLines->appendVertex(p2);
                        //}
                        //debugPts->setDataModified(true);
                        ////debugLines->setDataModified(true);

                        debugPts->clear();
                        debugLines->clear();
                        for (int i = 0; i < cd->getCollisionData()->TFVColData.getSize(); i++)
                        {
                            TriangleFixedVertexCollisionDataElement element = cd->getCollisionData()->TFVColData[i];
                            const Vec3d x0 = element.vertexPt;
                            const Vec3i& triVerts = indices[element.triIdx];
                            const Vec3d x1 = vertices[triVerts[0]];
                            const Vec3d x2 = vertices[triVerts[1]];
                            const Vec3d x3 = vertices[triVerts[2]];

                            debugPts->appendVertex(x0);

                            Vec3d        v0 = x2 - x1;
                            Vec3d        v1 = x3 - x1;
                            Vec3d        v2 = x0 - x1;
                            const double d00 = v0.dot(v0);
                            const double d01 = v0.dot(v1);
                            const double d11 = v1.dot(v1);
                            const double d20 = v2.dot(v0);
                            const double d21 = v2.dot(v1);
                            const double denom = d00 * d11 - d01 * d01;
                            const double v = (d11 * d20 - d01 * d21) / denom;
                            const double w = (d00 * d21 - d01 * d20) / denom;
                            const double u = 1.0 - v - w;

                            const Vec3d n = v0.cross(v1).normalized();
                            const double l = v2.dot(n);

                            //const double dist = m_configB->m_proximity;
                            /*if (l > dist)
                            {
                                c = 0.0;
                                return false;
                            }*/

                            debugLines->appendVertex(x1);
                            debugLines->appendVertex(x1 + -u * n);

                            debugLines->appendVertex(x2);
                            debugLines->appendVertex(x2 + -v * n);

                            debugLines->appendVertex(x3);
                            debugLines->appendVertex(x3 + -w * n);

                            //dcdxA[0] = Vec3d::Zero();

                            // Weight n out over the 3 verts (u,v,w sum to 1)
                            //dcdxB[0] = -u * n;
                            //dcdxB[1] = -v * n;
                            //dcdxB[2] = -w * n;

                            //c = l - dist;
                        }
                        debugPts->setDataModified(true);
                        debugLines->setDataModified(true);
                    }, "test"));

                /*scene->getTaskGraph()->insertAfter(clothObj->getPbdModel()->getTaskGraph()->getSink(), std::make_shared<TaskNode>([&]()
                    {
                        VecDataArray<double, 3>& velocities = *std::dynamic_pointer_cast<VecDataArray<double, 3>>(std::dynamic_pointer_cast<PointSet>(clothObj->getPhysicsGeometry())->getVertexAttribute("Velocities"));
                        if (velocities[0].norm() > 0.9)
                        {
                            test = true;
                        }
                    }, "hmm"));*/
            });

        sceneManager->init();

        connect<KeyPressEvent>(viewer->getKeyboardDevice(), EventType::KeyEvent, [&](KeyPressEvent* e)
            {
                if (e->m_keyPressType == KEY_PRESS && e->m_key == 'h')
                {
                    sceneManager->update();
                    printf("advance\n");
                }
            });
        /*connect<Event>(viewer, EventType::PostUpdate, [&](Event*)
            {
                if (test == false)
                {
                    sceneManager->update();
                }
            });*/

        // Add mouse and keyboard controls to the viewer
        {
            imstkNew<MouseSceneControl> mouseControl(viewer->getMouseDevice());
            mouseControl->setSceneManager(sceneManager);
            viewer->addControl(mouseControl);

            imstkNew<KeyboardSceneControl> keyControl(viewer->getKeyboardDevice());
            keyControl->setSceneManager(sceneManager);
            keyControl->setViewer(viewer);
            viewer->addControl(keyControl);
        }

        // Start viewer running, scene as paused
        //sceneManager->requestStatus(ThreadStatus::Paused);
        viewer->start();
    }
#else
    // Run the simulation
    {
        // Setup a viewer to render in its own thread
        imstkNew<VTKViewer> viewer("Viewer");
        viewer->setActiveScene(scene);

        // Setup a scene manager to advance the scene in its own thread
        imstkNew<SceneManager> sceneManager("Scene Manager");
        sceneManager->setActiveScene(scene);
        viewer->addChildThread(sceneManager); // SceneManager will start/stop with viewer

        // Add mouse and keyboard controls to the viewer
        {
            imstkNew<MouseSceneControl> mouseControl(viewer->getMouseDevice());
            mouseControl->setSceneManager(sceneManager);
            viewer->addControl(mouseControl);

            imstkNew<KeyboardSceneControl> keyControl(viewer->getKeyboardDevice());
            keyControl->setSceneManager(sceneManager);
            keyControl->setViewer(viewer);
            viewer->addControl(keyControl);
        }

        // Start viewer running, scene as paused
        sceneManager->requestStatus(ThreadStatus::Paused);
        viewer->start();
    }
#endif

    return 0;
}
