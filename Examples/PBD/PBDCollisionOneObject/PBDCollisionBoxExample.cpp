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
//
//#include "imstkCamera.h"
//#include "imstkCollisionGraph.h"
//#include "imstkKeyboardSceneControl.h"
//#include "imstkLight.h"
//#include "imstkLogger.h"
//#include "imstkMouseSceneControl.h"
//#include "imstkNew.h"
//#include "imstkObjectInteractionFactory.h"
//#include "imstkOneToOneMap.h"
//#include "imstkPbdModel.h"
//#include "imstkPbdObject.h"
//#include "imstkRenderMaterial.h"
//#include "imstkScene.h"
//#include "imstkSceneManager.h"
//#include "imstkSurfaceMesh.h"
//#include "imstkTetrahedralMesh.h"
//#include "imstkTetraTriangleMap.h"
//#include "imstkVisualModel.h"
//#include "imstkVTKViewer.h"
//#include "imstkSphere.h"
//#include "imstkImplicitGeometryToSurfaceMeshCD.h"
//#include "imstkCollisionData.h"
//
//#include "imstkVTKRenderer.h"
//
//using namespace imstk;
//
/////
///// \brief This example demonstrates the collision interaction
///// using Position based dynamics
/////
//int
//main()
//{
//    // Setup logger (write to file and stdout)
//    Logger::startLogger();
//
//    // Setup the scene
//    imstkNew<Scene> scene("PbdCollisionOneDragon");
//
//    scene->getActiveCamera()->setPosition(0, 1, -5.0);
//    scene->getActiveCamera()->setFocalPoint(0.0, 0.0, 0.0);
//
//    imstkNew<CollidingObject> testObj("SphereObj");
//    imstkNew<Sphere> sphere(Vec3d(0.0, 0.0, 0.0), 0.5);
//    testObj->setVisualGeometry(sphere);
//    testObj->setCollidingGeometry(sphere);
//    testObj->getVisualModel(0)->getRenderMaterial()->setOpacity(0.5);
//    scene->addSceneObject(testObj);
//
//    imstkNew<CollidingObject> optSphereObj("OptSphereObj");
//    imstkNew<Sphere> optSphere(Vec3d(0.0, 0.0, 0.0), 0.025);
//    optSphereObj->setVisualGeometry(optSphere);
//    optSphereObj->setCollidingGeometry(optSphere);
//    optSphereObj->getVisualModel(0)->getRenderMaterial()->setColor(Color::Red);
//    scene->addSceneObject(optSphereObj);
//
//    imstkNew<CollidingObject> triangleOptSphereObj("triangleOptsphereObj");
//    imstkNew<Sphere> tOptSphere(Vec3d(0.0, 0.0, 0.0), 0.025);
//    triangleOptSphereObj->setVisualGeometry(tOptSphere);
//    triangleOptSphereObj->setCollidingGeometry(tOptSphere);
//    triangleOptSphereObj->getVisualModel(0)->getRenderMaterial()->setColor(Color::Blue);
//    scene->addSceneObject(triangleOptSphereObj);
//
//    // Single equilaterial triangle (on the XZ plane, lying exactly on the sphere, centroid is nearest point)
//    std::shared_ptr<VecDataArray<double, 3>> verticesPtr = std::make_shared<VecDataArray<double, 3>>(3);
//    VecDataArray<double, 3>& vertices = *verticesPtr;
//   /* vertices[0] = Vec3d(-0.5, 0.9, -1.8);
//    vertices[1] = Vec3d(0.0, 0.5, 0.5);
//    vertices[2] = Vec3d(0.5, 0.5, -0.5);*/
//   /* vertices[0] = Vec3d(-0.5, 0.2, -0.5);
//    vertices[1] = Vec3d(0.0, 0.2, 0.5);
//    vertices[2] = Vec3d(0.5, 0.2, -0.5);*/
//    /*vertices[0] = Vec3d(1.2, 0.0, 0.0);
//    vertices[1] = Vec3d(0.0, -1.0, 0.0);
//    vertices[2] = Vec3d(0.0, 0.0, -1.0);*/
//    /*vertices[0] = Vec3d(-0.5, 0.9, 1.8);
//    vertices[1] = Vec3d(0.0, 0.5, 1.5);
//    vertices[2] = Vec3d(0.5, 0.5, 1.2);*/
//    vertices[0] = Vec3d(-0.5, 5.0, -0.5);
//    vertices[1] = Vec3d(0.0, 5.0, 0.5);
//    vertices[2] = Vec3d(0.5, 5.0, -0.5);
//    std::shared_ptr<VecDataArray<int, 3>> indicesPtr = std::make_shared<VecDataArray<int, 3>>(1);
//    VecDataArray<int, 3>& indices = *indicesPtr;
//    indices[0] = Vec3i(0, 1, 2);
//
//    imstkNew<SurfaceMesh> triangleMesh;
//    triangleMesh->initialize(verticesPtr, indicesPtr);
//
//    imstkNew<RenderMaterial> floorMaterial;
//    floorMaterial->setDisplayMode(RenderMaterial::DisplayMode::WireframeSurface);
//    imstkNew<VisualModel> floorVisualModel(triangleMesh.get());
//    floorVisualModel->setRenderMaterial(floorMaterial);
//
//    imstkNew<CollidingObject> floorObj("Floor");
//    floorObj->setCollidingGeometry(triangleMesh);
//    floorObj->addVisualModel(floorVisualModel);
//
//    scene->addSceneObject(floorObj);
//
//    //imstkNew<CollisionData> colData;
//    //imstkNew<ImplicitGeometryToSurfaceMeshCD> cd(sphere, triangleMesh, colData);
//    //cd->computeCollisionData();
//    //system("pause");
//
//
//    // Run the simulation
//    // Setup a viewer to render in its own thread
//    imstkNew<VTKViewer> viewer("Viewer");
//    viewer->setActiveScene(scene);
//
//    // Setup a scene manager to advance the scene in its own thread
//    imstkNew<SceneManager> sceneManager("Scene Manager");
//    sceneManager->setActiveScene(scene);
//    viewer->addChildThread(sceneManager); // SceneManager will start/stop with viewer
//
//    //{
//    //    ImplicitFunctionCentralGradient centralGrad;
//    //    centralGrad.setFunction(sphere);
//    //    centralGrad.setDx(Vec3d(0.05, 0.05, 0.05));
//    //    //double m_gradEpsilon = 0.01;
//    //    double m_stepLength = 0.01;
//
//    //    const int i1 = indices[0][0];
//    //    const int i2 = indices[0][1];
//    //    const int i3 = indices[0][2];
//
//    //    const Vec3d& v1 = vertices[i1];
//    //    const Vec3d& v2 = vertices[i2];
//    //    const Vec3d& v3 = vertices[i3];
//    //    const Vec3d centroid = (v1 + v2 + v3) / 3.0;
//
//    //    // Find the maximal point from centroid for radius
//    //    const double rSqr1 = (centroid - v1).squaredNorm();
//    //    const double rSqr2 = (centroid - v2).squaredNorm();
//    //    const double rSqr3 = (centroid - v3).squaredNorm();
//    //    double minRSqr = rSqr1;
//    //    if (minRSqr < rSqr2)
//    //    {
//    //        minRSqr = rSqr2;
//    //    }
//    //    if (minRSqr < rSqr3)
//    //    {
//    //        minRSqr = rSqr3;
//    //    }
//    //    const double r = std::sqrt(minRSqr);
//
//    //    // Get distance from sphere/triangle center to surface
//    //    const double distToSurf = sphere->getFunctionValue(centroid);
//
//    //    // First compute the starting point for optimization
//    //    // That would be the closest point to the implicit surface
//    //    const double dist1 = sphere->getFunctionValue(v1);
//    //    const double dist2 = sphere->getFunctionValue(v2);
//    //    const double dist3 = sphere->getFunctionValue(v3);
//
//    //    //Vec3d pos = centroid;
//    //    //double minDist = distToSurf;
//    //    Vec3d pos = vertices[indices[0][0]] + (vertices[indices[0][1]] - vertices[indices[0][0]]) * 0.378;
//    //    double minDist = sphere->getFunctionValue((vertices[indices[0][0]] + vertices[indices[0][1]]) * 0.5);
//    //   /* if (dist1 < minDist)
//    //    {
//    //        minDist = dist1;
//    //        pos = vertices[indices[0][0]];
//    //    }*/
//    //    if (dist2 < minDist)
//    //    {
//    //        minDist = dist2;
//    //        pos = vertices[indices[0][1]];
//    //    }
//    //    if (dist3 < minDist)
//    //    {
//    //        minDist = dist3;
//    //        pos = vertices[indices[0][2]];
//    //    }
//
//    //    Vec3d barycentricCoords = baryCentric(pos, v1, v2, v3); // Current iterate point BC coords
//    //    Vec3d grad = centralGrad(pos);
//    //    Vec3d barycentricGrad = Vec3d(grad.dot(v1), grad.dot(v2), grad.dot(v3));
//
//    //    printf("position:             %f, %f, %f\n", pos[0], pos[1], pos[2]);
//    //    printf("dist:                 %f\n", sphere->getFunctionValue(pos));
//    //    printf("barycentric Coords:   %f, %f, %f\n", barycentricCoords[0], barycentricCoords[1], barycentricCoords[2]);
//    //    printf("gradient:             %f, %f, %f\n", grad[0], grad[1], grad[2]);
//    //    printf("gradient magnitude:   %f\n", grad.norm());
//    //    printf("barycentric gradient: %f, %f, %f\n", barycentricGrad[0], barycentricGrad[1], barycentricGrad[2]);
//    //    printf("bc grad magnitude:    %f\n\n", barycentricGrad.norm());
//
//    //    optSphere->setPosition(pos);
//    //    optSphere->modified();
//    //    tOptSphere->setPosition(closestPointOnTriangle(pos, v1, v2, v3));
//    //    tOptSphere->modified();
//
//    //    int j = 0;
//    //    connect<Event>(sceneManager, EventType::PostUpdate, [&](Event*)
//    //    {
//    //        if (j++ % 100 == 0)
//    //        {
//    //            // Step along the gradient direction wrt barycentric coordinates (eq5)
//    //            barycentricCoords -= barycentricGrad * m_stepLength;
//
//    //            // Update point
//    //            pos = barycentricCoords[0] * v1 + barycentricCoords[1] * v2 + barycentricCoords[2] * v3;
//
//    //            // Project onto triangle
//    //            //if (barycentricCoords[0] < 0.0 || barycentricCoords[2] < 0.0 || (barycentricCoords[1] + barycentricCoords[2]) > 1.0)
//    //            //{
//    //               // pos = closestPointOnTriangle(pos, v1, v2, v3);
//    //            pos[1] = 0.6;
//
//    //                barycentricCoords = baryCentric(pos, v1, v2, v3); // Update bc coords
//    //            //}
//
//    //            //barycentricCoords = baryCentric(pos, v1, v2, v3); // Update bc coords
//
//    //            // Compute new gradient
//    //            grad = centralGrad(pos);
//    //            barycentricGrad = Vec3d(grad.dot(v1), grad.dot(v2), grad.dot(v3));
//    //            //printf("Test:                 %f, %f, %f\n");
//
//    //            printf("position:             %f, %f, %f\n", pos[0], pos[1], pos[2]);
//    //            printf("dist:                 %f\n", sphere->getFunctionValue(pos));
//    //            printf("barycentric Coords:   %f, %f, %f\n", barycentricCoords[0], barycentricCoords[1], barycentricCoords[2]);
//    //            printf("gradient:             %f, %f, %f\n", grad[0], grad[1], grad[2]);
//    //            printf("gradient magnitude:   %f\n", grad.norm());
//    //            printf("barycentric gradient: %f, %f, %f\n", barycentricGrad[0], barycentricGrad[1], barycentricGrad[2]);
//    //            printf("bc grad magnitude:    %f\n\n", barycentricGrad.norm());
//    //        }
//
//    //        optSphere->setPosition(pos);
//    //        optSphere->modified();
//    //        tOptSphere->setPosition(closestPointOnTriangle(pos, v1, v2, v3));
//    //        tOptSphere->modified();
//    //    });
//    //}
//    {
//        ImplicitFunctionCentralGradient centralGrad;
//        centralGrad.setFunction(sphere);
//        centralGrad.setDx(Vec3d(0.05, 0.05, 0.05));
//        //double m_gradEpsilon = 0.01;
//        double m_stepLength = 0.01;
//
//        const int i1 = indices[0][0];
//        const int i2 = indices[0][1];
//        const int i3 = indices[0][2];
//
//        const Vec3d& v1 = vertices[i1];
//        const Vec3d& v2 = vertices[i2];
//        const Vec3d& v3 = vertices[i3];
//
//        const double dist1 = sphere->getFunctionValue(vertices[indices[0][0]]);
//        const double dist2 = sphere->getFunctionValue(vertices[indices[0][1]]);
//        const double dist3 = sphere->getFunctionValue(vertices[indices[0][2]]);
//
//        Vec3d pos = vertices[indices[0][0]];
//        double minDist = dist1;
//        if (dist2 < minDist)
//        {
//            minDist = dist2;
//            pos = vertices[indices[0][1]];
//        }
//        if (dist3 < minDist)
//        {
//            minDist = dist3;
//            pos = vertices[indices[0][2]];
//        }
//
//        optSphere->setPosition(pos);
//        optSphere->modified();
//
//        int j = 0;
//        int k = 0;
//        connect<Event>(sceneManager, EventType::PostUpdate, [&](Event*)
//            {
//                if (j++ % 100 == 0)
//                {
//                    Vec3d grad = centralGrad(pos);
//
//                    // Find support point s
//                    Vec3d s = vertices[indices[0][0]];
//                    double minS = vertices[indices[0][0]].dot(grad);
//                    if (vertices[indices[0][1]].dot(grad) < minS)
//                    {
//                        minS = vertices[indices[0][1]].dot(grad);
//                        s = vertices[indices[0][1]];
//                    }
//                    if (vertices[indices[0][2]].dot(grad) < minS)
//                    {
//                        minS = vertices[indices[0][2]].dot(grad);
//                        s = vertices[indices[0][2]];
//                    }
//
//                    pos += (2.0 / (2.0 + k)) * (s - pos);
//
//                    optSphere->setPosition(pos);
//                    optSphere->modified();
//                    printf("Opt %f\n", minS);
//                    k++;
//                }
//            });
//    }
//
//
//    // Add mouse and keyboard controls to the viewer
//    {
//        imstkNew<MouseSceneControl> mouseControl(viewer->getMouseDevice());
//        mouseControl->setSceneManager(sceneManager);
//        viewer->addControl(mouseControl);
//
//        imstkNew<KeyboardSceneControl> keyControl(viewer->getKeyboardDevice());
//        keyControl->setSceneManager(sceneManager);
//        keyControl->setViewer(viewer);
//        viewer->addControl(keyControl);
//    }
//
//    // Start viewer running, scene as paused
//    sceneManager->requestStatus(ThreadStatus::Paused);
//    std::dynamic_pointer_cast<VTKRenderer>(viewer->getActiveRenderer())->setAxesVisibility(false);
//    viewer->start();
//
//    return 0;
//}

#include "imstkCamera.h"
#include "imstkCollisionGraph.h"
#include "imstkKeyboardSceneControl.h"
#include "imstkLight.h"
#include "imstkLogger.h"
#include "imstkMouseSceneControl.h"
#include "imstkNew.h"
#include "imstkObjectInteractionFactory.h"
#include "imstkOneToOneMap.h"
#include "imstkPbdModel.h"
#include "imstkPbdObject.h"
#include "imstkRenderMaterial.h"
#include "imstkScene.h"
#include "imstkSceneManager.h"
#include "imstkSurfaceMesh.h"
#include "imstkTetrahedralMesh.h"
#include "imstkTetraTriangleMap.h"
#include "imstkVisualModel.h"
#include "imstkVTKViewer.h"
#include "imstkPlane.h"
#include "imstkCollisionPair.h"
#include "imstkImplicitGeometryToSurfaceMeshCD.h"
#include "imstkDebugRenderGeometry.h"
#include "imstkCollisionData.h"
#include "imstkKeyboardDeviceClient.h"
#include "imstkMeshIO.h"
#include "imstkLineMesh.h"
#include "imstkImplicitGeometryToPointSetCD.h"
#include "imstkTaskGraph.h"
#include "imstkPBDCollisionHandling.h"

using namespace imstk;

// parameters to play with
const double youngModulus     = 1000.0;
const double poissonRatio     = 0.3;
const double timeStep         = 0.01;
const double contactStiffness = 0.1;
const int    maxIter = 5;

///
/// \brief Create a surface mesh
/// \param nRows number of vertices in x-direction
/// \param nCols number of vertices in y-direction
///
std::shared_ptr<SurfaceMesh>
createPlaneSurfaceMesh(const double width, const double height, const double z, const size_t nRows, const size_t nCols)
{
    const double dy = width / static_cast<double>(nCols - 1);
    const double dx = height / static_cast<double>(nRows - 1);

    imstkNew<VecDataArray<double, 3>> verticesPtr;
    VecDataArray<double, 3>&          vertices = *verticesPtr.get();
    vertices.resize(nRows * nCols);

    for (size_t i = 0; i < nRows; ++i)
    {
        for (size_t j = 0; j < nCols; j++)
        {
            const double y = static_cast<double>(dy * j);
            const double x = static_cast<double>(dx * i);
            vertices[i * nCols + j] = Vec3d(x - height * 0.5, z, y - width * 0.5);
        }
    }

    // c. Add connectivity data
    imstkNew<VecDataArray<int, 3>> trianglesPtr;
    VecDataArray<int, 3>&          triangles = *trianglesPtr.get();
    for (int i = 0; i < nRows - 1; ++i)
    {
        for (int j = 0; j < nCols - 1; j++)
        {
            triangles.push_back(Vec3i(i * nCols + j, i * nCols + j + 1, (i + 1) * nCols + j));
            triangles.push_back(Vec3i((i + 1) * nCols + j + 1, (i + 1) * nCols + j, i * nCols + j + 1));
        }
    }

    imstkNew<SurfaceMesh> surfMesh;
    surfMesh->initialize(verticesPtr, trianglesPtr);
    return surfMesh;
}

std::shared_ptr<TetrahedralMesh>
createCubeTetMesh(const double width, const double height, const double depth, const Vec3d pos)
{
    const double halfWidth  = width * 0.5;
    const double halfHeight = height * 0.5;
    const double halfDepth  = depth * 0.5;

    imstkNew<TetrahedralMesh> results;

    std::shared_ptr<VecDataArray<double, 3>> verticesPtr = std::make_shared<VecDataArray<double, 3>>(8);
    VecDataArray<double, 3>&                 vertices    = *verticesPtr;
    vertices[0] = Vec3d(-halfWidth, halfHeight, -halfDepth) + pos;
    vertices[1] = Vec3d(halfWidth, halfHeight, -halfDepth) + pos;
    vertices[2] = Vec3d(-halfWidth, halfHeight, halfDepth) + pos;
    vertices[3] = Vec3d(halfWidth, halfHeight, halfDepth) + pos;
    vertices[4] = Vec3d(-halfWidth, -halfHeight, -halfDepth) + pos;
    vertices[5] = Vec3d(halfWidth, -halfHeight, -halfDepth) + pos;
    vertices[6] = Vec3d(-halfWidth, -halfHeight, halfDepth) + pos;
    vertices[7] = Vec3d(halfWidth, -halfHeight, halfDepth) + pos;

    std::shared_ptr<VecDataArray<int, 4>> indicesPtr = std::make_shared<VecDataArray<int, 4>>(5);
    VecDataArray<int, 4>&                 indices    = *indicesPtr;
    indices[0] = Vec4i(4, 1, 2, 7);
    indices[1] = Vec4i(4, 5, 1, 7);
    indices[2] = Vec4i(4, 1, 0, 2);
    indices[3] = Vec4i(1, 7, 3, 2);
    indices[4] = Vec4i(6, 4, 2, 7);

    results->initialize(verticesPtr, indicesPtr);

    return results;
}

//std::shared_ptr<TetrahedralMesh> createCubeTetMesh(const double width, const double height, const double depth, const Vec3d pos)
//{
//    auto tetMesh = MeshIO::read<TetrahedralMesh>(iMSTK_DATA_ROOT "/asianDragon/asianDragon.veg");
//    tetMesh->translate(pos, Geometry::TransformType::ApplyToData);
//    return tetMesh;
//}

///
/// \brief This example demonstrates the collision interaction
/// using Position based dynamics
///
int
main()
{
    // Setup logger (write to file and stdout)
    Logger::startLogger();

    // Setup the scene
    imstkNew<Scene> scene("PbdCollisionOneDragon");

    scene->getActiveCamera()->setPosition(0, 3.0, -50.0);
    scene->getActiveCamera()->setFocalPoint(0.0, 0.0, 0.0);

    // set up the meshes
    std::shared_ptr<TetrahedralMesh> tetMesh = createCubeTetMesh(10.0, 10.0, 10.0, Vec3d(0.0, 6.0, 0.0));
    imstkNew<SurfaceMesh>            surfMesh;
    tetMesh->extractSurfaceMesh(surfMesh, true);
    surfMesh->flipNormals();

    // set up visual model based on high res mesh
    imstkNew<RenderMaterial> material;
    material->setDisplayMode(RenderMaterial::DisplayMode::WireframeSurface);
    material->setLineWidth(0.5);
    material->setColor(Color::Blue);
    material->setShadingModel(RenderMaterial::ShadingModel::Flat);
    material->setOpacity(0.5);
    imstkNew<VisualModel> surfMeshModel(surfMesh.get());
    surfMeshModel->setRenderMaterial(material);

    // configure the deformable object
    imstkNew<PbdObject> deformableObj("DeformableObj");
    deformableObj->addVisualModel(surfMeshModel);
    deformableObj->setCollidingGeometry(surfMesh);
    deformableObj->setPhysicsGeometry(tetMesh);
    deformableObj->setPhysicsToCollidingMap(std::make_shared<OneToOneMap>(tetMesh, surfMesh));

    // Create model and object
    imstkNew<PbdModel> pbdModel;
    pbdModel->setModelGeometry(tetMesh);

    // configure model
    imstkNew<PBDModelConfig> pbdParams;

    // FEM constraint
    pbdParams->m_femParams->m_YoungModulus = youngModulus;
    pbdParams->m_femParams->m_PoissonRatio = poissonRatio;
    pbdParams->enableFEMConstraint(PbdConstraint::Type::FEMTet,
        PbdFEMConstraint::MaterialType::Corotation);

    // Other parameters
    // \todo use lumped mass
    pbdParams->m_uniformMassValue = 1.0;
    pbdParams->m_gravity    = Vec3d(0, -10.0, 0);
    pbdParams->m_defaultDt  = timeStep;
    pbdParams->m_iterations = maxIter;
    pbdParams->m_collisionParams->m_proximity = 0.0;
    pbdParams->m_collisionParams->m_stiffness = 1.0;

    pbdModel->configure(pbdParams);
    deformableObj->setDynamicalModel(pbdModel);

    scene->addSceneObject(deformableObj);

    imstkNew<Plane> floorGeom(Vec3d(0.0, 0.0, 0.0), Vec3d(0.0, 1.0, 0.0));
    floorGeom->setWidth(100.0);

    imstkNew<CollidingObject> floorObj("Floor");
    floorObj->setCollidingGeometry(floorGeom);
    floorObj->setVisualGeometry(floorGeom);
    floorObj->getVisualModel(0)->getRenderMaterial()->setDisplayMode(RenderMaterial::DisplayMode::WireframeSurface);
    floorObj->getVisualModel(0)->getRenderMaterial()->setOpacity(0.5);

    scene->addSceneObject(floorObj);

    // Collision
    std::shared_ptr<CollisionPair> collisionPair = std::dynamic_pointer_cast<CollisionPair>(makeObjectInteractionPair(deformableObj, floorObj,
        InteractionType::PbdObjToCollidingObjCollision, CollisionDetection::Type::SurfaceMeshToImplicit));
    std::shared_ptr<ImplicitGeometryToSurfaceMeshCD> cd = std::dynamic_pointer_cast<ImplicitGeometryToSurfaceMeshCD>(collisionPair->getCollisionDetection());
    cd->setEpsilon(0.01);
    /* std::shared_ptr<ImplicitGeometryToPointSetCD> cd = std::dynamic_pointer_cast<ImplicitGeometryToPointSetCD>(collisionPair->getCollisionDetection());*/
    scene->getCollisionGraph()->addInteraction(collisionPair);

    imstkNew<DebugRenderPoints> debugPts("Debug Points");
    imstkNew<VisualModel>       debugPtsVsModel(debugPts.get());
    debugPtsVsModel->getRenderMaterial()->setDisplayMode(RenderMaterial::DisplayMode::WireframeSurface);
    debugPtsVsModel->getRenderMaterial()->setColor(Color::Yellow);
    debugPtsVsModel->getRenderMaterial()->setPointSize(8.0);
    debugPtsVsModel->getRenderMaterial()->setRenderPointsAsSpheres(true);
    scene->addDebugVisualModel(debugPtsVsModel);

    imstkNew<DebugRenderLines> debugLines("Debug Lines");
    imstkNew<VisualModel>      debugLinesVsModel(debugLines.get());
    debugLinesVsModel->getRenderMaterial()->setDisplayMode(RenderMaterial::DisplayMode::Wireframe);
    debugLinesVsModel->getRenderMaterial()->setColor(Color::Green);
    debugLinesVsModel->getRenderMaterial()->setLineWidth(4.0);
    debugLinesVsModel->getRenderMaterial()->setBackFaceCulling(false);
    scene->addDebugVisualModel(debugLinesVsModel);

//#define DEBUG

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
        sceneManager->init();

        connect<KeyPressEvent>(viewer->getKeyboardDevice(), EventType::KeyEvent, [&](KeyPressEvent* e)
        {
            if (e->m_keyPressType == KEY_PRESS && e->m_key == ' ')
            {
                sceneManager->update();

                VecDataArray<double, 3>& vertices    = *std::dynamic_pointer_cast<PointSet>(deformableObj->getPhysicsGeometry())->getVertexPositions();
                VecDataArray<int, 4>& indices        = *std::dynamic_pointer_cast<TetrahedralMesh>(deformableObj->getPhysicsGeometry())->getTetrahedraIndices();
                VecDataArray<int, 3>& surfIndices    = *surfMesh->getTriangleIndices();
                std::shared_ptr<GeometryMap> geomMap = deformableObj->getPhysicsToCollidingMap();

                debugPts->clear();
                debugLines->clear();
                for (int i = 0; i < cd->getCollisionData()->PDColData.getSize(); i++)
                {
                    const Vec3d& p1 = vertices[cd->getCollisionData()->PDColData[i].nodeIdx];
                    const Vec3d p2  = p1 + cd->getCollisionData()->PDColData[i].penetrationDepth * cd->getCollisionData()->PDColData[i].dirAtoB;

                    //debugPts->appendVertex(p1);

                    debugLines->appendVertex(p1);
                    debugLines->appendVertex(p2);
                }
                debugPts->setDataModified(true);
                debugLines->setDataModified(true);
            }
            });

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

        connect<Event>(scene, EventType::Configure, [&](Event*)
            {
                scene->getTaskGraph()->insertAfter(cd->getTaskNode(), std::make_shared<TaskNode>([&]()
                    {
                        VecDataArray<double, 3>& vertices = *std::dynamic_pointer_cast<PointSet>(deformableObj->getPhysicsGeometry())->getVertexPositions();
                        VecDataArray<int, 3>& indices = *std::dynamic_pointer_cast<SurfaceMesh>(deformableObj->getCollidingGeometry())->getTriangleIndices();

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
                            const Vec3d dir = element.dir;

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

                            //const Vec3d n = v0.cross(v1).normalized();
                            //const double l = v2.dot(dir);

                            //const double dist = m_configB->m_proximity;
                            /*if (l > dist)
                            {
                                c = 0.0;
                                return false;
                            }*/

                            Vec3d dcdxB[3];
                            dcdxB[0] = -u * dir;
                            dcdxB[1] = -v * dir;
                            dcdxB[2] = -w * dir;

                            double c = element.closestDistance;
                            double lambda = c / (dcdxB[0].squaredNorm() + dcdxB[1].squaredNorm() + dcdxB[2].squaredNorm());

                            Vec3d dx[3];
                            dx[0] = -lambda * dcdxB[0];
                            dx[1] = -lambda * dcdxB[1];
                            dx[2] = -lambda * dcdxB[2];

                            debugLines->appendVertex(x1);
                            //debugLines->appendVertex(x1 + -u * n);
                            debugLines->appendVertex(x1 + dx[0]);

                            debugLines->appendVertex(x2);
                            //debugLines->appendVertex(x2 + -v * n);
                            debugLines->appendVertex(x2 + dx[1]);

                            debugLines->appendVertex(x3);
                            //debugLines->appendVertex(x3 + -w * n);
                            debugLines->appendVertex(x3 + dx[2]);

                            /*for (int i = 0; i < 3; i++)
                            {
                                if (dx[i].norm() > 0.3)
                                {
                                    std::cout << "id: " << triVerts[i] << std::endl;
                                }
                            }*/
                        }
                        debugPts->setDataModified(true);
                        debugLines->setDataModified(true);
                        //std::cout << "Iter" << std::endl;
                    }, "test"));
            
                /*scene->getTaskGraph()->insertAfter(collisionPair->getCollisionHandlingA()->getTaskNode(), std::make_shared<TaskNode>([&]()
                    {
                        MeshIO::write(std::dynamic_pointer_cast<TetrahedralMesh>(deformableObj->getPhysicsGeometry()), "C:/Users/Andx_/Desktop/test.vtk");
                    }, "testa"));*/
            });

        // Add mouse and keyboard controls to the viewer
        {
            imstkNew<MouseSceneControl> mouseControl(viewer->getMouseDevice());
            mouseControl->setSceneManager(sceneManager);
            viewer->addControl(mouseControl);

            connect<KeyPressEvent>(viewer->getKeyboardDevice(), EventType::KeyEvent, [&](KeyPressEvent* e)
                {
                    if (e->m_key == 'h' && e->m_keyPressType == KEY_PRESS)
                    {
                        scene->advance();
                    }
                });

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
