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

#include "CameraVRControl.h"
#include "imstkCamera.h"
#include "imstkDirectionalLight.h"
#include "imstkLogger.h"
#include "imstkMeshIO.h"
#include "imstkNew.h"
#include "imstkRenderMaterial.h"
#include "imstkScene.h"
#include "imstkSceneManager.h"
#include "imstkSceneObject.h"
#include "imstkSceneObjectController.h"
#include "imstkSimulationManager.h"
#include "imstkSurfaceMesh.h"
#include "imstkVisualModel.h"
#include "imstkVisualObjectImporter.h"
#include "imstkVRDeviceClient.h"
#include "imstkVTKVRViewer.h"

using namespace imstk;

std::shared_ptr<SceneObject>
makeHandleObject()
{
    imstkNew<SceneObject> scalpelHandle("ScalpelHandle");
    auto                  toolHandleMesh = MeshIO::read<SurfaceMesh>(iMSTK_DATA_ROOT "/Surgical Instruments/Scalpel/Scalpel_Handle.dae");
    toolHandleMesh->translate(0.0, 0.0, 1.0, Geometry::TransformType::ApplyToData);
    toolHandleMesh->rotate(Vec3d(0.0, 1.0, 0.0), 3.14, Geometry::TransformType::ApplyToData);
    toolHandleMesh->rotate(Vec3d(1.0, 0.0, 0.0), -1.57, Geometry::TransformType::ApplyToData);
    toolHandleMesh->scale(0.06, Geometry::TransformType::ApplyToData);
    toolHandleMesh->computeVertexNormals(); // Recompute these as we have transformed the vertices

    scalpelHandle->setVisualGeometry(toolHandleMesh);

    std::shared_ptr<RenderMaterial> material = scalpelHandle->getVisualModel(0)->getRenderMaterial();
    material->setDisplayMode(RenderMaterial::DisplayMode::Surface);
    material->setShadingModel(RenderMaterial::ShadingModel::PBR);
    material->setMetalness(0.9);
    material->setRoughness(0.2);
    material->addTexture(std::make_shared<Texture>(
                iMSTK_DATA_ROOT "/Surgical Instruments/Scalpel/Scalpel_Albedo.png",
                Texture::Type::Diffuse));
    material->setIsDynamicMesh(false);

    return scalpelHandle;
}

std::shared_ptr<SceneObject>
makeBlade(std::string filename)
{
    imstkNew<SceneObject> scalpelBlade(filename);
    auto                  blade10Mesh = MeshIO::read<SurfaceMesh>(iMSTK_DATA_ROOT "/Surgical Instruments/Scalpel/" + filename + ".dae");
    blade10Mesh->translate(0.0, 0.0, 1.0, Geometry::TransformType::ApplyToData);
    blade10Mesh->rotate(Vec3d(0.0, 1.0, 0.0), 3.14, Geometry::TransformType::ApplyToData);
    blade10Mesh->rotate(Vec3d(1.0, 0.0, 0.0), -1.57, Geometry::TransformType::ApplyToData);
    blade10Mesh->scale(0.06, Geometry::TransformType::ApplyToData);
    blade10Mesh->computeVertexNormals(); // Recompute these as we have transformed the vertices

    scalpelBlade->setVisualGeometry(blade10Mesh);

    std::shared_ptr<RenderMaterial> material = scalpelBlade->getVisualModel(0)->getRenderMaterial();
    material->setDisplayMode(RenderMaterial::DisplayMode::Surface);
    material->setShadingModel(RenderMaterial::ShadingModel::PBR);
    material->setMetalness(0.9);
    material->setRoughness(0.2);
    material->addTexture(std::make_shared<Texture>(
                iMSTK_DATA_ROOT "/Surgical Instruments/Scalpel/Scalpel_Albedo.png",
                Texture::Type::Diffuse));
    material->setIsDynamicMesh(false);

    return scalpelBlade;
}

///
/// \brief This example demonstrates rendering and controlling a SceneObject with OpenVR
/// as well as swapping a tool
///
int
main()
{
    // Write log to stdout and file
    Logger::startLogger();

    // Setup the scene
    imstkNew<Scene> scene("VRControllerExample");

    std::shared_ptr<SceneObject> scalpelHandle = makeHandleObject();
    scene->addSceneObject(scalpelHandle);

    std::shared_ptr<SceneObject> scalpelBlade10 = makeBlade("Scalpel_Blade10");
    scene->addSceneObject(scalpelBlade10);

    std::shared_ptr<SceneObject> scalpelBlade15 = makeBlade("Scalpel_Blade15");
    scene->addSceneObject(scalpelBlade15);
    scalpelBlade15->getVisualGeometry()->setTranslation(0.2, 1.0, -0.8);

    std::shared_ptr<SceneObject> tableObj = ObjectIO::importSceneObject("Instrument Table",
        iMSTK_DATA_ROOT "/Surgical instruments/Instrument Table/Instrument_Table.dae",
        iMSTK_DATA_ROOT "/Surgical instruments/Instrument Table/");
    scene->addSceneObject(tableObj);

    // Lights
    imstkNew<DirectionalLight> dirLight;
    dirLight->setIntensity(4);
    dirLight->setColor(Color(1.0, 0.95, 0.8));
    scene->addLight("dirlight", dirLight);

    {
        // Add a module to run the viewer
        imstkNew<VTKVRViewer> viewer;
        viewer->setActiveScene(scene);

        // Add a module to run the scene
        imstkNew<SceneManager> sceneManager;
        sceneManager->setActiveScene(scene);

        imstkNew<SimulationManager> driver;
        driver->addModule(viewer);
        driver->addModule(sceneManager);
        driver->setDesiredDt(0.01); // Spend less time updating & more time rendering

        // Add a VR controller for the scalpel handle
        imstkNew<SceneObjectController> controller1(scalpelHandle, viewer->getVRDeviceClient(RIGHT_CONTROLLER));
        scene->addController(controller1);
        // Add a VR controller for the scalpel blade
        imstkNew<SceneObjectController> controller2(scalpelBlade10, viewer->getVRDeviceClient(RIGHT_CONTROLLER));
        scene->addController(controller2);

        imstkNew<CameraVRControl> camControl;
        camControl->setRotateDevice(viewer->getVRDeviceClient(RIGHT_CONTROLLER));
        camControl->setTranslateDevice(viewer->getVRDeviceClient(LEFT_CONTROLLER));
        camControl->setTranslateSpeedScale(1.0);
        camControl->setRotateSpeedScale(1.0);
        camControl->setCamera(scene->getActiveCamera());
        viewer->addControl(camControl); // Only needs to update every render

        bool                            blade10InHand = true;
        std::shared_ptr<VRDeviceClient> rightControllerDevice = viewer->getVRDeviceClient(RIGHT_CONTROLLER);
        connect<ButtonEvent>(rightControllerDevice, &VRDeviceClient::buttonStateChanged,
            [&](ButtonEvent* e)
            {
                // When any button pressed, swap blade
                if (e->m_buttonState == BUTTON_PRESSED)
                {
                    //rightControllerDevice->applyVibration(0.5f, 300000000.0f, 3000.0f);
                    rightControllerDevice->applyVibration(0.5f, -1.0f, 0.0f);

                    const Vec3d& posControl = rightControllerDevice->getPosition();
                    if (blade10InHand)
                    {
                        // Swap to blade 15 only if it's close in space
                        Vec3d min, max;
                        scalpelBlade15->getVisualGeometry()->computeBoundingBox(min, max);
                        const Vec3d posBlade = (min + max) * 0.5;
                        const double dist    = (posControl - posBlade).norm();
                        LOG(INFO) << "Dist: " << dist;
                        if (dist < 2.0)
                        {
                            const Vec3d t = scalpelBlade15->getVisualGeometry()->getTranslation();
                            const Mat3d r = scalpelBlade15->getVisualGeometry()->getRotation();

                            // Set the new blade to move
                            controller2->setControlledSceneObject(scalpelBlade15);
                            blade10InHand = false;

                            scalpelBlade10->getVisualGeometry()->setTranslation(t);
                            scalpelBlade10->getVisualGeometry()->setRotation(r);
                        }
                    }
                    else
                    {
                        // Swap to blade 10 only if it's close in space
                        Vec3d min, max;
                        scalpelBlade10->getVisualGeometry()->computeBoundingBox(min, max);
                        const Vec3d posBlade = (min + max) * 0.5;
                        const double dist    = (posControl - posBlade).norm();
                        LOG(INFO) << "Dist: " << dist;
                        if (dist < 2.0)
                        {
                            const Vec3d t = scalpelBlade10->getVisualGeometry()->getTranslation();
                            const Mat3d r = scalpelBlade10->getVisualGeometry()->getRotation();

                            controller2->setControlledSceneObject(scalpelBlade10);
                            blade10InHand = true;

                            // Swap transforms of the blades
                            scalpelBlade15->getVisualGeometry()->setTranslation(t);
                            scalpelBlade15->getVisualGeometry()->setRotation(r);
                        }
                    }
                }
        });

        driver->start();
    }

    return 0;
}