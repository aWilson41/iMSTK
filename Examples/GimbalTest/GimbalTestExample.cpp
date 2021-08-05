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
#include "imstkCylinder.h"
#include "imstkKeyboardSceneControl.h"
#include "imstkDirectionalLight.h"
#include "imstkLogger.h"
#include "imstkMouseSceneControl.h"
#include "imstkNew.h"
#include "imstkOrientedBox.h"
#include "imstkPlane.h"
#include "imstkSurfaceMesh.h"
#include "imstkMeshIO.h"
#include "imstkRenderMaterial.h"
#include "imstkScene.h"
#include "imstkSceneManager.h"
#include "imstkSceneObject.h"
#include "imstkSimulationManager.h"
#include "imstkVisualModel.h"
#include "imstkVTKViewer.h"
#include "imstkOrientedBox.h"
#include "imstkGimbalControl.h"
#include "imstkKeyboardDeviceClient.h"

using namespace imstk;

///
/// \brief This example demonstrates the geometry transforms in imstk
///
int
main()
{
    // Setup logger (write to file and stdout)
    Logger::startLogger();

    imstkNew<Scene>       scene("GeometryTransforms");
    imstkNew<SceneObject> geomObj("bunny");
    {
        geomObj->setVisualGeometry(MeshIO::read<SurfaceMesh>(iMSTK_DATA_ROOT"stanfordBunny/stanfordBunny.stl"));
        geomObj->getVisualGeometry()->scale(0.1, Geometry::TransformType::ApplyToData);
        geomObj->getVisualGeometry()->translate(-Vec3d(1.0, 1.0, 1.0), Geometry::TransformType::ApplyToData);
        geomObj->getVisualModel(0)->getRenderMaterial()->setOpacity(0.5);
        geomObj->getVisualModel(0)->getRenderMaterial()->setColor(Color::Red);
        scene->addSceneObject(geomObj);
    }

    // Set Camera configuration
    scene->getActiveCamera()->setPosition(Vec3d(0.0, 10.0, 10.0));

    // Light
    imstkNew<DirectionalLight> light;
    light->setDirection(Vec3d(5.0, -8.0, -5.0));
    light->setIntensity(1.0);
    scene->addLight("light", light);

    // Run the simulation
    {
        // Setup a viewer to render in its own thread
        imstkNew<VTKViewer> viewer("Viewer");
        viewer->setActiveScene(scene);
        viewer->setDebugAxesLength(0.0, 0.0, 0.0);

        // Setup a scene manager to advance the scene in its own thread
        imstkNew<SceneManager> sceneManager("Scene Manager");
        sceneManager->setExecutionType(Module::ExecutionType::ADAPTIVE);
        sceneManager->setActiveScene(scene);

        imstkNew<SimulationManager> driver;
        driver->addModule(viewer);
        driver->addModule(sceneManager);

        // Add mouse and keyboard controls to the viewer
        {
            imstkNew<MouseSceneControl> mouseControl1(viewer->getMouseDevice());
            mouseControl1->setSceneManager(sceneManager);
            viewer->addControl(mouseControl1);

            imstkNew<GimbalControl> gizmoControl(viewer->getMouseDevice());
            gizmoControl->setControlledGeometry(geomObj->getVisualGeometry());
            gizmoControl->setViewer(viewer);
            scene->addSceneObject(gizmoControl->getGimbalPanObject());
            //scene->addSceneObject(mouseControl->getGimbalScaleObject());
            //scene->addSceneObject(mouseControl->getGimbalRotateObject());
            viewer->addControl(gizmoControl);

            /*imstkNew<KeyboardSceneControl> keyControl(viewer->getKeyboardDevice());
            keyControl->setSceneManager(sceneManager);
            keyControl->setModuleDriver(driver);
            viewer->addControl(keyControl);*/

            connect<KeyEvent>(viewer->getKeyboardDevice(), &KeyboardDeviceClient::keyPress, [&](KeyEvent* e)
                {
                    if (e->m_key == 'w')
                    {
                        gizmoControl->setTransformMode(GimbalControl::TransformMode::PAN);
                    }
                    else if (e->m_key == 'e')
                    {
                        gizmoControl->setTransformMode(GimbalControl::TransformMode::ROTATE);
                    }
                    else if (e->m_key == 'r')
                    {
                        gizmoControl->setTransformMode(GimbalControl::TransformMode::SCALE);
                    }
                });
        }

        driver->start();
    }

    return 0;
}
