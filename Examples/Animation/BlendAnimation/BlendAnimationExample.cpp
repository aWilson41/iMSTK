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

#include "imstkAnimatedObject.h"
#include "imstkAnimationModel.h"
#include "imstkCamera.h"
#include "imstkKeyboardDeviceClient.h"
#include "imstkKeyboardSceneControl.h"
#include "imstkPointLight.h"
#include "imstkLogger.h"
#include "imstkMouseDeviceClient.h"
#include "imstkMouseSceneControl.h"
#include "imstkNew.h"
#include "imstkScene.h"
#include "imstkSceneManager.h"
#include "imstkSimulationManager.h"
#include "imstkVisualObjectImporter.h"
#include "imstkVTKRenderer.h"
#include "imstkVTKViewer.h"
#include "imstkVisualModel.h"
#include "imstkRenderMaterial.h"
#include "imstkLineMesh.h"
#include "imstkDirectionalLight.h"
#include "imstkMeshIO.h"

using namespace imstk;

///
/// \brief This example demonstrates linear skinned animation blending for hands.
/// Press the buttons 1, 2, & 3 to move each finger. Each plays a separate 
/// animation from the file that only keys animation on the bones of that finger.
/// This way each finger can be played separately.
///
int
main()
{
    // Setup logger (write to file and stdout)
    Logger::startLogger();

    // Setup the scene
    imstkNew<Scene> scene("BlendAnimation");
    // Use our own debug camera
    scene->getConfig()->debugCamBoundingBox = false;
    scene->getCamera("debug")->setFocalPoint(Vec3d(0.0, 0.0, 0.0));
    scene->getCamera("debug")->setPosition(Vec3d(0.0, 1.0, 1.0));

    // Import our animated object
    auto leftHand = ObjectIO::importAnimObject("leftHand",
        "C:/Users/Andrew/Desktop/hands/LeftHand.fbx",
        "C:/Users/Andrew/Desktop/hands/");
    leftHand->getVisualModel(0)->setIsVisible(false);
    leftHand->getVisualModel(0)->getRenderMaterial()->setShadingModel(RenderMaterial::ShadingModel::PBR);
    scene->addSceneObject(leftHand);

    // Light (white)
    imstkNew<DirectionalLight> dirLight;
    dirLight->setDirection(Vec3d(0.0, 0.0, -1.0).normalized());
    dirLight->setIntensity(1.0);
    scene->addLight("dirLight", dirLight);

    /*imstkNew<PointLight> pointLight;
    pointLight->setPosition(Vec3d(2.0, 0.0, 0.0));
    pointLight->setIntensity(1.0);
    scene->addLight("pointLight", pointLight);*/

    /*imstkNew<DirectionalLight> whiteLight2("whiteLight2");
    whiteLight2->setDirection(Vec3d(0.0, 0.0, -1.0));
    whiteLight2->setIntensity(1.0);
    scene->addLight(whiteLight2);*/

    scene->getActiveCamera()->setFocalPoint(Vec3d(0.0, 0.0, 0.0));
    scene->getActiveCamera()->setPosition(Vec3d(0.0, 0.0, 1.2));

    auto leftHandAnimModel =
        std::dynamic_pointer_cast<AnimationModel>(leftHand->getAnimationModel());
    leftHandAnimModel->setGenDebugSkeleton(true);
    leftHand->getVisualModel(0)->getRenderMaterial()->setOpacity(0.2); // Render the hand as slighty transparent

    // Add the debug skeleton
    std::shared_ptr<LineMesh>       debugGeom = leftHandAnimModel->getDebugSkeleton();
    imstkNew<VisualModel>           skeletonVisualModel;
    skeletonVisualModel->setGeometry(debugGeom);
    {
        auto material = skeletonVisualModel->getRenderMaterial();
        material->setBackFaceCulling(false);
        material->setColor(Color::Red);
        material->setLineWidth(2.0);
        material->setPointSize(6.0);
        material->setDisplayMode(RenderMaterial::DisplayMode::WireframeSurface);
    }
    leftHand->addVisualModel(skeletonVisualModel);

    // Set no animation (bind pose)
    //leftHandAnimModel->playAnimation("DefaultPose");

    // Run the simulation
    {
        // Setup a viewer to render
        imstkNew<VTKViewer> viewer;
        viewer->setActiveScene(scene);
        viewer->setDebugAxesLength(0.1, 0.1, 0.1);

        // Setup a scene manager to advance the scene in its own thread
        imstkNew<SceneManager> sceneManager;
        sceneManager->setActiveScene(scene);
        sceneManager->pause();

        imstkNew<SimulationManager> driver;
        driver->addModule(viewer);
        driver->addModule(sceneManager);

        double rotX = PI * 2.0;
        double rotY = PI;
        double rotZ = 5.117;
        Vec3d  trans = Vec3d(0.0, 0.015, -0.165);

        const Mat4d leftHandOffset =
            mat4dRotation(Rotd(rotX, Vec3d(1.0, 0.0, 0.0))) *
            mat4dRotation(Rotd(rotY, Vec3d(0.0, 1.0, 0.0))) *
            mat4dRotation(Rotd(rotZ, Vec3d(0.0, 0.0, 1.0))) *
            mat4dTranslate(trans) *
            mat4dScale(Vec3d(0.1, 0.1, 0.1));

        connect<Event>(viewer, &Viewer::postUpdate,
            [&](Event*)
            {
                const Vec2d& mousePos = viewer->getMouseDevice()->getPos();
                const Vec3d posLeft = Vec3d(mousePos[0] - 0.5, mousePos[1] - 0.5, 0.0);
                const Quatd quatLeft = Quatd(Rotd(PI * 0.5, Vec3d(0.0, 1.0, 0.0)));// viewer->getVRDeviceClient(OPENVR_LEFT_CONTROLLER)->getOrientation();

                leftHandAnimModel->setTransform(mat4dTranslate(posLeft) * mat4dRotation(quatLeft) * leftHandOffset);
                leftHand->getAnimationModel()->update();
            });
        connect<KeyEvent>(viewer->getKeyboardDevice(), &KeyboardDeviceClient::keyPress,
            [&](KeyEvent* e)
            {
                if (e->m_key == '1')
                {
                    leftHandAnimModel->stopAnimation("LeftHandArmature|LeftHand_3FingerUp");
                    leftHandAnimModel->playAnimation("LeftHandArmature|LeftHand_3FingerDown");
                }
                else if (e->m_key == '2')
                {
                    leftHandAnimModel->playAnimation("LeftHandArmature|LeftHand_IndexDown");
                }
                else if (e->m_key == '3')
                {
                    leftHandAnimModel->stopAnimation("LeftHandArmature|LeftHand_ThumbUp");
                    leftHandAnimModel->playAnimation("LeftHandArmature|LeftHand_ThumbDown");
                }
                else if (e->m_key == '4')
                {
                    MeshIO::write(debugGeom, "C:/Users/Andrew/Desktop/test1.vtk");
                    MeshIO::write(std::dynamic_pointer_cast<PointSet>(leftHand->getVisualGeometry()), "C:/Users/Andrew/Desktop/test2.vtk");
                }
            });
        connect<KeyEvent>(viewer->getKeyboardDevice(), &KeyboardDeviceClient::keyRelease,
            [&](KeyEvent* e)
            {
                if (e->m_key == '1')
                {
                    if (e->m_keyPressType == KEY_RELEASE)
                    {
                        leftHandAnimModel->stopAnimation("LeftHandArmature|LeftHand_3FingerDown");
                        leftHandAnimModel->playAnimation("LeftHandArmature|LeftHand_3FingerUp");
                    }
                }
                else if (e->m_key == '2')
                {
                    if (e->m_keyPressType == KEY_RELEASE)
                    {
                        leftHandAnimModel->stopAnimation("LeftHandArmature|LeftHand_IndexDown");
                    }
                }
                else if (e->m_key == '3')
                {
                    if (e->m_keyPressType == KEY_RELEASE)
                    {
                        leftHandAnimModel->stopAnimation("LeftHandArmature|LeftHand_ThumbDown");
                        leftHandAnimModel->playAnimation("LeftHandArmature|LeftHand_ThumbUp");
                    }
                }
            });

        // Add mouse and keyboard controls to the viewer
        {
            imstkNew<MouseSceneControl> mouseControl(viewer->getMouseDevice());
            mouseControl->setSceneManager(sceneManager);
            viewer->addControl(mouseControl);

            imstkNew<KeyboardSceneControl> keyControl(viewer->getKeyboardDevice());
            keyControl->setSceneManager(sceneManager);
            keyControl->setModuleDriver(driver);
            viewer->addControl(keyControl);
        }

        driver->start();
    }

    return 0;
}
