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
#include "imstkLight.h"
#include "imstkLogger.h"
#include "imstkNew.h"
#include "imstkOpenVRDeviceClient.h"
#include "imstkScene.h"
#include "imstkSceneManager.h"
#include "imstkSimulationManager.h"
#include "imstkVisualObjectImporter.h"
#include "imstkVTKOpenVRViewer.h"
#include "imstkVisualModel.h"
#include "imstkRenderMaterial.h"
#include "imstkDirectionalLight.h"

using namespace imstk;

///
/// \brief This example demonstrates linear skinned animation for hands in 
/// VR. It also demonstrates IK.
///
int
main()
{
    // Setup logger (write to file and stdout)
    Logger::startLogger();

    imstkNew<Scene> scene("SkinnedAnimation");

    auto leftHand = ObjectIO::importAnimObject("leftHand",
        "C:/Users/Andx_/Desktop/hands/LeftHand.fbx",
        "C:/Users/Andx_/Desktop/hands/");
    auto rightHand = ObjectIO::importAnimObject("rightHand",
        "C:/Users/Andx_/Desktop/hands/RightHand.fbx",
        "C:/Users/Andx_/Desktop/hands/");

    // Add in scene
    scene->addSceneObject(leftHand);
    scene->addSceneObject(rightHand);

    // Light (white)
    /*imstkNew<DirectionalLight> whiteLight("whiteLight");
    whiteLight->setDirection(Vec3d(5.0, -8.0, 5.0));
    whiteLight->setIntensity(1.0);
    scene->addLight(whiteLight);*/

    imstkNew<DirectionalLight> whiteLight;
    whiteLight->setDirection(Vec3d(0.0, 0.0, -1.0));
    whiteLight->setIntensity(1.0);
    scene->addLight("light", whiteLight);

    auto leftHandAnimModel = std::dynamic_pointer_cast<AnimationModel>(leftHand->getAnimationModel());
    auto rightHandAnimModel = std::dynamic_pointer_cast<AnimationModel>(rightHand->getAnimationModel());

    // Run the simulation
    {
        // Setup a viewer to render
        imstkNew<VTKOpenVRViewer> viewer("Viewer");
        viewer->setActiveScene(scene);

        // Setup a scene manager to advance the scene in its own thread
        imstkNew<SceneManager> sceneManager("Scene Manager");
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
        const Mat4d rightHandOffset =
            mat4dRotation(Rotd(-rotX, Vec3d(1.0, 0.0, 0.0))) *
            mat4dRotation(Rotd(-rotY, Vec3d(0.0, 1.0, 0.0))) *
            mat4dRotation(Rotd(-rotZ, Vec3d(0.0, 0.0, 1.0))) *
            mat4dTranslate(trans) *
            mat4dScale(Vec3d(0.1, 0.1, 0.1));

        double t = 0.0;
        connect<Event>(viewer, &Viewer::postUpdate,
            [&](Event*)
            {
                const Vec3d posLeft = viewer->getVRDeviceClient(OPENVR_LEFT_CONTROLLER)->getPosition();
                const Quatd quatLeft = viewer->getVRDeviceClient(OPENVR_LEFT_CONTROLLER)->getOrientation();

                const Vec3d posRight = viewer->getVRDeviceClient(OPENVR_RIGHT_CONTROLLER)->getPosition();
                const Quatd quatRight = viewer->getVRDeviceClient(OPENVR_RIGHT_CONTROLLER)->getOrientation();

                leftHandAnimModel->setTransform(mat4dTranslate(posLeft) * mat4dRotation(quatLeft) * leftHandOffset);
                rightHandAnimModel->setTransform(mat4dTranslate(posRight) * mat4dRotation(quatRight) * rightHandOffset);
                leftHand->getAnimationModel()->update();
                rightHand->getAnimationModel()->update();
            });
        connect<ButtonEvent>(viewer->getVRDeviceClient(OPENVR_LEFT_CONTROLLER), &OpenVRDeviceClient::buttonStateChanged,
            [&](ButtonEvent* e)
            {
                if (e->m_button == 3)
                {
                    if (e->m_buttonState == BUTTON_PRESSED)
                    {
                        leftHandAnimModel->stopAnimation("LeftHandArmature|LeftHand_3FingerUp");
                        leftHandAnimModel->playAnimation("LeftHandArmature|LeftHand_3FingerDown");
                    }
                    else if (e->m_buttonState == BUTTON_RELEASED)
                    {
                        leftHandAnimModel->stopAnimation("LeftHandArmature|LeftHand_3FingerDown");
                        leftHandAnimModel->playAnimation("LeftHandArmature|LeftHand_3FingerUp");
                    }
                }
                else if (e->m_button == 0)
                {
                    if (e->m_buttonState == BUTTON_PRESSED)
                    {
                        leftHandAnimModel->playAnimation("LeftHandArmature|LeftHand_IndexDown");
                    }
                    else if (e->m_buttonState == BUTTON_TOUCHED)
                    {
                        leftHandAnimModel->stopAnimation("LeftHandArmature|LeftHand_IndexDown");
                        leftHandAnimModel->playAnimation("LeftHandArmature|LeftHand_IndexTouch");
                    }
                    else if (e->m_buttonState == BUTTON_UNTOUCHED)
                    {
                        leftHandAnimModel->stopAnimation("LeftHandArmature|LeftHand_IndexDown");
                        leftHandAnimModel->stopAnimation("LeftHandArmature|LeftHand_IndexTouch");
                    }
                }
                else if (e->m_button == 4)
                {
                    if (e->m_buttonState == BUTTON_TOUCHED)
                    {
                        leftHandAnimModel->stopAnimation("LeftHandArmature|LeftHand_ThumbUp");
                        leftHandAnimModel->playAnimation("LeftHandArmature|LeftHand_ThumbDown");
                    }
                    else if (e->m_buttonState == BUTTON_UNTOUCHED)
                    {
                        leftHandAnimModel->stopAnimation("LeftHandArmature|LeftHand_ThumbDown");
                        leftHandAnimModel->playAnimation("LeftHandArmature|LeftHand_ThumbUp");
                    }
                }
            });
        connect<ButtonEvent>(viewer->getVRDeviceClient(OPENVR_RIGHT_CONTROLLER), &OpenVRDeviceClient::buttonStateChanged,
            [&](ButtonEvent* e)
            {
                if (e->m_button == 3)
                {
                    if (e->m_buttonState == BUTTON_PRESSED)
                    {
                        rightHandAnimModel->stopAnimation("LeftHandArmature|LeftHand_3FingerUp");
                        rightHandAnimModel->playAnimation("LeftHandArmature|LeftHand_3FingerDown");
                    }
                    else if (e->m_buttonState == BUTTON_RELEASED)
                    {
                        rightHandAnimModel->stopAnimation("LeftHandArmature|LeftHand_3FingerDown");
                        rightHandAnimModel->playAnimation("LeftHandArmature|LeftHand_3FingerUp");
                    }
                }
                else if (e->m_button == 0)
                {
                    if (e->m_buttonState == BUTTON_PRESSED)
                    {
                        rightHandAnimModel->playAnimation("LeftHandArmature|LeftHand_IndexDown");
                    }
                    else if (e->m_buttonState == BUTTON_TOUCHED)
                    {
                        rightHandAnimModel->stopAnimation("LeftHandArmature|LeftHand_IndexDown");
                        rightHandAnimModel->playAnimation("LeftHandArmature|LeftHand_IndexTouch");
                    }
                    else if (e->m_buttonState == BUTTON_UNTOUCHED)
                    {
                        rightHandAnimModel->stopAnimation("LeftHandArmature|LeftHand_IndexDown");
                        rightHandAnimModel->stopAnimation("LeftHandArmature|LeftHand_IndexTouch");
                    }
                }
                else if (e->m_button == 4)
                {
                    if (e->m_buttonState == BUTTON_TOUCHED)
                    {
                        rightHandAnimModel->stopAnimation("LeftHandArmature|LeftHand_ThumbUp");
                        rightHandAnimModel->playAnimation("LeftHandArmature|LeftHand_ThumbDown");
                    }
                    else if (e->m_buttonState == BUTTON_UNTOUCHED)
                    {
                        rightHandAnimModel->stopAnimation("LeftHandArmature|LeftHand_ThumbDown");
                        rightHandAnimModel->playAnimation("LeftHandArmature|LeftHand_ThumbUp");
                    }
                }
            });

        driver->start();
    }

    return 0;
}
