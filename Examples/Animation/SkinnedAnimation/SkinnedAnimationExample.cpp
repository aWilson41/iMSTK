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
#include "imstkAnimationNode.h"
#include "imstkCamera.h"
#include "imstkDirectionalLight.h"
#include "imstkKeyboardDeviceClient.h"
#include "imstkKeyboardSceneControl.h"
#include "imstkLight.h"
#include "imstkLineMesh.h"
#include "imstkLogger.h"
#include "imstkMouseSceneControl.h"
#include "imstkNew.h"
#include "imstkOpenVRDeviceClient.h"
#include "imstkPlane.h"
#include "imstkRenderMaterial.h"
#include "imstkScene.h"
#include "imstkSceneManager.h"
#include "imstkSceneObject.h"
#include "imstkSimulationManager.h"
#include "imstkAnimationModel.h"
#include "imstkSphere.h"
#include "imstkViewer.h"
#include "imstkVisualModel.h"
#include "imstkVisualObjectImporter.h"
#include "imstkVTKOpenVRViewer.h"
#include "imstkVTKRenderer.h"
#include "imstkVTKScreenCaptureUtility.h"
#include "imstkVTKViewer.h"
#include "imstkIKNodeAnimation.h"

using namespace imstk;

///
/// \brief This example demonstrates linear skinned animation with
/// a simple IK using a target position
///
int
main()
{
    // Setup logger (write to file and stdout)
    Logger::startLogger();

    imstkNew<Scene> scene("SkinnedAnimation");

    auto armObj = ObjectIO::importAnimObject("Arm",
        "C:/Users/Andx_/Documents/Blender/Arm Test/arm.fbx",
        "C:/Users/Andx_/Documents/Blender/Arm Test/");

    // Add in scene
    scene->addSceneObject(armObj);

    // Light (white)
    imstkNew<DirectionalLight> whiteLight;
    whiteLight->setDirection(Vec3d(5.0, -8.0, 5.0));
    whiteLight->setIntensity(1.0);
    scene->addLight("light", whiteLight);

    // Target position
    Vec3d targetPos = Vec3d(1.0, 0.0, 0.0);
    imstkNew<Sphere> sphere(targetPos, 0.1);
    imstkNew<SceneObject> sphereObj("TargetPosSphere");
    sphereObj->setVisualGeometry(sphere);
    scene->addSceneObject(sphereObj);

    // Update Camera
    scene->getActiveCamera()->setPosition(0.72, 2.01, 0.65);
    scene->getActiveCamera()->setFocalPoint(0.77, 1.48, 0.18);
    scene->getActiveCamera()->setViewUp(-0.04, 0.67, -0.74);

    std::shared_ptr<AnimationModel> armAnimModel  =
        std::dynamic_pointer_cast<AnimationModel>(armObj->getAnimationModel());
    armAnimModel->setGenDebugSkeleton(true);
    armObj->getVisualModel(0)->getRenderMaterial()->setOpacity(0.5);
    std::shared_ptr<LineMesh>       debugGeom = armAnimModel->getDebugSkeleton();
    imstkNew<VisualModel>           skeletonVisualModel(debugGeom);
    skeletonVisualModel->getRenderMaterial()->setDisplayMode(RenderMaterial::DisplayMode::Surface);
    skeletonVisualModel->getRenderMaterial()->setLineWidth(5.0);
    armObj->addVisualModel(skeletonVisualModel);

    // Setup IK
    {
        // This arm has 3 bones. Two for the arm, one for the end effector.
        // We will setup our target position
        imstkNew<InverseJacobianIKSystem> ikSystem;
        ikSystem->setTargetPos(&targetPos);
        ikSystem->setNumberOfIterations(5);

        // Setup two IK animations for the two bones
        imstkNew<IKNodeAnimation> bone1IKAnim;
        bone1IKAnim->setNode(armAnimModel->getNode("Bone1"));
        imstkNew<IKNodeAnimation> bone2IKAnim;
        bone2IKAnim->setNode(armAnimModel->getNode("Bone2"));

        // Add animations to the system
        ikSystem->addNodeAnimation(bone1IKAnim);
        ikSystem->addNodeAnimation(bone2IKAnim);

        // Then setup a track to play the animations
        imstkNew<AnimationTrack> track;
        track->channels.push_back(bone1IKAnim);
        track->channels.push_back(bone2IKAnim);

        armAnimModel->addAnimation("IK", track);

        // Play IK for the duration of the program
        armAnimModel->playAnimation("IK");
    }

    // Run the simulation
    {
        // Setup a viewer to render
        imstkNew<VTKViewer> viewer("Viewer");
        viewer->setActiveScene(scene);

        // Setup a scene manager to advance the scene in its own thread
        imstkNew<SceneManager> sceneManager("Scene Manager");
        sceneManager->setActiveScene(scene);
        sceneManager->pause();

        imstkNew<SimulationManager> driver;
        driver->addModule(viewer);
        driver->addModule(sceneManager);

        connect<Event>(viewer, &Viewer::postUpdate,
            [&](Event*)
            {
                armObj->getAnimationModel()->update();
            });
        connect<KeyEvent>(viewer->getKeyboardDevice(), &KeyboardDeviceClient::keyPress,
            [&](KeyEvent* e)
            {
                if (e->m_key == 'a' && e->m_keyPressType == KEY_PRESS)
                {
                    targetPos += Vec3d(0.1, 0.0, 0.0);
                    sphere->setPosition(targetPos);
                    sphere->modified();
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
