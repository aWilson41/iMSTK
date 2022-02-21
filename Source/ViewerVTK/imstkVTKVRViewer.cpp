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

#include "imstkVTKVRViewer.h"
#include "imstkCamera.h"
#include "imstkDeviceControl.h"
#include "imstkLogger.h"
#include "imstkScene.h"
#include "imstkVRDeviceClient.h"
#include "imstkVTKInteractorStyleVR.h"
#include "imstkVTKRenderer.h"

#include <vtkMatrix4x4.h>
#include <vtkRenderer.h>

// \todo: Should be removed upon upgrade of VTK. OPENXR preferred
#ifdef iMSTK_USE_OPENXR
#include <vtkOpenXRRenderWindow.h>
#include <vtkOpenXRRenderWindowInteractor.h>
#include <vtkVRModel.h>
#include <vtkOpenGLState.h>

using vtkImstkVRModel = vtkVRModel;
using vtkImstkVRRenderWindow = vtkOpenXRRenderWindow;
using vtkImstkVRRenderWindowInteractor = vtkOpenXRRenderWindowInteractor;
#else
#include <vtkOpenVRModel.h>
#include <vtkOpenVRRenderWindow.h>
#include "imstkVtkOpenVRRenderWindowInteractorImstk.h"

using vtkImstkVRModel = vtkOpenVRModel;
using vtkImstkVRRenderWindow = vtkOpenVRRenderWindow;
using vtkImstkVRRenderWindowInteractor = vtkOpenVRRenderWindowInteractorImstk;
#endif

namespace imstk
{
VTKVRViewer::VTKVRViewer(std::string name) : AbstractVTKViewer(name)
{
    // Create the interactor style
    auto vrInteractorStyle = vtkSmartPointer<vtkInteractorStyleVR>::New();
    m_vtkInteractorStyle = vrInteractorStyle;

    // Create the interactor
    auto iren = vtkSmartPointer<vtkImstkVRRenderWindowInteractor>::New();
    iren->SetInteractorStyle(m_vtkInteractorStyle);

    // Create the RenderWindow
    m_vtkRenderWindow = vtkSmartPointer<vtkImstkVRRenderWindow>::New();
    m_vtkRenderWindow->SetInteractor(iren);
    iren->SetRenderWindow(m_vtkRenderWindow);
    m_vtkRenderWindow->HideCursor();

    m_vrDeviceClients.push_back(vrInteractorStyle->getLeftControllerDeviceClient());
    m_vrDeviceClients.push_back(vrInteractorStyle->getRightControllerDeviceClient());
    m_vrDeviceClients.push_back(vrInteractorStyle->getHmdDeviceClient());
}

void
VTKVRViewer::setActiveScene(std::shared_ptr<Scene> scene)
{
    // If already current scene
    if (scene == m_activeScene)
    {
        LOG(WARNING) << scene->getName() << " already is the viewer current scene.";
        return;
    }

    // If the current scene has a renderer, remove it
    if (m_activeScene)
    {
        auto vtkRenderer = std::dynamic_pointer_cast<VTKRenderer>(this->getActiveRenderer())->getVtkRenderer();
        if (m_vtkRenderWindow->HasRenderer(vtkRenderer))
        {
            m_vtkRenderWindow->RemoveRenderer(vtkRenderer);
        }
    }

    // Update current scene
    m_activeScene = scene;

    // Create renderer if it doesn't exist
    if (!m_rendererMap.count(m_activeScene))
    {
        m_rendererMap[m_activeScene] = std::make_shared<VTKRenderer>(m_activeScene, true);
    }

    // Cast to VTK renderer
    auto vtkRenderer = std::dynamic_pointer_cast<VTKRenderer>(this->getActiveRenderer())->getVtkRenderer();

    // Set renderer to renderWindow
    m_vtkRenderWindow->AddRenderer(vtkRenderer);

    m_vtkInteractorStyle->SetCurrentRenderer(vtkRenderer);
}

void
VTKVRViewer::setPhysicalToWorldTransform(const Mat4d& physicalToWorldMatrix)
{
    auto                 renWin = vtkImstkVRRenderWindow::SafeDownCast(m_vtkRenderWindow);
    vtkNew<vtkMatrix4x4> mat;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            mat->SetElement(i, j, physicalToWorldMatrix(i, j));
        }
    }
    renWin->SetPhysicalToWorldMatrix(mat);
}

Mat4d
VTKVRViewer::getPhysicalToWorldTransform()
{
    auto                 renWin = vtkImstkVRRenderWindow::SafeDownCast(m_vtkRenderWindow);
    Mat4d                transform;
    vtkNew<vtkMatrix4x4> mat;
    renWin->GetPhysicalToWorldMatrix(mat);
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            transform(i, j) = mat->GetElement(i, j);
        }
    }
    return transform;
}

void
VTKVRViewer::setRenderingMode(const Renderer::Mode mode)
{
    if (!m_activeScene)
    {
        LOG(WARNING) << "Missing scene, can not set rendering mode.\n"
                     << "Use Viewer::setCurrentScene to setup scene.";
        return;
    }

    // Setup renderer
    this->getActiveRenderer()->setMode(mode, true);
}

void
VTKVRViewer::processEvents()
{
    // Custom call to only process input events, do not perform a render
    auto iren = vtkImstkVRRenderWindowInteractor::SafeDownCast(m_vtkRenderWindow->GetInteractor());
    auto ren  = std::dynamic_pointer_cast<imstk::VTKRenderer>(getActiveRenderer());
    iren->DoOneEvent(vtkImstkVRRenderWindow::SafeDownCast(m_vtkRenderWindow), ren->getVtkRenderer(), false);

    // Update all controls
    for (auto control : m_controls)
    {
        control->update(m_dt);
    }
}

bool
VTKVRViewer::initModule()
{
    if (!AbstractVTKViewer::initModule())
    {
        return false;
    }

    // VR interactor doesn't support timers, here we throw timer event every update
    // another option would be to conform VTKs VR interactor
    auto iren = vtkImstkVRRenderWindowInteractor::SafeDownCast(m_vtkRenderWindow->GetInteractor());
    //iren->Start(); // Cannot use
    if (iren->HasObserver(vtkCommand::StartEvent))
    {
        iren->InvokeEvent(vtkCommand::StartEvent, nullptr);
        return true;
    }

#ifdef iMSTK_USE_OPENXR
    // Actions must be added after initialization of interactor
    vtkInteractorStyleVR* iStyle = vtkInteractorStyleVR::SafeDownCast(m_vtkInteractorStyle);
    iStyle->addButtonActions();
    iStyle->addMovementActions();
    iStyle->addHapticAction();

    iren->Initialize();
#else
    renWin->Initialize();
    iren->Initialize();

    // Hide the device overlays
    // \todo: Display devices in debug mode
    // Must do one render to initialize vtkOpenVRModel's to then hide the devices
    renWin->Render();

    // Actions must be added after initialization of interactor for openvr
    vtkInteractorStyleVR* iStyle = vtkInteractorStyleVR::SafeDownCast(m_vtkInteractorStyle);
    iStyle->addButtonActions();
    iStyle->addMovementActions();
    setControllerVisibility(false);
#endif

    return true;
}

void
VTKVRViewer::updateModule()
{
    auto ren = std::dynamic_pointer_cast<imstk::VTKRenderer>(getActiveRenderer());
    if (ren == nullptr)
    {
        return;
    }

    // For the VR view we can't supply the a camera in the normal sense
    // we need to pre multiply a "user view"
    std::shared_ptr<Camera> cam  = getActiveScene()->getActiveCamera();
    const Mat4d&            view = cam->getView();
    setPhysicalToWorldTransform(view);

#ifdef iMSTK_USE_OPENXR
    auto renWin = vtkImstkVRRenderWindow::SafeDownCast(m_vtkRenderWindow);

    auto ostate = renWin->GetState();
    renWin->MakeCurrent();
    ostate->Reset();
    ostate->Push();

    if (vtkOpenXRManager::GetInstance()->IsSessionRunning())
    {
        vtkOpenXRManager* xrManager = vtkOpenXRManager::GetInstance();
        if (!xrManager->WaitAndBeginFrame())
        {
            return;
        }

        if (renWin->GetTrackHMD())
        {
            renWin->UpdateHMDMatrixPose();
        }

        if (xrManager->GetShouldRenderCurrentFrame())
        {
            // Call visual update on every scene object
            getActiveScene()->updateVisuals();
            // Update all the rendering delegates
            ren->updateRenderDelegates();

            // Start rendering
            renWin->Superclass::Render();

            // OpenXR does not begin until after the begin event in its event loop
            // Yet the models can't be hidden until first render
            if (!m_didFirstRender)
            {
                m_didFirstRender = true;
                setControllerVisibility(false);
            }
        }

        xrManager->EndFrame();
    }

    ostate->Pop();
#else
    // Call visual update on every scene object
    getActiveScene()->updateVisuals();
    // Update all the rendering delegates
    ren->updateRenderDelegates();

    // Render
    //m_vtkRenderWindow->GetInteractor()->Render();
    m_vtkRenderWindow->Render();
#endif
}

std::shared_ptr<VRDeviceClient>
VTKVRViewer::getVRDeviceClient(int deviceType)
{
    auto iter = std::find_if(m_vrDeviceClients.begin(), m_vrDeviceClients.end(),
        [&](const std::shared_ptr<VRDeviceClient>& deviceClient)
        {
            return static_cast<int>(deviceClient->getDeviceType()) == deviceType;
        });
    return (iter == m_vrDeviceClients.end()) ? nullptr : *iter;
}

void
VTKVRViewer::setControllerVisibility(const bool visible)
{
    auto renWin = vtkImstkVRRenderWindow::SafeDownCast(m_vtkRenderWindow);

    if (renWin == nullptr)
    {
        LOG(FATAL) << "Tried to set controller visibility before render window was initialized";
        return;
    }

    // Hide all controller models
#ifdef iMSTK_USE_OPENXR
    for (uint32_t hand :
         { vtkOpenXRManager::ControllerIndex::Left, vtkOpenXRManager::ControllerIndex::Right })
    {
        vtkVRModel* trackedDeviceModel = renWin->GetModelForDeviceHandle(hand);
        if (trackedDeviceModel != nullptr)
        {
            trackedDeviceModel->SetVisibility(visible);
        }
    }
#else
    for (uint32_t i = 0; i < vr::k_unMaxTrackedDeviceCount; i++)
    {
        vtkVRModel* trackedDeviceModel = renWin->GetTrackedDeviceModel(i);
        if (trackedDeviceModel != nullptr)
        {
            trackedDeviceModel->SetVisibility(false);
        }
    }
#endif
}
} // namespace imstk
