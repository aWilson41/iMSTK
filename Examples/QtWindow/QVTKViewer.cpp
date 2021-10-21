///*=========================================================================
//
//   Library: iMSTK
//
//   Copyright (c) Kitware, Inc. & Center for Modeling, Simulation,
//   & Imaging in Medicine, Rensselaer Polytechnic Institute.
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0.txt
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.
//
//=========================================================================*/

#include "QVTKViewer.h"
#include "imstkLogger.h"
#include "imstkScene.h"
#include "imstkVTKInteractorStyle.h"
#include "imstkVTKRenderer.h"

#include <QVTKOpenGLNativeWidget.h>
#include <vtkCamera.h>
#include <vtkOpenGLRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkTextActor.h>

#ifdef WIN32
#include <vtkWin32HardwareWindow.h>
#include <vtkWin32RenderWindowInteractor.h>
#else
#include "imstkVtkXRenderWindowInteractor2.h"
#endif

#include <qapplication.h>

using namespace imstk;

QVTKViewer::QVTKViewer(QVTKOpenGLNativeWidget* widget, std::string name) : AbstractVTKViewer(name)
{
    m_vtkRenderWindow = widget->renderWindow();

    // Create the interactor style
    m_interactorStyle    = std::make_shared<VTKInteractorStyle>();
    m_vtkInteractorStyle = std::dynamic_pointer_cast<vtkInteractorStyle>(m_interactorStyle);

    // Create the interactor
#ifdef WIN32
    vtkNew<vtkRenderWindowInteractor> iren;
    iren->SetInteractorStyle(m_vtkInteractorStyle.get());
#else
    vtkSmartPointer<vtkXRenderWindowInteractor2> iren = vtkSmartPointer<vtkXRenderWindowInteractor2>::New();
    iren->SetInteractorStyle(m_vtkInteractorStyle.get());
#endif

    // Create the RenderWindow
    m_vtkRenderWindow->SetInteractor(iren);
}

void
QVTKViewer::setActiveScene(std::shared_ptr<Scene> scene)
{
    // If already current scene
    if (scene == m_activeScene)
    {
        LOG(WARNING) << scene->getName() << " already is the viewer current scene.";
        return;
    }

    // If the current scene has a renderer, remove it
    if (m_activeScene != nullptr)
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
        m_rendererMap[m_activeScene] = std::make_shared<VTKRenderer>(m_activeScene, false);
    }

    // Cast to VTK renderer
    auto vtkRenderer = std::dynamic_pointer_cast<VTKRenderer>(this->getActiveRenderer())->getVtkRenderer();

    // Set renderer to renderWindow
    m_vtkRenderWindow->AddRenderer(vtkRenderer);
    m_vtkInteractorStyle->SetCurrentRenderer(vtkRenderer);

    // Update the camera
    std::shared_ptr<VTKRenderer> ren = std::dynamic_pointer_cast<VTKRenderer>(getActiveRenderer());
    if (ren != nullptr)
    {
        ren->updateCamera();
    }
}

void
QVTKViewer::setDebugAxesLength(double x, double y, double z)
{
    auto vtkRenderer = std::dynamic_pointer_cast<VTKRenderer>(getActiveRenderer());
    if (vtkRenderer != nullptr)
    {
        vtkRenderer->setAxesLength(x, y, z);
    }
}

void
QVTKViewer::setRenderingMode(const Renderer::Mode mode)
{
    if (!m_activeScene)
    {
        LOG(WARNING) << "Missing scene, can not set rendering mode.\n"
                     << "Use Viewer::setCurrentScene to setup scene.";
        return;
    }

    // Switch the renderer to the mode
    this->getActiveRenderer()->setMode(mode, false);

    updateModule();

    if (m_config->m_hideCurzor)
    {
        m_vtkRenderWindow->HideCursor();
    }

    if (m_config->m_hideBorder)
    {
        m_vtkRenderWindow->BordersOff();
    }

    if (m_config->m_fullScreen)
    {
        m_vtkRenderWindow->FullScreenOn();
    }
}

bool
QVTKViewer::initModule()
{
    if (!AbstractVTKViewer::initModule())
    {
        return false;
    }

    if (m_vtkRenderWindow->GetInteractor()->HasObserver(vtkCommand::StartEvent))
    {
        m_vtkRenderWindow->GetInteractor()->InvokeEvent(vtkCommand::StartEvent, nullptr);
    }
    m_vtkRenderWindow->Render();

    return true;
}

std::shared_ptr<KeyboardDeviceClient>
QVTKViewer::getKeyboardDevice() const
{
    return std::dynamic_pointer_cast<VTKInteractorStyle>(m_interactorStyle)->getKeyboardDeviceClient();
}

std::shared_ptr<MouseDeviceClient>
QVTKViewer::getMouseDevice() const
{
    return std::dynamic_pointer_cast<VTKInteractorStyle>(m_interactorStyle)->getMouseDeviceClient();
}

void
QVTKViewer::updateModule()
{
    // Process all Qt UI events (include render updates)
    if (m_executionType != ExecutionType::PARALLEL)
    {
        QApplication::processEvents(QEventLoop::ProcessEventsFlag::AllEvents, 7);
    }

    // Do render
    std::shared_ptr<VTKRenderer> ren = std::dynamic_pointer_cast<VTKRenderer>(getActiveRenderer());
    if (ren == nullptr)
    {
        return;
    }

    // Update Camera
    ren->updateCamera();

    // Call visual update on every scene object
    getActiveScene()->updateVisuals();
    // Update all the rendering delegates
    ren->updateRenderDelegates();

    // Automatically determine near and far planes (not used atm)
    //ren->getVtkRenderer()->ResetCameraClippingRange();

    // Render
    m_vtkRenderWindow->Render();
}