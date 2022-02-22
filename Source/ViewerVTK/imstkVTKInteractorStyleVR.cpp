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

#include "imstkVTKInteractorStyleVR.h"
#include "imstkLogger.h"
#include "imstkVRDeviceClient.h"

#include <vtkEventData.h>
#include <vtkMath.h>
#include <vtkObjectFactory.h>
#ifdef iMSTK_USE_OPENXR
#include <vtkOpenXRRenderWindowInteractor.h>
using vtkImstkVRRenderWindowInteractor = vtkOpenXRRenderWindowInteractor;
#else
#include <vtkOpenVRRenderWindowInteractor.h>
using vtkImstkVRRenderWindowInteractor = vtkOpenVRRenderWindowInteractor;
#endif

vtkStandardNewMacro(vtkInteractorStyleVR);

vtkInteractorStyleVR::vtkInteractorStyleVR()
{
    // Setup the VR device clients
    m_leftControllerDeviceClient  = imstk::VRDeviceClient::New(LEFT_CONTROLLER);
    m_rightControllerDeviceClient = imstk::VRDeviceClient::New(RIGHT_CONTROLLER);
    m_hmdDeviceClient = imstk::VRDeviceClient::New(VR_HMD);
}

void
vtkInteractorStyleVR::OnButtonPress(vtkEventData* data, int buttonId)
{
    vtkEventDataForDevice*   eventDataButton = data->GetAsEventDataForDevice();
    const vtkEventDataAction action = eventDataButton->GetAction();
    const vtkEventDataDevice device = eventDataButton->GetDevice();

    if (device == vtkEventDataDevice::LeftController)
    {
        if (action == vtkEventDataAction::Press)
        {
            m_leftControllerDeviceClient->emitButtonPress(buttonId);
        }
        else if (action == vtkEventDataAction::Release)
        {
            m_leftControllerDeviceClient->emitButtonRelease(buttonId);
        }
    }
    else if (device == vtkEventDataDevice::RightController)
    {
        if (action == vtkEventDataAction::Press)
        {
            m_rightControllerDeviceClient->emitButtonPress(buttonId);
        }
        else if (action == vtkEventDataAction::Release)
        {
            m_rightControllerDeviceClient->emitButtonRelease(buttonId);
        }
    }
}

void
vtkInteractorStyleVR::addMovementActions()
{
    auto iren = vtkImstkVRRenderWindowInteractor::SafeDownCast(GetInteractor());
    if (iren == nullptr)
    {
        LOG(FATAL) << "Cannot add movement actions, wrong interactor";
        return;
    }

#ifdef iMSTK_USE_OPENXR
    iren->AddAction("leftgripmovement",
#else
    CHECK(iren->GetInitialized()) << "Cannot addMovementActions to style until "
        "interactor has been initialized";
    iren->AddAction("/actions/vtk/in/LeftGripMovement", true,
#endif
        [this](vtkEventData * ed)
    {
        vtkEventDataDevice3D* edd = ed->GetAsEventDataDevice3D();
        const double* pos = edd->GetTrackPadPosition();
        m_leftControllerDeviceClient->setTrackpadPosition(imstk::Vec2d(pos[0], pos[1]));
        });
#ifdef iMSTK_USE_OPENXR
    iren->AddAction("rightgripmovement",
#else
    iren->AddAction("/actions/vtk/in/RightGripMovement", true,
#endif
        [this](vtkEventData * ed)
    {
        vtkEventDataDevice3D* edd = ed->GetAsEventDataDevice3D();
        const double* pos = edd->GetTrackPadPosition();
        m_rightControllerDeviceClient->setTrackpadPosition(imstk::Vec2d(pos[0], pos[1]));
        });
}

void
vtkInteractorStyleVR::addButtonActions()
{
    auto iren = vtkImstkVRRenderWindowInteractor::SafeDownCast(GetInteractor());
    if (iren == nullptr)
    {
        LOG(FATAL) << "Cannot add button actions, wrong interactor";
        return;
    }

    // Called when buttons are pressed/released
    std::array<std::string, 6> buttonActionNames =
    {
#ifdef iMSTK_USE_OPENXR
        "button0pressed",
        "button1pressed",
        "button2pressed",
        "button3pressed",
        "grippressed",
        "triggerpressed"
#else
        "/actions/vtk/in/Button0Pressed",
        "/actions/vtk/in/Button1Pressed",
        "/actions/vtk/in/Button2Pressed",
        "/actions/vtk/in/Button3Pressed",
        "/actions/vtk/in/GripPressed",
        "/actions/vtk/in/TriggerPressed"
#endif
    };
    for (int i = 0; i < 6; i++)
    {
#ifdef iMSTK_USE_OPENXR
        iren->AddAction(buttonActionNames[i],
#else
        CHECK(iren->GetInitialized()) << "Cannot addButtonActions to style until "
            "interactor has been initialized";
        iren->AddAction(buttonActionNames[i], false,
#endif
            [this, i](vtkEventData * ed)
        {
            OnButtonPress(ed, i);
            });
    }
}

void
vtkInteractorStyleVR::addHapticAction()
{
    auto iren = vtkImstkVRRenderWindowInteractor::SafeDownCast(GetInteractor());
    if (iren == nullptr)
    {
        LOG(FATAL) << "Cannot add haptic actions, wrong interactor";
        return;
    }

#ifdef iMSTK_USE_OPENXR
    // Allow the devices to call vibration without it having to hold the interactor
    m_leftControllerDeviceClient->setVibrationFunc(
        [iren](double amplitude, double duration, double frequency)
        {
            iren->ApplyVibration("left_haptic", vtkOpenXRManager::ControllerIndex::Left, amplitude, duration, frequency);
        });
    m_rightControllerDeviceClient->setVibrationFunc(
        [iren](double amplitude, double duration, double frequency)
        {
            iren->ApplyVibration("right_haptic", vtkOpenXRManager::ControllerIndex::Right, amplitude, duration, frequency);
        });
#endif
}

void
vtkInteractorStyleVR::OnMove3D(vtkEventData* eventData)
{
    if (eventData->GetType() != vtkCommand::Move3DEvent)
    {
        return;
    }
    auto eventDataDevice = static_cast<vtkEventDataDevice3D*>(eventData);

                                                              if (vtkEventDataDevice::LeftController == eventDataDevice->GetDevice())
    {
        imstk::Vec3d pos;
        eventDataDevice->GetWorldPosition(pos.data());
        double orientation[4];
        eventDataDevice->GetWorldOrientation(orientation);
        m_leftControllerDeviceClient->setPose(pos, imstk::Quatd(imstk::Rotd(vtkMath::RadiansFromDegrees(orientation[0]),
            imstk::Vec3d(orientation[1], orientation[2], orientation[3]))));
    }
                                                              else if (vtkEventDataDevice::RightController == eventDataDevice->GetDevice())
    {
        imstk::Vec3d pos;
        eventDataDevice->GetWorldPosition(pos.data());
        double orientation[4];
        eventDataDevice->GetWorldOrientation(orientation);
        m_rightControllerDeviceClient->setPose(pos, imstk::Quatd(imstk::Rotd(vtkMath::RadiansFromDegrees(orientation[0]),
            imstk::Vec3d(orientation[1], orientation[2], orientation[3]))));
    }
                                                              else if (vtkEventDataDevice::HeadMountedDisplay == eventDataDevice->GetDevice())
    {
        imstk::Vec3d pos;
        eventDataDevice->GetWorldPosition(pos.data());
        double orientation[4];
        eventDataDevice->GetWorldOrientation(orientation);
        m_hmdDeviceClient->setPose(pos, imstk::Quatd(imstk::Rotd(vtkMath::RadiansFromDegrees(orientation[0]),
            imstk::Vec3d(orientation[1], orientation[2], orientation[3]))));
    }
}