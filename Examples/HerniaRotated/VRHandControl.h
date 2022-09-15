/*
** This file is part of the Interactive Medical Simulation Toolkit (iMSTK)
** iMSTK is distributed under the Apache License, Version 2.0.
** See accompanying NOTICE for details.
*/

#include "imstkAnalyticalGeometry.h"
#include "imstkDeviceControl.h"
#include "imstkOpenVRDeviceClient.h"
#include "imstkPbdObject.h"
#include "imstkPbdObjectGrasping.h"
#include "imstkVisualModel.h"

namespace imstk
{
///
/// \brief Defines control scheme for VR controller hands and the associated
/// tool
/// 
class VRHandControl : public DeviceControl
{
public:
    VRHandControl(const std::string& name = "VRHandControl") : DeviceControl(name) { }
    ~VRHandControl() override = default;

    void setDevice(std::shared_ptr<DeviceClient> device) override
    {
        if (device == nullptr)
        {
            return;
        }

        // Remove old observer if it exists
        if (m_deviceClient != nullptr)
        {
            disconnect(m_deviceClient, this, &DeviceClient::buttonStateChanged);
        }

        // Set the new device
        m_deviceClient = std::dynamic_pointer_cast<OpenVRDeviceClient>(device);
        DeviceControl::setDevice(device);

        // Subscribe to the device clients events
        connect(device, &DeviceClient::buttonStateChanged, this, &VRHandControl::buttonStateChanged);
    }

    virtual void buttonStateChanged(ButtonEvent* e)
    {
        // Tool Grasping toggles on and off
        if (e->m_button == m_toolGraspButton)
        {
            if (e->m_buttonState == BUTTON_PRESSED)
            {
                bool hasConstraints = false;
                for (auto grasper : m_toolGraspers)
                {
                    if (grasper->hasConstraints())
                    {
                        hasConstraints = true;
                    }
                }

                if (hasConstraints)
                {
                    for (auto grasper : m_toolGraspers)
                    {
                        grasper->endGrasp();
                    }
                }
                else
                {
                    for (auto grasper : m_toolGraspers)
                    {
                        auto graspGeom = std::dynamic_pointer_cast<AnalyticalGeometry>(
                            m_handObj->getCollidingGeometry());
                        grasper->beginCellGrasp(graspGeom);
                    }
                }
            }
        }
        // Stomach grasping must be held to grasp, doesn't lock
        if (e->m_button == m_stomachGraspButton)
        {
            if (e->m_buttonState == BUTTON_PRESSED)
            {
                int toolHeldIndex = -1;
                for (int i = 0; i < m_toolGraspers.size(); i++)
                {
                    if (m_toolGraspers[i]->hasConstraints())
                    {
                        toolHeldIndex = i;
                    }
                }

                if (toolHeldIndex != -1)
                {
                    std::shared_ptr<PbdObject> toolObj =
                        m_toolGraspers[toolHeldIndex]->getObjectToGrasp();
                    // Use a slighty larger capsule at the tip
                    auto graspCapsule = std::dynamic_pointer_cast<AnalyticalGeometry>(
                        toolObj->getVisualModel(1)->getGeometry());

                    for (auto grasper : m_organGraspers)
                    {
                        grasper->beginCellGrasp(graspCapsule);
                    }
                }
            }
            else if (e->m_buttonState == BUTTON_RELEASED)
            {
                for (auto grasper : m_organGraspers)
                {
                    grasper->endGrasp();
                }
            }
        }
    }

    void addToolGrasper(std::shared_ptr<PbdObjectGrasping> grasper)
    {
        m_toolGraspers.push_back(grasper);
    }

    void addOrganGrasper(std::shared_ptr<PbdObjectGrasping> grasper)
    {
        m_organGraspers.push_back(grasper);
    }

    void setGraspButtons(const int toolGraspButton, const int stomachGraspButton)
    {
        m_toolGraspButton = toolGraspButton;
        m_stomachGraspButton = stomachGraspButton;
    }

    void setHandObject(std::shared_ptr<PbdObject> handObj) { m_handObj = handObj; }

protected:
    std::shared_ptr<OpenVRDeviceClient> m_deviceClient = nullptr;
    std::vector<std::shared_ptr<PbdObjectGrasping>> m_toolGraspers;
    std::vector<std::shared_ptr<PbdObjectGrasping>> m_organGraspers;
    std::shared_ptr<PbdObject> m_handObj = nullptr; // The hand to grasp the tool
    // The tool the hand is grasping is the one that
    std::shared_ptr<PbdObject> m_currentToolObj = nullptr;
    int m_toolGraspButton = 6;
    int m_stomachGraspButton = 4;
};
}