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

#include "imstkVRDeviceClient.h"
#include "imstkLogger.h"

namespace imstk
{
std::shared_ptr<VRDeviceClient>
VRDeviceClient::New(DeviceType deviceType)
{
    return std::shared_ptr<VRDeviceClient>(new VRDeviceClient(deviceType));
}

void
VRDeviceClient::emitButtonTouched(const int buttonId)
{
    const int prevButtonState = m_buttons[buttonId];
    m_buttons[buttonId] = BUTTON_TOUCHED;
    if (prevButtonState != BUTTON_TOUCHED)
    {
        this->postEvent(ButtonEvent(VRDeviceClient::buttonStateChanged(), buttonId, BUTTON_TOUCHED));
    }
}

void
VRDeviceClient::emitButtonUntouched(const int buttonId)
{
    const int prevButtonState = m_buttons[buttonId];
    m_buttons[buttonId] = BUTTON_UNTOUCHED;
    if (prevButtonState != BUTTON_UNTOUCHED)
    {
        this->postEvent(ButtonEvent(VRDeviceClient::buttonStateChanged(), buttonId, BUTTON_UNTOUCHED));
    }
}

void
VRDeviceClient::emitButtonPress(const int buttonId)
{
    const int prevButtonState = m_buttons[buttonId];
    m_buttons[buttonId] = BUTTON_PRESSED;
    if (prevButtonState != BUTTON_PRESSED)
    {
        this->postEvent(ButtonEvent(VRDeviceClient::buttonStateChanged(), buttonId, BUTTON_PRESSED));
    }
}

void
VRDeviceClient::emitButtonRelease(const int buttonId)
{
    const int prevButtonState = m_buttons[buttonId];
    m_buttons[buttonId] = BUTTON_RELEASED;
    if (prevButtonState != BUTTON_RELEASED)
    {
        this->postEvent(ButtonEvent(VRDeviceClient::buttonStateChanged(), buttonId, BUTTON_RELEASED));
    }
}

void
VRDeviceClient::applyVibration(const double amplitude, const double duration, const double frequency)
{
    // No way to ask the interactor to apply this here
    if (m_vibrationFunc != nullptr)
    {
        m_vibrationFunc(amplitude, duration, frequency);
    }
    else
    {
        LOG(FATAL) << "Tried to apply vibration on device without vibration function";
    }
}
} // namespace imstk