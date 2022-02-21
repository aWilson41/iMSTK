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

#include <vtkInteractorStyle3D.h>
#include <memory>

namespace imstk
{
class VRDeviceClient;
} // namespace imstk

///
/// \class vtkInteractorStyleVR
///
/// \brief VTK Interactor style for VR
///
class vtkInteractorStyleVR : public vtkInteractorStyle3D
{
public:
    static vtkInteractorStyleVR* New();
    vtkTypeMacro(vtkInteractorStyleVR, vtkInteractorStyle3D);

    void OnMove3D(vtkEventData* edata) override;

    // Alternatively could be used, switch to when dropping VTK 9.1 support
    // void SetupActions(vtkRenderWindowInteractor* iren) override;

    ///
    /// \brief Adds button actions
    ///
    void addButtonActions();

    ///
    /// \brief Adds thumbstick movement actions
    ///
    void addMovementActions();

#ifdef iMSTK_USE_OPENXR
    ///
    /// \brief Add haptics action
    ///
    void addHapticAction();
#endif

    std::shared_ptr<imstk::VRDeviceClient> getLeftControllerDeviceClient() const { return m_leftControllerDeviceClient; }
    std::shared_ptr<imstk::VRDeviceClient> getRightControllerDeviceClient() const { return m_rightControllerDeviceClient; }
    std::shared_ptr<imstk::VRDeviceClient> getHmdDeviceClient() const { return m_hmdDeviceClient; }

    vtkInteractorStyleVR();

protected:
    void OnButtonPress(vtkEventData* data, int buttonId);

public:
    std::shared_ptr<imstk::VRDeviceClient> m_leftControllerDeviceClient;
    std::shared_ptr<imstk::VRDeviceClient> m_rightControllerDeviceClient;
    std::shared_ptr<imstk::VRDeviceClient> m_hmdDeviceClient;
};
