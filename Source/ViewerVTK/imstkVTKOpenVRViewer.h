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

#pragma once

#include "imstkAbstractVTKViewer.h"

namespace imstk
{
class OpenVRDeviceClient;

///
/// \class VTKViewer
///
/// \brief Subclasses viewer for the VTK rendering back-end
/// Creates vtk renderer for each scene. Forwards mouse and keyboard events
/// to the vtk renderwindow
///
class VTKOpenVRViewer : public AbstractVTKViewer
{
public:
    VTKOpenVRViewer(std::string name = "VTKOpenVRViewer");
    ~VTKOpenVRViewer() override      = default;

    ///
    /// \brief Destructor
    ///
    void setRenderingMode(const Renderer::Mode mode) override;

    ///
    /// \brief Set scene to be rendered
    ///
    void setActiveScene(std::shared_ptr<Scene> scene) override;

    ///
    /// \brief Transform to physical space
    ///
    void setPhysicalToWorldTransform(const Mat4d& physicalToWorldMatrix);

    ///
    /// \brief Get transform to physical space
    ///
    Mat4d getPhysicalToWorldTransform();

    ///
    /// \brief Get one of the device clients for VR
    ///
    std::shared_ptr<OpenVRDeviceClient> getVRDeviceClient(int deviceType);

    ///
    /// \brief Acquire the full list of VR devices tied to this viewer
    ///
    const std::vector<std::shared_ptr<OpenVRDeviceClient>>& getVRDeviceClients() const { return m_vrDeviceClients; }

    ///
    /// \brief VTKOpenVRViewer overrides to provide a non-rendering
    /// event processing loop (to deal with vsync blockage)
    ///
    void processEvents() override;

    ///
    /// \brief Toggle the controller model representation. Default off.
    ///
    void setControllerVisibility(const bool visible);

protected:
    bool initModule() override;

    void updateModule() override;

    std::vector<std::shared_ptr<OpenVRDeviceClient>> m_vrDeviceClients; ///> The VR controllers are tied to the view
#ifdef iMSTK_USE_OPENXR
    bool m_didFirstRender = false;                                  ///> Used for lazy initialization of models
#endif
};
} // imstk
