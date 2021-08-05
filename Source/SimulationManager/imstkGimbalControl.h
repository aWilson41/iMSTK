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

#include "imstkMouseControl.h"

namespace imstk
{
class Geometry;
class MouseDeviceClient;
class SceneObject;
class Viewer;

///
/// \class GimbalControl
///
/// \brief
///
class GimbalControl : public MouseControl
{
public:
    enum class TransformMode
    {
        PAN,
        ROTATE,
        SCALE
    };
    ///
    /// \brief Defines space to do the transform in
    /// Global defines fixed basis for which to rotate, pan, scale with
    /// Local uses the local basis
    /// ex: (1,0,0) is always x in global, but in local, it changes after you rotate
    ///
    enum class SpaceMode
    {
        GLOBAL,
        LOCAL
    };

public:
    GimbalControl() = default;
    GimbalControl(std::shared_ptr<MouseDeviceClient> device);
    ~GimbalControl() override = default;

public:
    ///
    /// \brief Sets the viewer, the full viewer is needed for picking
    /// as window resolution is needed
    /// 
    void setViewer(std::shared_ptr<Viewer> viewer) { m_viewer = viewer; }

    ///
    /// \brief Enable the mouse control, independent of the debug mode
    ///
    void setEnabled(bool enable);

    ///
    /// \return true if the controls are enabled, either explicitly or debug is on in the scenecontrol
    ///
    bool getEnabled() const;

    void setTransformMode(const TransformMode& mode) { m_transformMode = mode; }
    const TransformMode& getTransformMode() const { return m_transformMode; }

    void setControlledGeometry(std::shared_ptr<Geometry> geom) { m_controlledGeometry = geom; }
    std::shared_ptr<Geometry> getControlledGeometry() const { return m_controlledGeometry; }

    std::shared_ptr<SceneObject> getGimbalPanObject() const { return m_gimbalPanWidgetObject; }
    std::shared_ptr<SceneObject> getGimbalScaleObject() const { return m_gimbalScaleWidgetObject; }
    std::shared_ptr<SceneObject> getGimbalRotateObject() const { return m_gimbalRotateWidgetObject; }

public:
    void printControls() override;

    ///
    /// \brief On the mouse scene control button press
    ///
    void OnButtonPress(const int key) override;
    void OnButtonRelease(const int key) override;
    void OnMouseMove(const Vec2d& pos) override;

    void update(const double dt) override;

protected:
    ///
    /// \brief Returns the mouse world position by nearest point projection onto the axes
    /// 
    Vec3d getMouseWorldPosOnAxes(Vec3d pos, Vec3d axes, double& dist);

protected:
    std::shared_ptr<Geometry>    m_controlledGeometry = nullptr;
    std::shared_ptr<SceneObject> m_gimbalPanWidgetObject = nullptr;
    std::shared_ptr<SceneObject> m_gimbalScaleWidgetObject = nullptr;
    std::shared_ptr<SceneObject> m_gimbalRotateWidgetObject = nullptr;
    std::shared_ptr<Viewer> m_viewer = nullptr;

    TransformMode m_transformMode = TransformMode::SCALE;
    SpaceMode     m_spaceMode = SpaceMode::GLOBAL;

    Vec3d m_initMouseWorldPos = Vec3d::Zero(); // Mouse world position at the time of clicking
    Mat4d m_initTransform = Mat4d::Identity(); ///> Transform of the object before

    bool m_enabled = true; ///> Whether the widget is enabled
    bool m_buttonDown = false;

    int m_axesActivated = -1;
};
}