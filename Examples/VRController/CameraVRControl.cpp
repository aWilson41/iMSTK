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

#include "CameraVRControl.h"
#include "imstkCamera.h"
#include "imstkLogger.h"
#include "imstkVRDeviceClient.h"

using namespace imstk;

void
CameraVRControl::printControls()
{
    LOG(INFO) << "Mouse Scene Controls: Only usable in debug mode";
    LOG(INFO) << "----------------------------------------------------------------------";
    LOG(INFO) << " | Left Trackpad   - rotate view";
    LOG(INFO) << " | Right Trakcpad  - translate view";
    LOG(INFO) << "----------------------------------------------------------------------";
}

void
CameraVRControl::update(const double dt)
{
    // We may switch cameras on the controller
    m_deltaTransform = Mat4d::Identity();
    if (m_camera == nullptr)
    {
        return;
    }

    if (m_rotateDevice != nullptr)
    {
        m_rotateDevice->update();

        const Vec2d& pos  = m_rotateDevice->getTrackpadPosition();
        const Mat4d& view = m_camera->getView();
        m_camera->setView(
            view * mat4dRotation(Rotd(-pos[0] * m_rotateSpeedScale * dt, Vec3d(0.0, 1.0, 0.0))));
    }
    if (m_translateDevice != nullptr)
    {
        m_translateDevice->update();

        const Vec2d& pos = m_translateDevice->getTrackpadPosition();

        /*double dy = 0.0;
        if (m_translateDevice->getButton(2))
        {
            dy = m_translateVerticalSpeedScale;
        }
        if (m_translateDevice->getButton(3))
        {
            dy = -m_translateVerticalSpeedScale;
        }*/

        // The forward direction should be defined by the current view orientation around y

        const Mat4d& view = m_camera->getView();

        Vec3d t = Vec3d::Zero();
        Mat3d r = Mat3d::Identity();
        Vec3d s = Vec3d::Ones();
        mat4dTRS(view, t, r, s);
        // Decompose the x, y, z rotations
        r.col(0) = Vec3d(1.0, 0.0, 0.0);
        r.col(2) = Vec3d(0.0, 0.0, 1.0);

        const Vec3d forwardMovement = r * Vec3d(pos[0], 0.0, -pos[1]);
        //const Vec3d verticalMovement = Vec3d(0.0, dy, 0.0);

        m_deltaTransform =
            mat4dTranslate((forwardMovement /*+ verticalMovement*/) * m_translateSpeedScale * dt);
        m_camera->setView(view * m_deltaTransform);
    }
}