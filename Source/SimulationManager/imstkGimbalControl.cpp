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

#include "imstkGimbalControl.h"
#include "imstkCamera.h"
#include "imstkLogger.h"
#include "imstkMouseDeviceClient.h"
#include "imstkScene.h"
#include "imstkSceneManager.h"
#include "imstkKeyboardDeviceClient.h"
#include "imstkGeometry.h"
#include "imstkOrientedBox.h"
#include "imstkLineMesh.h"
#include "imstkVisualModel.h"
#include "imstkSceneObject.h"
#include "imstkVecDataArray.h"
#include "imstkRenderMaterial.h"
#include "imstkViewer.h"

#include "imstkSphere.h"
#include "imstkCollisionUtils.h"

namespace imstk
{
static Vec3d globalAxes[3] =
{
    Vec3d(1.0, 0.0, 0.0),
    Vec3d(0.0, 1.0, 0.0),
    Vec3d(0.0, 0.0, 1.0)
};

static std::shared_ptr<LineMesh> createAxesMesh()
{
    auto verticesPtr = std::make_shared<VecDataArray<double, 3>>(4);
    (*verticesPtr)[0] = Vec3d::Zero();
    (*verticesPtr)[1] = Vec3d(0.5, 0.0, 0.0);
    (*verticesPtr)[2] = Vec3d(0.0, 0.5, 0.0);
    (*verticesPtr)[3] = Vec3d(0.0, 0.0, 0.5);

    auto indicesPtr = std::make_shared<VecDataArray<int, 2>>(3);
    (*indicesPtr)[0] = Vec2i(0, 1);
    (*indicesPtr)[1] = Vec2i(0, 2);
    (*indicesPtr)[2] = Vec2i(0, 3);

    auto colorsPtr = std::make_shared<VecDataArray<unsigned char, 3>>(3);
    (*colorsPtr)[0] = Color(1.0, 0.0, 0.0);
    (*colorsPtr)[1] = Color(0.0, 1.0, 0.0);
    (*colorsPtr)[2] = Color(0.0, 0.0, 1.0);

    auto axesLineMesh = std::make_shared<LineMesh>();
    axesLineMesh->initialize(verticesPtr, indicesPtr);
    axesLineMesh->setCellScalars("Colors", colorsPtr);
    return axesLineMesh;
}

GimbalControl::GimbalControl(std::shared_ptr<MouseDeviceClient> device) : MouseControl(device)
{
    // Scale widget geometry
    {
        m_gimbalScaleWidgetObject = std::make_shared<SceneObject>("Scale");

        std::shared_ptr<Geometry> geometries[4] = {
            std::make_shared<OrientedBox>(Vec3d(0.5, 0.0, 0.0), Vec3d(0.05, 0.05, 0.05)),
            std::make_shared<OrientedBox>(Vec3d(0.0, 0.5, 0.0), Vec3d(0.05, 0.05, 0.05)),
            std::make_shared<OrientedBox>(Vec3d(0.0, 0.0, 0.5), Vec3d(0.05, 0.05, 0.05)),
            createAxesMesh()
        };
        for (int i = 0; i < 4; i++)
        {
            auto visualModel = std::make_shared<VisualModel>();
            visualModel->setGeometry(geometries[i]);
            m_gimbalScaleWidgetObject->addVisualModel(visualModel);
        }
        m_gimbalScaleWidgetObject->getVisualModel(0)->getRenderMaterial()->setColor(Color::Red);
        m_gimbalScaleWidgetObject->getVisualModel(1)->getRenderMaterial()->setColor(Color::Green);
        m_gimbalScaleWidgetObject->getVisualModel(2)->getRenderMaterial()->setColor(Color::Blue);
        m_gimbalScaleWidgetObject->getVisualModel(3)->getRenderMaterial()->setLineWidth(6.0);
    }
    // Pan widget geometry
    {
        m_gimbalPanWidgetObject = std::make_shared<SceneObject>("Pan");
        m_gimbalPanWidgetObject->setVisualGeometry(createAxesMesh());
        m_gimbalPanWidgetObject->getVisualModel(0)->getRenderMaterial()->setLineWidth(6.0);
    }
    // Rotate widget geometry
    {
        m_gimbalRotateWidgetObject = std::make_shared<SceneObject>("Rotate");

        auto circleMeshY = std::make_shared<LineMesh>();

        const int numVerts = 35;
        auto verticesPtr = std::make_shared<VecDataArray<double, 3>>(numVerts);
        auto indicesPtr = std::make_shared<VecDataArray<int, 2>>(numVerts);
        for (int i = 0; i < numVerts; i++)
        {
            const double rad = static_cast<double>(i) / numVerts * PI * 2.0;
            Vec3d pos = Vec3d(cos(rad), 0.0, sin(rad)) * 0.5;
            (*verticesPtr)[i] = pos;
            (*indicesPtr)[i] = Vec2i(i, (i + 1) % numVerts);
        }
        circleMeshY->initialize(verticesPtr, indicesPtr);
        auto visualModelY = std::make_shared<VisualModel>();
        visualModelY->setGeometry(circleMeshY);
        visualModelY->getRenderMaterial()->setColor(Color::Green);
        m_gimbalRotateWidgetObject->addVisualModel(visualModelY);
        m_gimbalRotateWidgetObject->getVisualModel(0)->getRenderMaterial()->setLineWidth(6.0);

        auto circleMeshZ = std::make_shared<LineMesh>();
        circleMeshZ->initialize(
            std::make_shared<VecDataArray<double, 3>>(*verticesPtr),
            std::make_shared<VecDataArray<int, 2>>(*indicesPtr));
        circleMeshZ->rotate(Vec3d(1.0, 0.0, 0.0), PI_2, Geometry::TransformType::ApplyToData);
        auto visualModelZ = std::make_shared<VisualModel>();
        visualModelZ->setGeometry(circleMeshZ);
        visualModelZ->getRenderMaterial()->setColor(Color::Blue);
        m_gimbalRotateWidgetObject->addVisualModel(visualModelZ);
        m_gimbalRotateWidgetObject->getVisualModel(1)->getRenderMaterial()->setLineWidth(6.0);

        auto circleMeshX = std::make_shared<LineMesh>();
        circleMeshX->initialize(
            std::make_shared<VecDataArray<double, 3>>(*verticesPtr),
            std::make_shared<VecDataArray<int, 2>>(*indicesPtr));
        circleMeshX->rotate(Vec3d(0.0, 0.0, 1.0), PI_2, Geometry::TransformType::ApplyToData);
        auto visualModelX = std::make_shared<VisualModel>();
        visualModelX->setGeometry(circleMeshX);
        visualModelX->getRenderMaterial()->setColor(Color::Red);
        m_gimbalRotateWidgetObject->addVisualModel(visualModelX);
        m_gimbalRotateWidgetObject->getVisualModel(2)->getRenderMaterial()->setLineWidth(6.0);
    }
}

void
GimbalControl::printControls()
{
    LOG(INFO) << "Gimbal Controls:";
    LOG(INFO) << "----------------------------------------------------------------------";
    LOG(INFO) << " | Left click drag - rotate, pan, or scale";
    LOG(INFO) << "----------------------------------------------------------------------";
}

void
GimbalControl::OnButtonPress(const int button)
{
    if (button == LEFT_BUTTON && getEnabled())
    {
        // Now we can test intersections with the widget depending on the mode
        if (m_transformMode == TransformMode::ROTATE)
        {
            // Compute distance to each rotational circle to see which you clicked
        }
        else if (m_transformMode == TransformMode::PAN || m_transformMode == TransformMode::SCALE)
        {
            // Measure closest distance to each axes from ray and pick the closest one
            int minAxes = 0;
            double minDist = IMSTK_DOUBLE_MAX;
            Vec3d clickPtOnAxes;
            for (int i = 0; i < 3; i++)
            {
                const Vec3d pos = m_controlledGeometry->getTranslation();

                // Position along axes
                double dist;
                const Vec3d ptOnAxes = getMouseWorldPosOnAxes(pos, globalAxes[i], dist);
                printf("dist %d: %f\n", i, dist);
                if (dist < minDist)
                {
                    minDist = dist;
                    minAxes = i;
                    clickPtOnAxes = ptOnAxes;
                }
            }

            // At some angles two axes may be very close to the ray, its then desirable
            // to choose the one closest to the view

            // Store the current transform of the object as well as the location on the axes you clicked
            m_initMouseWorldPos = clickPtOnAxes;
            m_initTransform = m_controlledGeometry->getTransform();
            m_axesActivated = minAxes;
        }

        m_buttonDown = true;
    }
}

void
GimbalControl::OnButtonRelease(const int button)
{
    if (button == LEFT_BUTTON)
    {
        m_buttonDown = false;
    }
}

Vec3d
GimbalControl::getMouseWorldPosOnAxes(Vec3d pos, Vec3d axes, double& dist)
{
    std::shared_ptr<Camera> cam = m_viewer->getActiveScene()->getActiveCamera();

    // This mouse positions is bot left is 0,0. Ranges to 1.0, 1.0 to right and up
    const Vec2d mousePos = m_mouseDeviceClient->getPos();

    // We need ndc coordinates (-1, 1)
    const Vec3d rayDir = cam->getEyeRayDir(Vec2d(mousePos[0] * 2.0 - 1.0, mousePos[1] * 2.0 - 1.0));
    const Vec3d rayStart = cam->getPosition();

    double t, s;
    CollisionUtils::axesToAxesClosestPoint(
        pos, axes,
        rayStart, rayDir,
        t, s);
    Vec3d pA = pos + t * axes;
    Vec3d pB = rayStart + s * rayDir;
    dist = (pB - pA).norm();
    return pA;
}

void
GimbalControl::OnMouseMove(const Vec2d& pos)
{
    // Don't display widgets when not enabled & don't mose move
    if (!getEnabled())
    {
        m_gimbalPanWidgetObject->getVisualModel(0)->setIsVisible(false);
        return;
    }
    // Display the widget if enabled but not currently being clicked
    else if (!m_buttonDown)
    {
        m_gimbalPanWidgetObject->getVisualModel(0)->setIsVisible(true);
        return;
    }
    // Don't display the widget while being clicked
    else
    {
        m_gimbalPanWidgetObject->getVisualModel(0)->setIsVisible(false);
    }
    
    const Vec3d axes = globalAxes[m_axesActivated];

    Vec3d t, s;
    Mat3d r;
    mat4dTRS(m_initTransform, t, r, s);

    if (m_transformMode == TransformMode::ROTATE)
    {
    }
    else if (m_transformMode == TransformMode::PAN)
    {
        double dist;
        const Vec3d currMouseWorldPos =
            getMouseWorldPosOnAxes(m_controlledGeometry->getTranslation(), axes, dist);
        const Vec3d translation = (currMouseWorldPos - m_initMouseWorldPos).cwiseProduct(axes);

        t += translation;

        m_gimbalPanWidgetObject->getVisualModel(0)->getGeometry()->setTranslation(m_controlledGeometry->getTranslation());
        m_gimbalPanWidgetObject->getVisualModel(0)->getGeometry()->postModified();
    }
    else if (m_transformMode == TransformMode::SCALE)
    {
        double dist;
        const Vec3d currMouseWorldPos =
            getMouseWorldPosOnAxes(m_controlledGeometry->getTranslation(), axes, dist);
        const Vec3d deltaScale = (currMouseWorldPos - m_initMouseWorldPos).cwiseProduct(axes);

        s += deltaScale;
    }

    const Mat4d postTransform = mat4dTranslate(t) * mat4dRotation(r) * mat4dScale(s);
    m_controlledGeometry->setTransform(postTransform);
    m_controlledGeometry->postModified();
}

void
GimbalControl::update(const double dt)
{
    // Scale the axes with the view
    /*if (getEnabled())
    {
        if (m_transformMode == TransformMode::PAN)
        {
            const Mat4d invView = m_viewer->getActiveScene()->getActiveCamera()->getView().transpose();
            const double maxScale = std::max(
                std::max(invView.col(0).norm(), invView.col(1).norm()),
                invView.col(2).norm());

            m_gimbalPanWidgetObject->getVisualModel(0)->getGeometry()->setScaling(maxScale);
            m_gimbalPanWidgetObject->getVisualModel(0)->postModified();
        }
    }*/
}

void
GimbalControl::setEnabled(bool enable)
{
    m_enabled = enable;
}

bool
GimbalControl::getEnabled() const
{
    return m_enabled;
}
}