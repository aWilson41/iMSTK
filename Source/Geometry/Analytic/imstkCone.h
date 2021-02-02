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

#include "imstkAnalyticalGeometry.h"

namespace imstk
{
///
/// \class Cone
///
/// \brief Cone geometry
///
class Cone : public AnalyticalGeometry
{
public:
    explicit Cone(const Vec3d& pos = Vec3d(0.0, 0.0, 0.0), const double height = 1.0, const double angle = PI / 4.0, const Vec3d& orientationAxis = Vec3d(0.0, 1.0, 0.0),
        const std::string& name = std::string("defaultCone")) : AnalyticalGeometry(Type::Cone, name)
    {
        setPosition(pos);
        setOrientationAxis(orientationAxis);
        setHeight(height);
        setAngle(angle);
    }
    ~Cone() override = default;

public:
    ///
    /// \brief Print the cube info
    ///
    void print() const override;

    ///
    /// \brief Returns the volume of the cone
    ///
    double getVolume() override { return PI * std::pow(getRadius(), 2.0) * m_height / 3.0; }

    ///
    /// \brief Returns the height of the cone
    ///
    double getHeight(DataType type = DataType::PostTransform);

    ///
    /// \brief Returns the angle (radians) of the nose of the cone
    /// 
    double getAngle() const { return m_angle; }

    ///
    /// \brief Returns the radius of bottom of the cone
    /// 
    double getRadius() const { return std::tan(m_angle * 0.5) * m_height; }

    std::string getTypeName() const override { return "Cone"; }

    ///
    /// \brief Sets the height of the cone
    ///
    void setHeight(const double h)
    {
        if (m_height == h)
        {
            return;
        }
        m_height = h;
        m_dataModified = true;
        m_transformApplied = false;
    }

    ///
    /// \brief Sets the angle of the cone (radians)
    /// 
    void setAngle(const double angle)
    {
        if (m_angle == angle)
        {
            return;
        }
        m_angle = angle;
        m_dataModified = true;
        m_transformApplied = false;
    }

    ///
    /// \brief Returns signed distance to surface at pos
    /// \todo Doesn't support orientation yet
    ///
    double getFunctionValue(const Vec3d& pos) const override
    {
        // Rotate sample point into aligned frame (with center of cone base/m_position at origin)
        const Mat3d invRotate = Quatd::FromTwoVectors(m_orientationAxisPostTransform, UP_VECTOR).matrix();
        //const Vec3d orientedPos = invRotate * (pos - m_position);

        const Vec3d centeredPos = invRotate * (pos - m_position) - Vec3d(0.0, m_height, 0.0);

        const Vec2d q = m_height * Vec2d(tan(m_angle), -1.0);

        const Vec2d w = Vec2d(Vec2d(centeredPos[0], centeredPos[2]).norm(), centeredPos[1]);
        const Vec2d a = w - q * std::min(std::max(w.dot(q) / q.dot(q), 0.0), 1.0);
        const Vec2d b = w - q.cwiseProduct(Vec2d(std::min(std::max(w[0] / q[0], 0.0), 1.0), 1.0));
        const double k = ((0.0 < q[1]) - (q[1] < 0.0));
        const double d = std::min(a.dot(a), b.dot(b));
        const double s = std::max(k * (w[0] * q[1] - w[1] * q[0]), k * (w[1] - q[1]));
        return sqrt(d) * ((0.0 < s) - (s < 0.0));
    }

protected:
    friend class VTKConeRenderDelegate;

    void applyScaling(const double s) override;
    void updatePostTransformData() const override;

    double m_height = 1.0;              ///> Height of the cone
    mutable double m_heightPostTransform = 1.0; ///> Height of cone after scaling
    double m_angle = 1.57 * 0.5;        ///> Angle of the cone
};
}
