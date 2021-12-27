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

#include "imstkGeometry.h"

namespace imstk
{
///
/// \class AbstractAnimationModel
///
/// \brief Abstract class for animations in imstk, defines basic animation control
/// interface
///
class AbstractAnimationModel
{
public:
    enum class InterpolationType
    {
        NEAREST,
        LINEAR
    };

protected:
    AbstractAnimationModel() = default;

public:
    virtual ~AbstractAnimationModel() = default;

public:
    InterpolationType getInterpolationType() const { return m_interpolationType; }
    void setInterpolationType(InterpolationType type) { m_interpolationType = type; }

    double getStartTime() const { return m_startTime; }
    void setStartTime(const double startTime) { m_startTime = startTime; }

    double getDt() const { return m_dt; }
    void setDt(const double dt) { m_dt = dt; }

public:
    ///
    /// \brief Update animation in time according to internal advance properties
    ///
    virtual void update()
    {
        m_t += m_dt;
        setTime(m_t);
    }

    ///
    /// \brief Reset animation to start time
    ///
    virtual void reset() { setTime(m_startTime); }

    ///
    /// \brief Set the time of the animation
    ///
    virtual void setTime(const double t) = 0;

protected:
    std::shared_ptr<Geometry> m_geometry  = nullptr;
    InterpolationType m_interpolationType = InterpolationType::LINEAR;
    double m_startTime = 0.0;
    double m_dt = 2.5; ///> ms to advance every update
    double m_t  = 0.0;  ///> Current time of the animation
};
} // imstk
