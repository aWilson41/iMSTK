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

#include "imstkNodeAnimation.h"

namespace imstk
{
///
/// \class VectorNodeAnimation
/// 
/// \brief Provides transforms of a node given time via samples its give
/// (usually imported from another software, such as Blender)
/// 
class KeyFrameNodeAnimation : public NodeAnimation
{
public:
    virtual ~KeyFrameNodeAnimation() override = default;

public:
    virtual void update(const double t) override
    {
        // todo: In some importers its possible not all 3 keys be present,
        // in which case the bind pose should be returned
        m_transform = (mat4dTranslate(getPosition(t)) *
            mat4dRotation(getRotation(t)) *
            mat4dScale(getScaling(t)));
    }

protected:
    Vec3d getPosition(double t) const
    {
        if (t <= m_positions.begin()->first)
        {
            return m_positions.begin()->second;
        }
        if (t >= m_positions.rbegin()->first)
        {
            return m_positions.rbegin()->second;
        }

        for (auto i = m_positions.begin(); i != std::prev(m_positions.end()); i++)
        {
            const double t1 = i->first;
            const double t2 = std::next(i)->first;
            if (t >= t1 && t < t2)
            {
                Vec3d  pos1 = i->second;
                Vec3d  pos2 = std::next(i)->second;
                double lerpT = (t - t1) / (t2 - t1);
                return pos1 + lerpT * (pos2 - pos1);
            }
        }
        return Vec3d(0.0, 0.0, 0.0);
    }

    Quatd getRotation(double t) const
    {
        if (t <= m_orientations.begin()->first)
        {
            return m_orientations.begin()->second;
        }
        if (t >= m_orientations.rbegin()->first)
        {
            return m_orientations.rbegin()->second;
        }

        for (auto i = m_orientations.begin(); i != std::prev(m_orientations.end()); i++)
        {
            const double t1 = i->first;
            const double t2 = std::next(i)->first;
            if (t >= t1 && t < t2)
            {
                Quatd  rot1 = i->second;
                Quatd  rot2 = std::next(i)->second;
                double lerpT = (t - t1) / (t2 - t1);
                return rot1.slerp(lerpT, rot2).normalized();
            }
        }
        return Quatd::Identity();
    }

    Vec3d getScaling(double t) const
    {
        if (t <= m_scalings.begin()->first)
        {
            return m_scalings.begin()->second;
        }
        if (t >= m_scalings.rbegin()->first)
        {
            return m_scalings.rbegin()->second;
        }

        for (auto i = m_scalings.begin(); i != std::prev(m_scalings.end()); i++)
        {
            const double t1 = i->first;
            const double t2 = std::next(i)->first;
            if (t >= t1 && t < t2)
            {
                Vec3d  s1 = i->second;
                Vec3d  s2 = std::next(i)->second;
                double lerpT = (t - t1) / (t2 - t1);
                return s1 + lerpT * (s2 - s1);
            }
        }
        return Vec3d(0.0, 0.0, 0.0);
    }

    bool hasScalings() const { return m_scalings.size() > 0; }
    bool hasOrientations() const { return m_orientations.size() > 0; }
    bool hasPositions() const { return m_positions.size() > 0; }

protected:
    friend class ObjectIO;

    std::map<double, Vec3d> m_positions;    ///> Positions keyed by time
    std::map<double, Quatd> m_orientations; ///> Orientations keyed by time
    std::map<double, Vec3d> m_scalings;     ///> Scalings keyed by time
};
}