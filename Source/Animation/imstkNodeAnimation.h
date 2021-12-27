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

#include "imstkMath.h"

namespace imstk
{
class AnimationNode;

///
/// \class NodeAnimation
/// 
/// \brief A node animation is something that provides a transform for a node
/// 
class NodeAnimation
{
protected:
    NodeAnimation() = default;
public:
    virtual ~NodeAnimation() = default;

public:
    virtual void update(const double t) = 0;

    const Mat4d& getTransform() const { return m_transform; }

    std::shared_ptr<AnimationNode> getNode() const { return m_animationNode; }
    virtual void setNode(std::shared_ptr<AnimationNode> node) { m_animationNode = node; }

protected:
    Mat4d m_transform;
    std::shared_ptr<AnimationNode> m_animationNode;
};
}
