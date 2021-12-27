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
class Attachment;

///
/// \brief
///
class AnimationNode
{
public:
    std::string name = "";

    Mat4d localBindTransform = Mat4d::Identity(); ///> Transform of node in relation to parent
    Mat4d boneOffset = Mat4d::Identity();         // Assumes one bone per node

    AnimationNode* parent = nullptr;
    std::vector<std::shared_ptr<AnimationNode>> children;

    std::vector<std::shared_ptr<Attachment>> attachments;
};
}