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
class PointSet;

struct VertexWeight
{
    int vertexId;
    double weight;
};

///
/// \brief An attachment is something that recieves/applies the pose of a bone
/// 
class Attachment
{
protected:
    Attachment() = default;
public:
    virtual ~Attachment() = default;

public:
    ///
    /// \brief Update the state of the attachment using the global pose
    ///
    virtual void update(const Mat4d& pose) = 0;

public:
    AnimationNode* parent = nullptr; ///> The node this attachment is on
};

///
/// \Brief This attachment applies weighted transforms on a set
/// of vertices on PointSet geometry given by weights
///
class MeshWeightAttachment : public Attachment
{
public:
    ~MeshWeightAttachment() override = default;

public:
    ///
    /// \brief Update the state of the attachment using the global pose
    /// Updates the vertices by bone weights
    ///
    void update(const Mat4d& pose) override;

public:
    std::vector<VertexWeight> weights;
    std::shared_ptr<PointSet> geometry;
};

///
/// \brief This attachment applies a transform to a PointSet
///
class MeshAttachment : public Attachment
{
public:
    ~MeshAttachment() override = default;

public:
    ///
    /// \brief Update the state of the attachment using the global pose
    /// Updates the transform of the mesh
    ///
    void update(const Mat4d& pose) override;

public:
    std::shared_ptr<PointSet> geometry;
};
}