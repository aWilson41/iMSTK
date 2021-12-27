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

#include "imstkAbstractAnimationModel.h"

#include <unordered_map>
#include <unordered_set>

namespace imstk
{
class AnimationNode;
class NodeAnimation;
class Geometry;
class LineMesh;

///
/// \brief Animation of a set of nodes
///
struct AnimationTrack
{
public:
    //double duration = 0.0;
    double timeScale = 1.0;
    double startTime = 0.0;
    double t = 0.0; // Local space t

    std::vector<std::shared_ptr<NodeAnimation>> channels;                                                               //> Each node would use a channel
    std::unordered_map<std::shared_ptr<AnimationNode>, std::vector<std::shared_ptr<NodeAnimation>>> m_nodeToAnimations; //> Each node can have a set of nodal animations on it
};

///
/// \class SkinnedAnimationModel
///
/// \brief node based animation of a pointset
///
class AnimationModel : public AbstractAnimationModel
{
public:
    AnimationModel();
    virtual ~AnimationModel() override = default;

public:
    ///
    /// \brief Add an animation to the set of possible animations (the stash) for this model
    /// 
    void addAnimation(std::string name, std::shared_ptr<AnimationTrack> anim);

    ///
    /// \brief Get/Set root node
    ///
    void setRootNode(std::shared_ptr<AnimationNode> rootNode);
    std::shared_ptr<AnimationNode> getRootNode() const { return m_rootNode; }

    ///
    /// \brief Get/Set the transform tacked onto the root
    ///
    void setTransform(const Mat4d& transform) { m_transform = transform; }
    const Mat4d& getTransform() const { return m_transform; }

    ///
    /// \brief Get bone by name
    ///
    std::shared_ptr<AnimationNode> getNode(std::string name) const
    {
        return (m_namedBoneMap.count(name) == 0) ? nullptr : m_namedBoneMap.at(name);
    }

    ///
    /// \brief Set the time in the animation
    ///
    virtual void setTime(const double t) override;

    ///
    /// \brief Add an animation to the active set
    /// \param unique name/key of the animation
    /// \param starting time of animation, if < current time it will play
    ///  immediately (default)
    ///
    void playAnimation(std::string animationName, const double startTime = -1.0, const double timeScale = 1.0);

    ///
    /// \brief Remove an animation from the active set
    /// 
    void stopAnimation(std::string animationName, const bool stopWhenComplete = false);

public:
    ///
    /// \brief Get the line mesh of debug bones if debug skeleton is on
    /// Leaf nodes are not present (info not available assimp)
    ///
    std::shared_ptr<LineMesh> getDebugSkeleton() { return m_genDebugSkeleton ? m_debugSkeletonGeometry : nullptr; }

    ///
    /// \brief If on, will generate debug skeleton everytime the bones are updated
    /// default false
    ///
    void setGenDebugSkeleton(const bool genDebugSkeleton) { m_genDebugSkeleton = genDebugSkeleton; }

    ///
    /// \brief If debug skeleton on, will update the skeleton to the current bone configuration
    ///
    void updateDebugSkeleton();

protected:
    // The current animations
    std::unordered_map<std::string, std::shared_ptr<AnimationTrack>> m_animationTracks;

    // The available animations for this model
    std::unordered_map<std::string, std::shared_ptr<AnimationTrack>> m_animationStash;

    // Nodes
    std::shared_ptr<AnimationNode> m_rootNode;
    std::vector<std::shared_ptr<AnimationNode>> m_nodes; ///> All the nodes (root and all children)
    Mat4d m_transform = Mat4d::Identity(); ///> Top level transform applied to root

    // Debug Skeleton
    std::shared_ptr<LineMesh> m_debugSkeletonGeometry;
    bool m_genDebugSkeleton = false;

    // Bones
    std::unordered_map<AnimationNode*, Mat4d> m_globalTransforms; ///> Global transforms of all the bones
    std::unordered_map<std::string, std::shared_ptr<AnimationNode>> m_namedBoneMap;

    // Misc
    std::unordered_set<std::shared_ptr<Geometry>> m_resetGeometries;    ///> Geometries effected
    std::unordered_set<std::shared_ptr<Geometry>> m_modifiedGeometries; ///> Geometries effected
};
} // imstk
