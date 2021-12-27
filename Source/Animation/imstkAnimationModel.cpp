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

#include "imstkAnimationModel.h"
#include "imstkAnimationNode.h"
#include "imstkAttachment.h"
#include "imstkLineMesh.h"
#include "imstkLogger.h"
#include "imstkNodeAnimation.h"
#include "imstkVecDataArray.h"

#include <stack>

namespace imstk
{
static void
updateBones(std::shared_ptr<AnimationNode> node, const Mat4d& parentTransform,
            std::unordered_map<AnimationNode*, Mat4d>& globalTransforms,
            std::unordered_map<std::string, std::shared_ptr<AnimationTrack>>& anim)
{
    Mat4d globalTransform;

    // If animated use the channels
    bool isAnimated = false;
    Mat4d localTransform = Mat4d::Identity();
    for (auto i : anim)
    {
        std::shared_ptr<AnimationTrack> track = i.second;
        for (auto nodeAnim : track->channels)
        {
            // If the animation effects this node (switch to map later)
            if (nodeAnim->getNode() == node)
            {
                // Update and get the transform of the animation
                nodeAnim->update(track->t);
                // No blending yet, just overwrites
                localTransform = nodeAnim->getTransform();
                isAnimated = true;
            }
        }
    }
    // If not animated use the local bind transform
    if (!isAnimated)
    {
        localTransform = node->localBindTransform;
    }

    globalTransform = parentTransform * localTransform;
    globalTransforms[node.get()] = globalTransform;

    for (int i = 0; i < node->children.size(); i++)
    {
        std::shared_ptr<AnimationNode> childNode = node->children[i];
        updateBones(childNode, globalTransform, globalTransforms, anim);
    }
}

AnimationModel::AnimationModel() : AbstractAnimationModel(),
    m_debugSkeletonGeometry(std::make_shared<LineMesh>())
{
    srand(static_cast<unsigned int>(time(NULL)));
}

void
AnimationModel::addAnimation(std::string name, std::shared_ptr<AnimationTrack> anim)
{
    m_animationStash[name] = anim;
}

void
AnimationModel::setRootNode(std::shared_ptr<AnimationNode> rootNode)
{
    this->m_rootNode = rootNode;

    // Collect all the nodes and geometries
    m_nodes.clear();
    std::shared_ptr<AnimationNode>             currNode = rootNode;
    std::stack<std::shared_ptr<AnimationNode>> nodeStack;
    nodeStack.push(rootNode);
    while (!nodeStack.empty())
    {
        std::shared_ptr<AnimationNode> node = nodeStack.top();
        m_namedBoneMap[node->name] = node;
        nodeStack.pop();

        // Gather weighted geometries to zero out
        for (auto attachment : node->attachments)
        {
            if (auto weightedMeshAttachment = std::dynamic_pointer_cast<MeshWeightAttachment>(attachment))
            {
                m_resetGeometries.insert(weightedMeshAttachment->geometry);
                m_modifiedGeometries.insert(weightedMeshAttachment->geometry);
                //m_geometryRoots[weightedMeshAttachment->geometry] = node.get();
            }
            else if (auto meshAttachment = std::dynamic_pointer_cast<MeshAttachment>(attachment))
            {
                m_modifiedGeometries.insert(meshAttachment->geometry);
                //m_geometryRoots[meshAttachment->geometry] = node.get();
            }
        }

        m_nodes.push_back(node);
        for (auto child : node->children)
        {
            nodeStack.push(child);
        }
    }
}

void
AnimationModel::setTime(const double t)
{
    // Update times of each track
    for (auto i : m_animationTracks)
    {
        i.second->t = (t - i.second->startTime) * i.second->timeScale;
    }

    // Update the global transforms of every bone
    updateBones(m_rootNode, m_transform, m_globalTransforms, m_animationTracks);

    // Update debug skeleton (if debug skeleton on)
    updateDebugSkeleton();

    // zero out the vertices, WeightedMeshAttachments will need to sum into them
    for (auto geom : m_resetGeometries)
    {
        std::shared_ptr<PointSet> pointSet = std::dynamic_pointer_cast<PointSet>(geom);
        if (pointSet != nullptr)
        {
            std::shared_ptr<VecDataArray<double, 3>> verticesPtr = pointSet->getVertexPositions(Geometry::DataType::PostTransform);
            verticesPtr->fill(Vec3d(0.0, 0.0, 0.0));
            std::shared_ptr<VecDataArray<double, 3>> vertexNormalsPtr = pointSet->getVertexNormals();
            vertexNormalsPtr->fill(Vec3d(0.0, 0.0, 0.0));
        }
    }

    // Update the attachments of the bones
    for (auto node : m_nodes)
    {
        const Mat4d& transform = m_globalTransforms[node.get()];

        for (auto attachment : node->attachments)
        {
            attachment->update(transform);
        }
    }

    for (auto geom : m_modifiedGeometries)
    {
        if (auto pointSet = std::dynamic_pointer_cast<PointSet>(geom))
        {
            std::shared_ptr<VecDataArray<double, 3>> vertexNormalsPtr = pointSet->getVertexNormals();
            VecDataArray<double, 3>& vertexNormals    = *vertexNormalsPtr;
            for (int i = 0; i < vertexNormals.size(); i++)
            {
                vertexNormals[i].normalize();
            }
            vertexNormalsPtr->postModified();
        }
        geom->postModified();
    }
}

void
AnimationModel::playAnimation(std::string animationName, const double startTime, const double timeScale)
{
    if (m_animationStash.count(animationName) == 0)
    {
        LOG(WARNING) << "AnimationModel::playAnimation, tried to play " << animationName << " which does not exist.";
        return;
    }
    m_animationTracks[animationName] = m_animationStash[animationName];
    // If the desired start time has already passed, start now
    const double st = std::max(startTime, m_t);
    m_animationTracks[animationName]->startTime = st;
    m_animationTracks[animationName]->timeScale = timeScale;
}

void
AnimationModel::stopAnimation(std::string animationName, bool stopWhenComplete)
{
    if (!stopWhenComplete)
    {
        // Remove immediately
        if (m_animationTracks.count(animationName) != 0)
        {
            m_animationTracks.erase(animationName);
        }
    }
    else
    {
        // \todo
    }
}

void
AnimationModel::updateDebugSkeleton()
{
    if (m_genDebugSkeleton)
    {
        std::shared_ptr<VecDataArray<double, 3>> verticesPtr = std::make_shared<VecDataArray<double, 3>>();
        std::shared_ptr<VecDataArray<int, 2>>    indicesPtr  = std::make_shared<VecDataArray<int, 2>>();

        VecDataArray<double, 3>& vertices = *verticesPtr;
        VecDataArray<int, 2>&    indices  = *indicesPtr;

        vertices.reserve(static_cast<int>(m_nodes.size() * 2));
        indices.reserve(static_cast<int>(m_nodes.size()));

        int i = 0;
        for (auto node : m_nodes)
        {
            bool hasSkin = false;
            for (auto attachment : node->attachments)
            {
                if (auto meshWeightAttachment = std::dynamic_pointer_cast<MeshWeightAttachment>(attachment))
                {
                    hasSkin = true;
                    break;
                }
            }
            if (hasSkin)
            {
                Vec3d tipPos = Vec3d(0.0, 0.0, 0.0);
                if (node->parent != nullptr)
                {
                    tipPos = (m_globalTransforms[node->parent] * Vec4d(0.0, 0.0, 0.0, 1.0)).head<3>();
                }
                const Vec3d srcPos = (m_globalTransforms[node.get()] * Vec4d(0.0, 0.0, 0.0, 1.0)).head<3>();

                vertices.push_back(srcPos);
                vertices.push_back(tipPos);
                indices.push_back(Vec2i(i, i + 1));
                i += 2;
            }
        }

        vertices.squeeze();
        indices.squeeze();

        m_debugSkeletonGeometry->setInitialVertexPositions(std::make_shared<VecDataArray<double, 3>>(*verticesPtr));
        m_debugSkeletonGeometry->setVertexPositions(verticesPtr);
        m_debugSkeletonGeometry->setLinesIndices(indicesPtr);

        if (!m_debugSkeletonGeometry->hasVertexAttribute("colors"))
        {
            std::shared_ptr<VecDataArray<unsigned char, 3>> scalarsPtr = std::make_shared<VecDataArray<unsigned char, 3>>();
            VecDataArray<unsigned char, 3>&                 scalars    = *scalarsPtr;
            scalars.resize(vertices.size());

            for (int j = 0; j < vertices.size(); j += 2)
            {
                Eigen::Matrix<unsigned char, 3, 1> srcColor;
                srcColor[0] = static_cast<unsigned char>(rand() % 255);
                srcColor[1] = static_cast<unsigned char>(rand() % 255);
                srcColor[2] = static_cast<unsigned char>(rand() % 255);
                scalars[j]  = srcColor;

                Eigen::Matrix<unsigned char, 3, 1> tipColor;
                tipColor[0]    = 255;
                tipColor[1]    = 255;
                tipColor[2]    = 255;
                scalars[j + 1] = tipColor;
            }

            m_debugSkeletonGeometry->setVertexAttribute("colors", scalarsPtr);
            m_debugSkeletonGeometry->setVertexScalars("colors");
        }

        m_debugSkeletonGeometry->postModified();
    }
}
}
