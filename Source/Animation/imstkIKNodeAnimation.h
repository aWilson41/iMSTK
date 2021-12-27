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

#include <unordered_map>

namespace imstk
{
class AnimationNode;

///
/// \class IKNodeAnimation
/// 
/// \brief Provides transforms of a node given time via an IKSystem
/// 
/// 
class IKNodeAnimation : public NodeAnimation
{
public:
    virtual ~IKNodeAnimation() = default;

public:
    virtual void update(const double t) override
    {
        m_transform = (mat4dTranslate(m_nodeTranslation) *
            mat4dRotation(Rotd(m_nodeRotation[0], Vec3d(1.0, 0.0, 0.0))) *
            mat4dRotation(Rotd(m_nodeRotation[1], Vec3d(0.0, 1.0, 0.0))) *
            mat4dRotation(Rotd(m_nodeRotation[2], Vec3d(0.0, 0.0, 1.0))) *
            mat4dScale(m_nodeScale));
    }

    virtual void setNode(std::shared_ptr<AnimationNode> node) override
    {
        NodeAnimation::setNode(node);

        // Record translation, rotation, and scale
    }

protected:
    friend class InverseJacobianIKSystem;
    friend class FABRIKSystem;

    // todo: Init these to the bind transform
    // this would give things like the length of a bone before any solving, which
    // is useful if you're only solving for rotation
    Vec3d m_nodeTranslation;
    Vec3d m_nodeRotation;
    Vec3d m_nodeScale;
};

///
/// \class IKSystem
/// 
/// \brief Solves the orientations of a set of nodes given a target end
/// effector position via the jacobian transpose method
/// \todo: Later split into base class, or generalize with constraints
/// 
class InverseJacobianIKSystem
{
public:
    InverseJacobianIKSystem() = default;
    virtual ~InverseJacobianIKSystem() = default;

public:
    void solve()
    {
        if (m_targetPos == nullptr)
        {
            return;
        }

        // Instead of using the inverse for the system we use the transpose
        // Effectively the goal is to use a set of gradients to advance
        // the rotations until the end effector reaches the target position
        // \todo: Generalize into a solver

        // This is done with simple cross products on (bone pos - end effector)
        // and (end effector - targetPos)

        std::vector<Vec3d> globalPos; // todo
        std::vector<Mat4d> globalTransforms; // todo

        // \todo: Compute end effector
        // The tip/translation component of the last bones transform
        // 
        Vec3d endEffector = Vec3d(0.0, 0.0, 0.0);

        // Calculate the jacobian (transpose)
        Eigen::Matrix<double, -1, 3> jT = Eigen::Matrix<double, -1, 3>();
        jT.resize(nodeAnimations.size() * 3, 3);
        
        for (int i = 0; i < nodeAnimations.size(); i++)
        {
            const Vec3d dist = endEffector - globalPos[i];
            const Mat4d& globalTransform = globalTransforms[i];

            jT.row(i * 3) = globalTransform.block<3, 1>(0, 0).cross(dist);
            jT.row(i * 3 + 1) = globalTransform.block<3, 1>(0, 1).cross(dist);
            jT.row(i * 3 + 2) = globalTransform.block<3, 1>(0, 2).cross(dist);
        }

        Eigen::Matrix<double, 3, -1> j = j.transpose();
        const Vec3d dEndEffector = *m_targetPos - endEffector;
        const double stepSize = (dEndEffector * jT).squaredNorm() * m_stepRatio / (dEndEffector * jT * j).squaredNorm();
        Eigen::Matrix<double, -1, 3> dTheta = dEndEffector * jT * stepSize;

        // Apply solution to bones
        for (int i = 0; i < nodeAnimations.size(); i++)
        {
            Vec3d& eulerRot = nodeAnimations[i]->m_nodeRotation;
            eulerRot[0] += dTheta(i * 3, 0);
            eulerRot[1] += dTheta(i * 3, 1);
            eulerRot[2] += dTheta(i * 3, 2);
        }
    }

    void addNodeAnimation(std::shared_ptr<IKNodeAnimation> node) { nodeAnimations.push_back(node); }

    void setTargetPos(Vec3d* targetPos) { m_targetPos = targetPos; }

    void setNumberOfIterations(const int numIter) { m_numberOfIterations = numIter; }
    const int getNumberOfIterations() const { return m_numberOfIterations; }

protected:
    void computeNodePosition(std::unordered_map<std::shared_ptr<AnimationNode>, Vec3d>& positions)
    {

    }

protected:
    std::vector<std::shared_ptr<IKNodeAnimation>> nodeAnimations;
    Vec3d* m_targetPos = nullptr;
    double m_stepRatio = 1.0;
    int m_numberOfIterations = 5;
};

///
/// \class FABRIKSystem
///
/// \brief Solves the orientations of a set of nodes given a target end
/// effector positions via the "Forward and Backward Reaching IK" method.
/// The current implementation only supports single strands systems. You
/// may setup multiple systems on one skeleton though. 
/// 
class FABRIKSystem
{
public:
    FABRIKSystem(Vec3d& targetPos) : m_targetPos(targetPos)
    {

    }
    virtual ~FABRIKSystem() = default;

public:
    void solve()
    {
    }

protected:
    std::vector<std::shared_ptr<IKNodeAnimation>> nodeAnimations;
    Vec3d& m_targetPos;
    double m_stepRatio = 1.0;
    int m_numberOfIterations = 5;
};
}