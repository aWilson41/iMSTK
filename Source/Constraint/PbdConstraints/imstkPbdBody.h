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
#include "imstkVecDataArray.h"

#include <unordered_map>

namespace
{
template<typename T>
void
copyOrAllocate(const std::shared_ptr<T>& src, std::shared_ptr<T>& dest)
{
    if (src != nullptr)
    {
        if (dest == nullptr)
        {
            dest = std::make_shared<T>();
        }
        *dest = *src;
    }
}
} // namespace

namespace imstk
{
///
/// \struct PbdBody
///
/// \brief Represents a pbd body in the model. This is a data
/// only object. It does no function.
///
struct PbdBody
{
    public:
        enum class Type
        {
            DEFORMABLE,
            DEFORMABLE_ORIENTED,
            RIGID
        };

    public:
        PbdBody() : bodyHandle(-1) { }
        PbdBody(const int handle) : bodyHandle(handle) { }

        ///
        /// \brief The body should have orientations if its DEFORMABLE_ORIENTED or RIGID
        ///
        bool getOriented() const
        {
            return (bodyType == Type::DEFORMABLE_ORIENTED || bodyType == Type::RIGID);
        }

        ///
        /// \brief Deep copy from src, copying dynamic allocations by value
        ///
        void deepCopy(const PbdBody& src)
        {
            fixedNodeInvMass = src.fixedNodeInvMass;
            bodyHandle       = src.bodyHandle;

            copyOrAllocate(src.prevVertices, prevVertices);
            copyOrAllocate(src.vertices, vertices);
            copyOrAllocate(src.velocities, velocities);
            copyOrAllocate(src.accelerations, accelerations);
            copyOrAllocate(src.masses, masses);
            copyOrAllocate(src.invMasses, invMasses);

            bodyType = src.bodyType;
            if (getOriented())
            {
                copyOrAllocate(src.prevOrientations, prevOrientations);
                copyOrAllocate(src.orientations, orientations);
                copyOrAllocate(src.angularVelocities, angularVelocities);
                copyOrAllocate(src.angularAccel, angularAccel);
                copyOrAllocate(src.inertias, inertias);
                copyOrAllocate(src.invInertias, invInertias);
            }

            fixedNodeIds     = src.fixedNodeIds;
            uniformMassValue = src.uniformMassValue;
        }

    public:
        int bodyHandle; ///< Id in the system

        // For Deformables
        std::shared_ptr<VecDataArray<double, 3>> prevVertices;
        std::shared_ptr<VecDataArray<double, 3>> vertices;

        std::shared_ptr<VecDataArray<double, 3>> velocities;
        std::shared_ptr<VecDataArray<double, 3>> accelerations;

        std::shared_ptr<DataArray<double>> masses;
        std::shared_ptr<DataArray<double>> invMasses;

        // For orientated pbd and rbd
        Type bodyType = Type::DEFORMABLE; // Flag to avoid extra allocation
        std::shared_ptr<StdVectorOfQuatd> prevOrientations;
        std::shared_ptr<StdVectorOfQuatd> orientations;

        std::shared_ptr<VecDataArray<double, 3>> angularVelocities;
        std::shared_ptr<VecDataArray<double, 3>> angularAccel;

        std::shared_ptr<StdVectorOfMat3d> inertias;
        std::shared_ptr<StdVectorOfMat3d> invInertias;

        ///< Nodal/vertex IDs of the nodes that are fixed
        std::vector<int> fixedNodeIds;
        ///< Mass properties, not used if per vertex masses are given in geometry attributes
        double uniformMassValue = 1.0;

        Vec3d initPosTest = Vec3d::Zero();
        Quatd initOrientationTest = Quatd::Identity();
        Mat3d initInertiaTest     = Mat3d::Identity();

        Vec3d externalForce  = Vec3d::Zero();
        Vec3d externalTorque = Vec3d::Zero();

        ///< Map for archiving fixed nodes' mass.
        std::unordered_map<int, double> fixedNodeInvMass;
};

///
/// \brief Provides interface for accessing particles from a 2d array of PbdBody,Particles
///
struct PbdState
{
    public:
        inline Vec3d& getPosition(const std::pair<int, int>& bodyParticleId) const { return (*m_bodies[bodyParticleId.first]->vertices)[bodyParticleId.second]; }
        inline Vec3d& getVelocity(const std::pair<int, int>& bodyParticleId) const { return (*m_bodies[bodyParticleId.first]->velocities)[bodyParticleId.second]; }
        inline Quatd& getOrientation(const std::pair<int, int>& bodyParticleId) const { return (*m_bodies[bodyParticleId.first]->orientations)[bodyParticleId.second]; }
        inline Vec3d& getAngularVelocity(const std::pair<int, int>& bodyParticleId) const { return (*m_bodies[bodyParticleId.first]->angularVelocities)[bodyParticleId.second]; }

        inline double getInvMass(const std::pair<int, int>& bodyParticleId) const { return (*m_bodies[bodyParticleId.first]->invMasses)[bodyParticleId.second]; }
        inline Mat3d& getInvInertia(const std::pair<int, int>& bodyParticleId) const { return (*m_bodies[bodyParticleId.first]->invInertias)[bodyParticleId.second]; }

        std::vector<std::shared_ptr<PbdBody>> m_bodies;
};
} // namespace imstk