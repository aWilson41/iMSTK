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
        PbdBody() : bodyHandle(-1) { }
        PbdBody(const int handle) : bodyHandle(handle) { }

        ///
        /// \brief Deep copy from src, copying dynamic allocations by value
        ///
        void deepCopy(const PbdBody& src)
        {
            fixedNodeInvMass = src.fixedNodeInvMass;
            bodyHandle       = src.bodyHandle;

            if (src.prevVertices != nullptr)
            {
                if (prevVertices == nullptr)
                {
                    prevVertices = std::make_shared<VecDataArray<double, 3>>();
                }
                * prevVertices = *src.prevVertices;
            }
            if (src.vertices != nullptr)
            {
                if (vertices == nullptr)
                {
                    vertices = std::make_shared<VecDataArray<double, 3>>();
                }
                * vertices = *src.vertices;
            }
            if (src.velocities != nullptr)
            {
                if (velocities == nullptr)
                {
                    velocities = std::make_shared<VecDataArray<double, 3>>();
                }
                * velocities = *src.velocities;
            }
            if (src.accelerations != nullptr)
            {
                if (accelerations == nullptr)
                {
                    accelerations = std::make_shared<VecDataArray<double, 3>>();
                }
                * accelerations = *src.accelerations;
            }
            if (src.masses != nullptr)
            {
                if (masses == nullptr)
                {
                    masses = std::make_shared<DataArray<double>>();
                }
                * masses = *src.masses;
            }
            if (src.invMasses != nullptr)
            {
                if (invMasses == nullptr)
                {
                    invMasses = std::make_shared<DataArray<double>>();
                }
                * invMasses = *src.invMasses;
            }

            fixedNodeIds     = src.fixedNodeIds;
            uniformMassValue = src.uniformMassValue;
        }

    protected:
        friend class PbdState;

        ///< Map for archiving fixed nodes' mass.
        std::unordered_map<int, double> fixedNodeInvMass;

    public:
        int bodyHandle; ///< Id in the system

        std::shared_ptr<VecDataArray<double, 3>> prevVertices;
        std::shared_ptr<VecDataArray<double, 3>> vertices;
        std::shared_ptr<VecDataArray<double, 3>> velocities;
        std::shared_ptr<VecDataArray<double, 3>> accelerations;
        std::shared_ptr<DataArray<double>> masses;
        std::shared_ptr<DataArray<double>> invMasses;

        ///< Nodal/vertex IDs of the nodes that are fixed
        std::vector<int> fixedNodeIds;
        ///< Mass properties, not used if per vertex masses are given in geometry attributes
        double uniformMassValue = 1.0;
};
} // namespace imstk