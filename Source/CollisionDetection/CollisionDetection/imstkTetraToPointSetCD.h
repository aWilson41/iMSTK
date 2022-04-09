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

#include "imstkCollisionDetectionAlgorithm.h"
#include "imstkMacros.h"
#include "imstkSpatialHashTableSeparateChaining.h"

namespace imstk
{
///
/// \class TetraToPointSetCD
///
/// \brief Computes if points lie in tetrahedrons using spatial hashing
/// Generates tetra-point contact data.
/// By default only generates contact data for both sides.
///
class TetraToPointSetCD : public CollisionDetectionAlgorithm
{
public:
    TetraToPointSetCD();
    ~TetraToPointSetCD() override = default;

    IMSTK_TYPE_NAME(TetraToPointSetCD)

public:
    ///
    /// \brief Compute collision data for both sides simultaneously
    ///
    void computeCollisionDataAB(
        std::shared_ptr<Geometry>      geomA,
        std::shared_ptr<Geometry>      geomB,
        std::vector<CollisionElement>& elementsA,
        std::vector<CollisionElement>& elementsB) override;

protected:
    SpatialHashTableSeparateChaining m_hashTableA; ///> Spatial hash table
    SpatialHashTableSeparateChaining m_hashTableB; ///> Spatial hash table
};
} // namespace imstk