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

#include "imstkCollisionDetection.h"
#include "imstkImplicitFunctionFiniteDifferenceFunctor.h"

namespace imstk
{
class ImplicitGeometry;
class SurfaceMesh;
struct CollisionData;

///
/// \class ImplicitGeometryToSurfaceMeshCD
///
/// \brief ImplicitGeometry to SuraceMesh collision detection. Implemented
/// from "Local Optimization for Robust Signed Distance Field Collision" by
/// Miles Macklin et al.
/// In particular it uses the gradient descent method, however, it should be noted
/// that convergence is sensitive to step size
///
class ImplicitGeometryToSurfaceMeshCD : public CollisionDetection
{
public:
    ///
    /// \brief
    /// \param ImplicitGeometry
    /// \param PointSet to test collision with
    /// \param CollisionData to write too
    ///
    ImplicitGeometryToSurfaceMeshCD(std::shared_ptr<ImplicitGeometry> implicitGeomA,
                                    std::shared_ptr<SurfaceMesh>      surfMeshB,
                                    std::shared_ptr<CollisionData>    colData);
    virtual ~ImplicitGeometryToSurfaceMeshCD() override = default;

public:
    const bool getUseSphereCulling() const { return m_useSphereCulling; }
    const double getEpsilon() const { return m_epsilon; }

    void setUseSphereCulling(const bool useSphereCulling) { m_useSphereCulling = useSphereCulling; }
    void setEpsilon(const double epsilon) { m_epsilon = epsilon; }

    ///
    /// \brief Detect collision and compute collision data
    ///
    void computeCollisionData() override;

private:
    std::shared_ptr<ImplicitGeometry> m_implicitGeomA;
    std::shared_ptr<SurfaceMesh>      m_surfMeshB;
    ImplicitFunctionCentralGradient   centralGrad;
    bool   m_useSphereCulling = false;
    double m_epsilon = 0.01;
};
}