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

#include "imstkImplicitGeometryToSurfaceMeshCD.h"
#include "imstkCollisionData.h"
#include "imstkImplicitGeometry.h"
#include "imstkSurfaceMesh.h"
#include "imstkVecDataArray.h"

#include "imstkImageData.h"
#include "imstkSignedDistanceField.h"

namespace imstk
{
ImplicitGeometryToSurfaceMeshCD::ImplicitGeometryToSurfaceMeshCD(std::shared_ptr<ImplicitGeometry> implicitGeomA,
                                                                 std::shared_ptr<SurfaceMesh>      surfMeshB,
                                                                 std::shared_ptr<CollisionData>    colData) :
    CollisionDetection(CollisionDetection::Type::SurfaceMeshToImplicit, colData),
    m_implicitGeomA(implicitGeomA),
    m_surfMeshB(surfMeshB)
{
    centralGrad.setFunction(m_implicitGeomA);
    if (m_implicitGeomA->getType() == Geometry::Type::SignedDistanceField)
    {
        centralGrad.setDx(std::dynamic_pointer_cast<SignedDistanceField>(m_implicitGeomA)->getImage()->getSpacing());
    }
    else
    {
        // For analytical gradients it doesn't matter as much
        centralGrad.setDx(Vec3d(0.05, 0.05, 0.05));
    }
}

void
ImplicitGeometryToSurfaceMeshCD::computeCollisionData()
{
    m_colData->clearAll();

    // First consider a contact for every single triangle in the surfacemesh
    // Quickly reject by triangle bounding sphere vs implicit tests via centroids
    std::shared_ptr<VecDataArray<double, 3>> verticesPtr = m_surfMeshB->getVertexPositions();
    std::shared_ptr<VecDataArray<int, 3>>    indicesPtr  = m_surfMeshB->getTriangleIndices();
    const VecDataArray<double, 3>&           vertices    = *verticesPtr;
    const VecDataArray<int, 3>&              indices     = *indicesPtr;

    // We try to find the closest point on each triangle to the implicit surface
    // This is done by gradient descent on the surface of a triangle (given the
    // gradient of the SDF over the barycentric coordinates of the triangle)
    // Additionally:
    // We immediately reject those that are not within bounding spheres
    // We terminate when gradient steps are small
    // The starting point for each triangle is computed by choosing the min among the 3 vertices

    std::vector<TriangleFixedVertexCollisionDataElement> tfvColData;
    //std::vector<Vec3d> test;

    // For every triangle
    for (int i = 0; i < indices.size(); i++)
    {
        const int i1 = indices[i][0];
        const int i2 = indices[i][1];
        const int i3 = indices[i][2];

        const Vec3d& v1       = vertices[i1];
        const Vec3d& v2       = vertices[i2];
        const Vec3d& v3       = vertices[i3];
        const Vec3d  centroid = (v1 + v2 + v3) / 3.0;

        // Find the maximal point from centroid for radius
        const double rSqr1   = (centroid - v1).squaredNorm();
        const double rSqr2   = (centroid - v2).squaredNorm();
        const double rSqr3   = (centroid - v3).squaredNorm();
        double       minRSqr = rSqr1;
        if (minRSqr < rSqr2)
        {
            minRSqr = rSqr2;
        }
        if (minRSqr < rSqr3)
        {
            minRSqr = rSqr3;
        }
        const double r = std::sqrt(minRSqr);

        // Get distance from sphere/triangle centroid to the surface
        const double distToSurf = m_implicitGeomA->getFunctionValue(centroid);

        // If within the bounding sphere
        if (distToSurf < r)
        {
            // First compute the starting point for optimization
            // That would be the closest point to the implicit surface
            const double dist1 = m_implicitGeomA->getFunctionValue(v1);
            const double dist2 = m_implicitGeomA->getFunctionValue(v2);
            const double dist3 = m_implicitGeomA->getFunctionValue(v3);

            Vec3d  pos     = v1; // Current iterate point
            double minDist = dist1;
            if (dist2 < minDist)
            {
                minDist = dist2;
                pos     = v2;
            }
            if (dist3 < minDist)
            {
                minDist = dist3;
                pos     = v3;
            }

            int    iter  = 0;
            double minSG = std::numeric_limits<double>::max();
            double sG1   = std::numeric_limits<double>::max();
            double sG2   = std::numeric_limits<double>::max();
            double sG3   = std::numeric_limits<double>::max();
            Vec3d  grad  = Vec3d::Zero();
            Vec3d  s     = Vec3d::Zero();
            // Max iterations of 20, should converge well if close to the surface
            // If far, it will not converge to epsilon, in which case we use max iterations
            // to bail. Due to culling there should be few such cases
            while (minSG > m_epsilon && iter < 100)
            {
                grad = centralGrad(pos);

                // Find support point s
                s     = v1;
                minSG = sG1 = v1.dot(grad);

                sG2 = v2.dot(grad);
                if (sG2 < minSG)
                {
                    minSG = sG2;
                    s     = v2;
                }
                sG3 = v3.dot(grad);
                if (sG3 < minSG)
                {
                    minSG = sG3;
                    s     = v3;
                }

                pos += (2.0 / (2.0 + static_cast<double>(iter))) * (s - pos);
                iter++;
                //printf("Eps: %f\n", minSG);
            }

            // After converging on the closest (deepest) point on the triangle to the implicit surface
            // Check the sign
            const double finalDist = m_implicitGeomA->getFunctionValue(pos);

            // This would imply that some part (the part at pos) of the triangle is inside the surface
            if (finalDist < 0.0)
            {
                // Contact normal
                const Vec3d n = -centralGrad(pos).normalized();
                // Penetration depth
                const double finalAbsDist = std::abs(finalDist);

                // Will this mess up on coincident edges?
                TriangleFixedVertexCollisionDataElement element;
                element.closestDistance = -finalAbsDist;
                element.dir = n;
                //std::cout << "Dist: " << finalAbsDist << std::endl;
                //element.vertexPt = pos;// -n * finalAbsDist;
                element.vertexPt = pos -n * finalAbsDist;
                element.triIdx = i;
                tfvColData.push_back(element);
                //test.push_back(pos);
            }
        }
    }

    // Don't add duplicates to pdColData
    /*for (int i = 0; i < tfvColData.size(); i++)
    {
        for (int j = i + 1; j < tfvColData.size(); j++)
        {
            if (tfvColData[i].vertexPt == tfvColData[j].vertexPt)
            {
                tfvColData.erase(tfvColData.begin() + j);
                j--;
            }
        }
    }*/

    /*for (int i = 0; i < colData.size(); i++)
    {
        m_colData->PDColData.safeAppend(colData[i]);
    }*/
    for (int i = 0; i < tfvColData.size(); i++)
    {
        m_colData->TFVColData.safeAppend(tfvColData[i]);
    }
}
}