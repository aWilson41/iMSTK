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

#include "imstkAttachment.h"
#include "imstkAnimationNode.h"
#include "imstkPointSet.h"
#include "imstkVecDataArray.h"

namespace imstk
{

void
MeshWeightAttachment::update(const Mat4d& pose)
{
    // The source vertices are already deformed, ie: if a bone bind pose says its 45 degrees rotated
    // the vertices will already be in bind pose. Applying 45 degrees rotation again would double up.

    // So we need the boneOffset. Which undeforms/inverts them such that pose would be all identity.
    // Then the pose can be applied
    const Mat4d transform = pose * parent->boneOffset;

    Vec3d t, s;
    Mat3d r;
    mat4dTRS(transform, t, r, s);
    const bool isNegative = (std::signbit(s[0]) && std::signbit(s[1]) && std::signbit(s[2])); // If all negative

    std::shared_ptr<VecDataArray<double, 3>> preVerticesPtr = geometry->getVertexPositions(Geometry::DataType::PreTransform);
    const VecDataArray<double, 3>&           preVertices    = *preVerticesPtr;

    std::shared_ptr<VecDataArray<double, 3>> postVerticesPtr = geometry->getVertexPositions(Geometry::DataType::PostTransform);
    VecDataArray<double, 3>&                 postVertices    = *postVerticesPtr;

    const bool hasNormals = geometry->hasVertexAttribute("preNormals");
    if (hasNormals)
    {
        std::shared_ptr<VecDataArray<double, 3>> preNormalsPtr = std::dynamic_pointer_cast<VecDataArray<double, 3>>(geometry->getVertexAttribute("preNormals"));
        const VecDataArray<double, 3>& preNormals = *preNormalsPtr;

        std::shared_ptr<VecDataArray<double, 3>> postNormalsPtr = geometry->getVertexNormals();
        VecDataArray<double, 3>& postNormals = *postNormalsPtr;

        for (size_t i = 0; i < weights.size(); i++)
        {
            VertexWeight& vertexWeight = weights[i];
            const int     vertexId = vertexWeight.vertexId;

            // PrePos is bind pose vertex in root space
            const Vec3d& prePos = preVertices[vertexId];
            postVertices[vertexId] += (transform * Vec4d(prePos[0], prePos[1], prePos[2], 1.0)).head<3>() * vertexWeight.weight;

            // Commonly done and usable but technically not correct, norm happens later
            const Vec3d& preNormal = preNormals[vertexId];
            if (!isNegative)
            {
                postNormals[vertexId] -= (r * preNormal) * vertexWeight.weight;
            }
            else
            {
                postNormals[vertexId] += (r * preNormal) * vertexWeight.weight;
            }
        }
    }
    else
    {
        for (size_t i = 0; i < weights.size(); i++)
        {
            VertexWeight& vertexWeight = weights[i];
            const int     vertexId = vertexWeight.vertexId;

            // PrePos is bind pose vertex in root space
            const Vec3d& prePos = preVertices[vertexId];
            postVertices[vertexId] += (transform * Vec4d(prePos[0], prePos[1], prePos[2], 1.0)).head<3>() * vertexWeight.weight;
        }
    }
}

void
MeshAttachment::update(const Mat4d& pose)
{
    geometry->setTransform(pose);
}
}