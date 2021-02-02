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

#include "imstkPbdPointTriangleConstraint.h"

#include <iostream>

namespace imstk
{
void
PbdPointTriangleConstraint::initConstraint(const Side side, Vec3d* ptA1, double* invMassA1,
    Vec3d* ptB1, double* invMassB1, Vec3d* ptB2, double* invMassB2, Vec3d* ptB3, double* invMassB3,
    std::shared_ptr<PbdCollisionConstraintConfig> configA,
    std::shared_ptr<PbdCollisionConstraintConfig> configB,
    Vec3d dir, double depth)
{
    m_side = side;

    m_bodiesFirst[0].vertex = ptA1;
    m_bodiesFirst[0].invMass = invMassA1;

    m_bodiesSecond[0].vertex = ptB1;
    m_bodiesSecond[0].invMass = invMassB1;
    m_bodiesSecond[1].vertex = ptB2;
    m_bodiesSecond[1].invMass = invMassB2;
    m_bodiesSecond[2].vertex = ptB3;
    m_bodiesSecond[2].invMass = invMassB3;

    m_configA = configA;
    m_configB = configB;

    m_dir = dir;
    m_depth = depth;
}

bool
PbdPointTriangleConstraint::computeValueAndGradient(double& c,
                                                    VecDataArray<double, 3>& dcdxA,
                                                    VecDataArray<double, 3>& dcdxB) const
{
    const Vec3d& x0 = *m_bodiesFirst[0].vertex;
    const Vec3d& x1 = *m_bodiesSecond[0].vertex;
    const Vec3d& x2 = *m_bodiesSecond[1].vertex;
    const Vec3d& x3 = *m_bodiesSecond[2].vertex;

    Vec3d        v0 = x2 - x1;
    Vec3d        v1 = x3 - x1;
    Vec3d        v2 = x0 - x1;
    const double d00 = v0.dot(v0);
    const double d01 = v0.dot(v1);
    const double d11 = v1.dot(v1);
    const double d20 = v2.dot(v0);
    const double d21 = v2.dot(v1);
    const double denom = d00 * d11 - d01 * d01;
    const double v = (d11 * d20 - d01 * d21) / denom;
    const double w = (d00 * d21 - d01 * d20) / denom;
    const double u = 1.0 - v - w;

    const Vec3d n = v0.cross(v1).normalized();
    const double l = v2.dot(n);

    //std::cout << "x0: " << x0[0] << ", " << x0[1] << ", " << x0[2] << std::endl;
    //std::cout << "n: " << n[0] << ", " << n[1] << ", " << n[2] << std::endl;
    //std::cout << "l: " << l << std::endl;

    //if (v < 0.0 || w < 0.0 || (v + w) > 1.0)
    //{
    //    //LOG(WARNING) << "Projection point not inside the triangle";
    //    c = 0.0;
    //    return false;
    //}

    if (m_side == Side::AB)
    {
        const double dist = m_configA->m_proximity + m_configB->m_proximity;
        if (l > dist)
        {
            c = 0.0;
            return false;
        }

        dcdxA[0] = n;
        dcdxB[0] = -u * n;
        dcdxB[1] = -v * n;
        dcdxB[2] = -w * n;

        c = l - dist;
    }
    else if (m_side == Side::A)
    {
        const double dist = m_configA->m_proximity;
        if (l > dist)
        {
            c = 0.0;
            return false;
        }

        dcdxA[0] = n;
        dcdxB[0] = Vec3d::Zero();
        dcdxB[1] = Vec3d::Zero();
        dcdxB[2] = Vec3d::Zero();

        c = l - dist;
    }
    else
    {
        const double dist = m_configB->m_proximity;
        /*if (l > dist)
        {
            c = 0.0;
            return false;
        }*/

        dcdxA[0] = Vec3d::Zero();

        // Weight n out over the 3 verts (u,v,w sum to 1)
        dcdxB[0] = -u * n;
        dcdxB[1] = -v * n;
        dcdxB[2] = -w * n;

        c = l - dist;
    }




    //const Vec3d& x0 = *m_bodiesFirst[0].vertex;
    //const Vec3d& x1 = *m_bodiesSecond[0].vertex;
    //const Vec3d& x2 = *m_bodiesSecond[1].vertex;
    //const Vec3d& x3 = *m_bodiesSecond[2].vertex;

    //Vec3d        v0 = x2 - x1;
    //Vec3d        v1 = x3 - x1;
    //Vec3d        v2 = x0 - x1;
    //const double d00 = v0.dot(v0);
    //const double d01 = v0.dot(v1);
    //const double d11 = v1.dot(v1);
    //const double d20 = v2.dot(v0);
    //const double d21 = v2.dot(v1);
    //const double denom = d00 * d11 - d01 * d01;
    //const double v = (d11 * d20 - d01 * d21) / denom;
    //const double w = (d00 * d21 - d01 * d20) / denom;
    //const double u = 1.0 - v - w;

    //const double l = v2.dot(m_dir);

    ////std::cout << "x0: " << x0[0] << ", " << x0[1] << ", " << x0[2] << std::endl;
    ////std::cout << "n: " << n[0] << ", " << n[1] << ", " << n[2] << std::endl;
    ////std::cout << "l: " << l << std::endl;

    ////if (v < 0.0 || w < 0.0 || (v + w) > 1.01)
    ////{
    ////    //LOG(WARNING) << "Projection point not inside the triangle";
    ////    c = 0.0;
    ////    return false;
    ////}

    ///*if (m_side == Side::AB)
    //{
    //    const double dist = m_configA->m_proximity + m_configB->m_proximity;
    //    if (l > dist)
    //    {
    //        c = 0.0;
    //        return false;
    //    }

    //    dcdxA[0] = n;
    //    dcdxB[0] = -u * n;
    //    dcdxB[1] = -v * n;
    //    dcdxB[2] = -w * n;

    //    c = l - dist;
    //}
    //else */if (m_side == Side::A)
    //{
    //    const double dist = m_configA->m_proximity;
    //    if (l > dist)
    //    {
    //        c = 0.0;
    //        return false;
    //    }

    //    dcdxA[0] = m_dir;
    //    dcdxB[0] = Vec3d::Zero();
    //    dcdxB[1] = Vec3d::Zero();
    //    dcdxB[2] = Vec3d::Zero();

    //    c = l - dist;
    //}
    //else
    //{
    //    const double dist = m_configB->m_proximity;
    //    /*if (l > dist)
    //    {
    //        c = 0.0;
    //        return false;
    //    }*/

    //    dcdxA[0] = Vec3d::Zero();

    //    // Weight n out over the 3 verts (u,v,w sum to 1)
    //    dcdxB[0] = -u * m_dir;
    //    dcdxB[1] = -v * m_dir;
    //    dcdxB[2] = -w * m_dir;

    //    std::cout << "depth: " << m_depth << std::endl;
    //    c = m_depth;
    //}

    return true;
}
} // imstk
