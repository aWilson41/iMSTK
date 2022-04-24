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

#include "imstkPbdFemTetConstraint.h"

namespace imstk
{
bool
PbdFemTetConstraint::initConstraint(
    const Vec3d& p0, const Vec3d& p1, const Vec3d& p2, const Vec3d& p3,
    const BodyVertexId& pIdx0, const BodyVertexId& pIdx1,
    const BodyVertexId& pIdx2, const BodyVertexId& pIdx3,
    std::shared_ptr<PbdFemConstraintConfig> config)
{
    m_bodyVertexIds[0] = pIdx0;
    m_bodyVertexIds[1] = pIdx1;
    m_bodyVertexIds[2] = pIdx2;
    m_bodyVertexIds[3] = pIdx3;

    m_elementVolume = (1.0 / 6.0) * (p3 - p0).dot((p1 - p0).cross(p2 - p0));
    m_config     = config;
    m_compliance = 1.0 / (config->m_lambda + 2 * config->m_mu);

    Mat3d m;
    m.col(0) = p0 - p3;
    m.col(1) = p1 - p3;
    m.col(2) = p2 - p3;

    const double det = m.determinant();
    if (fabs(det) > m_epsilon)
    {
        m_invRestMat = m.inverse();
        return true;
    }

    return false;
}

bool
PbdFemTetConstraint::computeValueAndGradient(
    std::vector<PbdBody>& bodies,
    double&               cval,
    std::vector<Vec3d>&   dcdx) const
{
    const BodyVertexId& i0 = m_bodyVertexIds[0];
    const BodyVertexId& i1 = m_bodyVertexIds[1];
    const BodyVertexId& i2 = m_bodyVertexIds[2];
    const BodyVertexId& i3 = m_bodyVertexIds[3];

    const Vec3d& p0 = (*bodies[i0.first].vertices)[i0.second];
    const Vec3d& p1 = (*bodies[i1.first].vertices)[i1.second];
    const Vec3d& p2 = (*bodies[i2.first].vertices)[i2.second];
    const Vec3d& p3 = (*bodies[i3.first].vertices)[i3.second];

    Mat3d m;
    m.col(0) = p0 - p3;
    m.col(1) = p1 - p3;
    m.col(2) = p2 - p3;

    // deformation gradient
    const Mat3d F = m * m_invRestMat;
    // First Piola-Kirchhoff tensor
    Mat3d P;
    // energy constraint
    double C = 0;

    const auto mu     = m_config->m_mu;
    const auto lambda = m_config->m_lambda;

    switch (m_material)
    {
    // P(F) = F*(2*mu*E + lambda*tr(E)*I)
    // E = (F^T*F - I)/2
    case MaterialType::StVK:
    {
        Mat3d E;
        E(0, 0) = 0.5 * (F(0, 0) * F(0, 0) + F(1, 0) * F(1, 0) + F(2, 0) * F(2, 0) - 1.0);                  // xx
        E(1, 1) = 0.5 * (F(0, 1) * F(0, 1) + F(1, 1) * F(1, 1) + F(2, 1) * F(2, 1) - 1.0);                  // yy
        E(2, 2) = 0.5 * (F(0, 2) * F(0, 2) + F(1, 2) * F(1, 2) + F(2, 2) * F(2, 2) - 1.0);                  // zz
        E(0, 1) = 0.5 * (F(0, 0) * F(0, 1) + F(1, 0) * F(1, 1) + F(2, 0) * F(2, 1));                        // xy
        E(0, 2) = 0.5 * (F(0, 0) * F(0, 2) + F(1, 0) * F(1, 2) + F(2, 0) * F(2, 2));                        // xz
        E(1, 2) = 0.5 * (F(0, 1) * F(0, 2) + F(1, 1) * F(1, 2) + F(2, 1) * F(2, 2));                        // yz
        E(1, 0) = E(0, 1);
        E(2, 0) = E(0, 2);
        E(2, 1) = E(1, 2);

        P = 2 * mu * E;
        double tr = E.trace();
        double lt = lambda * tr;
        P(0, 0) += lt;
        P(1, 1) += lt;
        P(2, 2) += lt;
        P        = F * P;

        C = E(0, 0) * E(0, 0) + E(0, 1) * E(0, 1) + E(0, 2) * E(0, 2)
            + E(1, 0) * E(1, 0) + E(1, 1) * E(1, 1) + E(1, 2) * E(1, 2)
            + E(2, 0) * E(2, 0) + E(2, 1) * E(2, 1) + E(2, 2) * E(2, 2);
        C = mu * C + 0.5 * lambda * tr * tr;

        break;
    }

    // P(F) = (2*mu*(F-R) + lambda*(J-1)*J*F^-T
    case MaterialType::Corotation:
    {
        Eigen::JacobiSVD<Mat3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Mat3d                   R = svd.matrixU() * svd.matrixV().adjoint();
        Vec3d                   Sigma(svd.singularValues());
        Mat3d                   invFT = svd.matrixU();
        invFT.col(0) /= Sigma(0);
        invFT.col(1) /= Sigma(1);
        invFT.col(2) /= Sigma(2);
        invFT *= svd.matrixV().adjoint();
        double J  = Sigma(0) * Sigma(1) * Sigma(2);
        Mat3d  FR = F - R;

        P = 2 * mu * FR + lambda * (J - 1) * J * invFT;

        C = FR(0, 0) * FR(0, 0) + FR(0, 1) * FR(0, 1) + FR(0, 2) * FR(0, 2)
            + FR(1, 0) * FR(1, 0) + FR(1, 1) * FR(1, 1) + FR(1, 2) * FR(1, 2)
            + FR(2, 0) * FR(2, 0) + FR(2, 1) * FR(2, 1) + FR(2, 2) * FR(2, 2);
        C = mu * C + 0.5 * lambda * (J - 1) * (J - 1);

        break;
    }
    // P(F) = mu*(F - mu*F^-T) + lambda*log(J)F^-T;
    case MaterialType::NeoHookean:
    {
        // This is a modified NeoHookean to deal with inversions
        // - for which log produces nans
        // - log also gets drastically large as we near 0, so an epsilon is used

        // Modified NeoHookean with flip
        Mat3d  invFT = F.inverse().transpose();
        double det   = F.determinant();
        double logJ  = log(std::abs(det)) * static_cast<double>(!std::signbit(det));

        // Modified NeoHookean with epsilon
        /*Mat3d  invFT = F.inverse().transpose();
        double det = F.determinant();
        double logJ = 0.0;
        if (det < 0.00001)
        {
            logJ = log(0.00001);
        }
        else
        {
            logJ = log(det);
        }*/

        // Original NeoHookean
        /* Mat3d  invFT = F.inverse().transpose();
         double logJ  = log(F.determinant());*/

        P = mu * (F - invFT) + lambda * logJ * invFT;

        C = F(0, 0) * F(0, 0) + F(0, 1) * F(0, 1) + F(0, 2) * F(0, 2)
            + F(1, 0) * F(1, 0) + F(1, 1) * F(1, 1) + F(1, 2) * F(1, 2)
            + F(2, 0) * F(2, 0) + F(2, 1) * F(2, 1) + F(2, 2) * F(2, 2);

        C = 0.5 * mu * (C - 3) - mu * logJ + 0.5 * lambda * logJ * logJ;

        break;
    }

    case MaterialType::Linear:
    {
        break;
    }
    default:
        break;
    }

    Mat3d gradC = m_elementVolume * P * m_invRestMat.transpose();
    cval    = C;
    cval   *=  m_elementVolume;
    dcdx[0] = gradC.col(0);
    dcdx[1] = gradC.col(1);
    dcdx[2] = gradC.col(2);
    dcdx[3] = -dcdx[0] - dcdx[1] - dcdx[2];

    return true;
}
} // namespace imstk