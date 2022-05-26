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

#include "EmbeddingConstraint.h"
#include "imstkCollisionUtils.h"

namespace imstk
{
void
EmbeddingConstraint::initConstraint(
    PbdState& bodies,
    const PbdParticleId& ptA1,
    const PbdParticleId& ptB1, const PbdParticleId& ptB2, const PbdParticleId& ptB3,
    Vec3d* p, Vec3d* q)
{
    // Set the triangle
    m_particles[0] = ptB1;
    m_particles[1] = ptB2;
    m_particles[2] = ptB3;
    const Vec3d& x1 = bodies.getPosition(m_particles[0]);
    const Vec3d& x2 = bodies.getPosition(m_particles[1]);
    const Vec3d& x3 = bodies.getPosition(m_particles[2]);

    // Compute intersection point & interpolant on triangle
    CollisionUtils::testSegmentTriangle(*p, *q, x1, x2, x3, m_uvw);
    m_iPt    = x1 * m_uvw[0] + x2 * m_uvw[1] + x3 * m_uvw[2];
    m_iPtVel = Vec3d::Zero();

    m_particles[3] = ptA1;

    // Completely rigid for PBD
    setStiffness(1.0);

    // Compute the interpolant on the line
    {
        m_p = p;
        m_q = q;
        const Vec3d pq = (*p - *q).normalized();
        const Vec3d d  = m_iPt - *q;
        m_t = pq.dot(d);
    }
}

Vec3d
EmbeddingConstraint::computeInterpolantDifference(const PbdState& bodies) const
{
    //const Vec3d& x0 = *m_bodiesFirst[0].vertex;
    const Vec3d& x1 = bodies.getPosition(m_particles[0]);
    const Vec3d& x2 = bodies.getPosition(m_particles[1]);
    const Vec3d& x3 = bodies.getPosition(m_particles[2]);

    Vec3d*      p    = m_p;
    Vec3d*      q    = m_q;
    const Vec3d pq   = (*p - *q);
    const Vec3d pq_n = pq.normalized();

    // Compute the location of the intersection point on both elements
    const Vec3d triPos  = x1 * m_uvw[0] + x2 * m_uvw[1] + x3 * m_uvw[2];
    const Vec3d linePos = (*q) + pq_n * m_t;

    // Compute the transform to align the triangle to the line
    return triPos - linePos;
}

bool
EmbeddingConstraint::computeValueAndGradient(PbdState&           bodies,
                                             double&             c,
                                             std::vector<Vec3d>& n,
                                             std::vector<Vec3d>& r) const
{
    // Triangle
    const Vec3d& x0 = bodies.getPosition(m_particles[0]);
    const Vec3d& x1 = bodies.getPosition(m_particles[1]);
    const Vec3d& x2 = bodies.getPosition(m_particles[2]);
    // Body center of mass
    const Vec3d& x3 = bodies.getPosition(m_particles[3]);

    // Compute the normal/axes of the line
    const Vec3d pq   = *m_p - *m_q;
    const Vec3d pq_n = pq.normalized();

    // Compute the difference between the two interpolated points on the elements
    Vec3d diff = computeInterpolantDifference(bodies);

    // Remove any normal movement (remove only fraction for sort of friction)
    // Frees normal movement
    diff = diff - diff.dot(pq_n) * pq_n * (1.0 - m_normalFriction);
    const Vec3d normal = diff.normalized();

    // Puncture point
    const Vec3d triPos = x0 * m_uvw[0] + x1 * m_uvw[1] + x2 * m_uvw[2];

    n[0] = normal;
    n[1] = normal;
    n[2] = normal;
    n[3] = -normal;
    r[3] = triPos - x3;

    c = -diff.norm() * (1.0 - m_compliance);

    return true;
}
} // namespace imstk