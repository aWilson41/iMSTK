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

#include "imstkEdgeEdgeCCDState.h"
#include "imstkPbdEdgeEdgeCCDConstraint.h"
#include "imstkLineMeshToLineMeshCCD.h"

namespace imstk
{
void
PbdEdgeEdgeCCDConstraint::initConstraint(
    const PbdParticleId& prevPtA0, const PbdParticleId& prevPtA1,
    const PbdParticleId& prevPtB0, const PbdParticleId& prevPtB1,
    const PbdParticleId& ptA0, const PbdParticleId& ptA1,
    const PbdParticleId& ptB0, const PbdParticleId& ptB1,
    double stiffnessA, double stiffnessB)
{
    m_particles[0] = prevPtA0;
    m_particles[1] = prevPtA1;
    m_particles[2] = ptA0;
    m_particles[3] = ptA1;

    m_particles[4] = prevPtB0;
    m_particles[5] = prevPtB1;
    m_particles[6] = ptB0;
    m_particles[7] = ptB1;

    m_stiffness[0] = stiffnessA;
    m_stiffness[1] = stiffnessB;
}

bool
PbdEdgeEdgeCCDConstraint::computeValueAndGradient(PbdState& bodies,
                                                  double& c, std::vector<Vec3d>& dcdx) const
{
    const Vec3d& currPt0 = bodies.getPosition(m_particles[2]);
    const Vec3d& currPt1 = bodies.getPosition(m_particles[3]);
    const Vec3d& currPt2 = bodies.getPosition(m_particles[6]);
    const Vec3d& currPt3 = bodies.getPosition(m_particles[7]);

    const Vec3d& prevPt0 = bodies.getPosition(m_particles[0]);
    const Vec3d& prevPt1 = bodies.getPosition(m_particles[1]);
    const Vec3d& prevPt2 = bodies.getPosition(m_particles[4]);
    const Vec3d& prevPt3 = bodies.getPosition(m_particles[5]);

    EdgeEdgeCCDState prevState(prevPt0, prevPt1, prevPt2, prevPt3);
    EdgeEdgeCCDState currState(currPt0, currPt1, currPt2, currPt3);

    double timeOfImpact  = 0.0;
    int    collisionType = EdgeEdgeCCDState::testCollision(prevState, currState, timeOfImpact);
    if (collisionType == 0)
    {
        c = 0.0;
        return false;
    }

    double s  = currState.si();
    double t  = currState.sj();
    Vec3d  n0 = prevState.pi() - prevState.pj();
    Vec3d  n1 = currState.pi() - currState.pj();

    Vec3d n = n1;
    // invert the normal if lines are crossing:
    bool crossing = false;
    if (n0.dot(n1) < 0)
    {
        n *= -1.0;
        crossing = true;
    }

    const double d = n.norm();
    if (d <= 0.0)
    {
        c = 0.0;
        return false;
    }
    n /= d;

    // keep the prev values static by assigning zero vector as solution gradient.
    // This can also be done by assigning invMass as zero for prev timestep vertices.
    dcdx[0] = Vec3d::Zero();
    dcdx[1] = Vec3d::Zero();
    dcdx[4] = Vec3d::Zero();
    dcdx[5] = Vec3d::Zero();

    dcdx[2] = (1 - s) * n;
    dcdx[3] = s * n;

    dcdx[6] = -(1 - t) * n;
    dcdx[7] = -t * n;

    if (crossing)
    {
        c = d + currState.thickness();
    }
    else
    {
        c = std::abs(d - currState.thickness());
    }

    return true;
}
} // namespace imstk