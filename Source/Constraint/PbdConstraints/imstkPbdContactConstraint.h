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

#include "imstkPbdConstraint.h"

#include <iostream>

namespace imstk
{
///
/// \class PbdContactConstraint
///
/// \brief A constraint on a rigid body that defines dtheta through dx applied
/// at a local position r on the body.
///
class PbdContactConstraint : public PbdConstraint
{
public:
    enum class ContactType
    {
        DEFORM,
        RIGID
    };

protected:
    PbdContactConstraint(const int numParticles, std::vector<ContactType> contactTypes) :
        PbdConstraint(numParticles), m_contactTypes(contactTypes)
    {
        CHECK(numParticles == contactTypes.size());
        m_n.resize(numParticles);
        m_r.resize(numParticles);
    }

public:
    ~PbdContactConstraint() override = default;

    ///
    /// \brief Unused
    ///
    bool computeValueAndGradient(PbdState&,
                                 double&, std::vector<Vec3d>&) const override
    {
        return false;
    }

    ///
    /// \brief Compute value and gradient of the constraint
    ///
    virtual bool computeValueAndGradient(PbdState&           bodies,
                                         double&             c,
                                         std::vector<Vec3d>& n,
                                         std::vector<Vec3d>& r) const = 0;

    ///
    /// \brief Update positions by projecting constraints.
    ///
    void projectConstraint(PbdState& bodies,
                           const double dt, const SolverType& type) override;

protected:
    std::vector<Vec3d>       m_n;
    std::vector<Vec3d>       m_r;
    std::vector<ContactType> m_contactTypes;
};

///
/// \class PbdTriangleToBodyConstraint
///
/// \brief Resolves a point defined local to the body to the surface of the triangle
/// Resolves the triangle to the local point
///
class PbdTriangleToBodyConstraint : public PbdContactConstraint
{
public:
    PbdTriangleToBodyConstraint() : PbdContactConstraint(4,
                                                         { ContactType::RIGID, ContactType::DEFORM, ContactType::DEFORM, ContactType::DEFORM })
    {
    }

    ///
    /// \param Body particle
    /// \param Local position on the particle
    /// \param Point0 on triangle
    /// \param Point1 on triangle
    /// \param Point2 on triangle
    ///
    void initConstraint(
        const PbdState& state,
        const PbdParticleId& bodyId,
        const Vec3d contactPtOnBody,
        const PbdParticleId& x0, const PbdParticleId& x1, const PbdParticleId& x2,
        const double stiffness)
    {
        m_particles[0] = bodyId;
        // Compute local, un-transformed position by applying inverse translation and rotation
        m_r[0] = state.getOrientation(bodyId).inverse()._transformVector(contactPtOnBody - state.getPosition(bodyId));
        m_particles[1] = x0;
        m_particles[2] = x1;
        m_particles[3] = x2;

        //setStiffness(stiffness);
        // Infinite stiffness/completely rigid
        m_stiffness  = 0.0;
        m_compliance = 0.0;
    }

    bool computeValueAndGradient(PbdState&           bodies,
                                 double&             c,
                                 std::vector<Vec3d>& n,
                                 std::vector<Vec3d>& r) const override
    {
        const Vec3d& bodyPos = bodies.getPosition(m_particles[0]);
        const Quatd& bodyOrientation = bodies.getOrientation(m_particles[0]);
        const Vec3d& x1 = bodies.getPosition(m_particles[1]);
        const Vec3d& x2 = bodies.getPosition(m_particles[2]);
        const Vec3d& x3 = bodies.getPosition(m_particles[3]);

        // Global position
        const Vec3d p = bodyPos + bodyOrientation._transformVector(r[0]);

        // Compute barycentric coordinates u,v,w
        const Vec3d  v0    = x2 - x1;
        const Vec3d  v1    = x3 - x1;
        const Vec3d  v2    = p - x1;
        const double d00   = v0.dot(v0);
        const double d01   = v0.dot(v1);
        const double d11   = v1.dot(v1);
        const double d20   = v2.dot(v0);
        const double d21   = v2.dot(v1);
        const double denom = d00 * d11 - d01 * d01;
        if (fabs(denom) < 1e-12)
        {
            c = 0.0;
            return false;
        }
        const double v = (d11 * d20 - d01 * d21) / denom;
        const double w = (d00 * d21 - d01 * d20) / denom;
        const double u = 1.0 - v - w;

        // This constraint becomes invalid if moved out of the triangle
        if (u < 0.0 || v < 0.0 || w < 0.0)
        {
            c = 0.0;
            return false;
        }

        // Triangle normal (pointing up on standard counter clockwise triangle)
        const Vec3d normal = v0.cross(v1).normalized();
        // Point could be on either side of triangle, we want to resolve to the triangles plane
        const double l = v2.dot(normal);

        // A
        n[0] = normal;
        // B
        n[1] = -u * normal;
        n[2] = -v * normal;
        n[3] = -w * normal;

        c = l;

        return true;
    }
};

///
/// \class PbdPointToBodyConstraint
///
/// \brief
///
class PbdVertexToBodyConstraint : public PbdContactConstraint
{
public:
    PbdVertexToBodyConstraint() : PbdContactConstraint(2,
                                                       { ContactType::RIGID, ContactType::DEFORM })
    {
    }

    ///
    /// \param Body particle
    /// \param Local position on the particle
    /// \param Point0 on triangle
    /// \param Point1 on triangle
    /// \param Point2 on triangle
    ///
    void initConstraint(
        const PbdState&      state,
        const PbdParticleId& bodyId,
        const Vec3d          contactPtOnBody,
        const PbdParticleId& x0,
        const double         stiffness)
    {
        m_particles[0] = bodyId;
        // Compute local, un-transformed position by applying inverse translation and rotation
        m_r[0] = state.getOrientation(bodyId).inverse()._transformVector(contactPtOnBody - state.getPosition(bodyId));
        m_particles[1] = x0;

        // Infinite stiffness/completely rigid
        m_stiffness  = 0.0;
        m_compliance = 0.0;
    }

    bool computeValueAndGradient(PbdState&           bodies,
                                 double&             c,
                                 std::vector<Vec3d>& n,
                                 std::vector<Vec3d>& r) const override
    {
        const Vec3d& bodyPos = bodies.getPosition(m_particles[0]);
        const Quatd& bodyOrientation = bodies.getOrientation(m_particles[0]);

        // Global position
        const Vec3d p = bodyPos + bodyOrientation._transformVector(r[0]);

        // Current position during solve
        const Vec3d& x1 = bodies.getPosition(m_particles[1]);

        const Vec3d diff = x1 - p;
        c = diff.norm();

        if (c == 0.0)
        {
            return false;
        }

        const Vec3d normal = diff / c;

        // A
        n[0] = -normal;
        // B
        n[1] = normal;

        return true;
    }
};

///
/// \class PbdEdgeToBodyConstraint
///
/// \brief
///
class PbdEdgeToBodyConstraint : public PbdContactConstraint
{
public:
    PbdEdgeToBodyConstraint() : PbdContactConstraint(3,
                                                     { ContactType::RIGID, ContactType::DEFORM, ContactType::DEFORM })
    {
    }

    void initConstraint(
        const PbdState& state,
        const PbdParticleId& bodyId,
        const Vec3d contactPtOnBody,
        const PbdParticleId& x0, const PbdParticleId& x1,
        const double stiffness)
    {
        m_particles[0] = bodyId;
        // Compute local, un-transformed position by applying inverse translation and rotation
        m_r[0] = state.getOrientation(bodyId).inverse()._transformVector(contactPtOnBody - state.getPosition(bodyId));
        m_particles[1] = x0;
        m_particles[2] = x1;

        // Infinite stiffness/completely rigid
        m_stiffness  = 0.0;
        m_compliance = 0.0;
    }

    bool computeValueAndGradient(PbdState&           bodies,
                                 double&             c,
                                 std::vector<Vec3d>& n,
                                 std::vector<Vec3d>& r) const override
    {
        const Vec3d& bodyPos = bodies.getPosition(m_particles[0]);
        const Quatd& bodyOrientation = bodies.getOrientation(m_particles[0]);

        // Global position
        const Vec3d p = bodyPos + bodyOrientation._transformVector(r[0]);

        // Just project p onto x3-x2. Get the normal component for distance to line
        const Vec3d& x1 = bodies.getPosition(m_particles[1]);
        const Vec3d& x2 = bodies.getPosition(m_particles[2]);

        const Vec3d  ab     = x2 - x1;
        const double length = ab.norm();
        if (length == 0.0)
        {
            // There is no distance between the edge, can't do anything
            c = 0.0;
            return false;
        }
        const Vec3d dir1 = ab / length;

        // Project onto the line
        const Vec3d  diff = p - x1;
        const double p1   = dir1.dot(diff);
        if (p1 < 0.0 || p1 > length)
        {
            c = 0.0;
            return false;
        }
        // Remove tangent component to get normal
        const Vec3d  diff1 = diff - p1 * dir1;
        const double l     = diff1.norm();
        if (l == 0.0)
        {
            // The point is on the line
            c = 0.0;
            return false;
        }
        const Vec3d  normal = diff1 / l;
        const double u      = p1 / length;

        // A
        n[0] = normal;
        // B
        n[1] = -(1.0 - u) * normal;
        n[2] = -u * normal;

        c = l;

        return true;
    }
};

///
/// \class PbdBodyToBodyConstraint
///
/// \brief Resolves contact with contact plane reprojection between two bodies
/// Resolves distance between two points on two bodies by given direction
///
class PbdBodyToBodyConstraint : public PbdContactConstraint
{
public:
    PbdBodyToBodyConstraint() : PbdContactConstraint(2,
                                                     { ContactType::RIGID, ContactType::RIGID })
    {
    }

    ///
    /// \param Body particle
    /// \param Local position on the particle
    /// \param Point0 on triangle
    /// \param Point1 on triangle
    /// \param Point2 on triangle
    ///
    void initConstraint(
        const PbdState&      state,
        const PbdParticleId& bodyId0,
        const Vec3d&         contactPtOnBody0,
        const Vec3d&         contactNormal,
        const PbdParticleId& bodyId1,
        const Vec3d          contactPtOnBody1,
        const double         stiffness)
    {
        m_particles[0] = bodyId0;
        // Compute local, un-transformed position by applying inverse translation and rotation
        m_r[0] = state.getOrientation(bodyId0).inverse()._transformVector(contactPtOnBody0 - state.getPosition(bodyId0));
        m_particles[1] = bodyId1;
        m_r[1] = state.getOrientation(bodyId1).inverse()._transformVector(contactPtOnBody1 - state.getPosition(bodyId1));

        m_contactNormal = contactNormal.normalized();

        // Infinite stiffness/completely rigid
        m_stiffness  = 0.0;
        m_compliance = 0.0;
    }

    bool computeValueAndGradient(PbdState&           bodies,
                                 double&             c,
                                 std::vector<Vec3d>& n,
                                 std::vector<Vec3d>& r) const override
    {
        const Vec3d& bodyPos0 = bodies.getPosition(m_particles[0]);
        const Quatd& bodyOrientation0 = bodies.getOrientation(m_particles[0]);
        const Vec3d  p0 = bodyPos0 + bodyOrientation0._transformVector(r[0]);

        const Vec3d& bodyPos1 = bodies.getPosition(m_particles[1]);
        const Quatd& bodyOrientation1 = bodies.getOrientation(m_particles[1]);
        const Vec3d  p1 = bodyPos1 + bodyOrientation1._transformVector(r[1]);

        const Vec3d diff = p1 - p0;

        c = diff.dot(m_contactNormal);

        // A
        n[0] = -m_contactNormal;
        // B
        n[1] = m_contactNormal;

        return true;
    }

protected:
    Vec3d m_contactNormal = Vec3d::Zero();
};

///
/// \class
///
/// \brief Body contact constraint for one way resolution of a local point on body
/// to the given contact plane.
///
class PbdBodyContactNormalConstraint : public PbdContactConstraint
{
public:
    PbdBodyContactNormalConstraint() : PbdContactConstraint(1,
                                                            { ContactType::RIGID })
    {
    }

    ///
    /// \param Body particle
    /// \param Local position on the particle
    /// \param Point0 on triangle
    /// \param Point1 on triangle
    /// \param Point2 on triangle
    ///
    void initConstraint(
        const PbdState&      state,
        const PbdParticleId& bodyId0,
        const Vec3d&         contactPtOnBody0,
        const Vec3d&         contactNormal,
        const Vec3d&         contactPtOnBody2,
        const double         stiffness)
    {
        m_particles[0] = bodyId0;
        // Compute local, un-transformed position by applying inverse translation and rotation
        m_r[0] = state.getOrientation(bodyId0).inverse()._transformVector(contactPtOnBody0 - state.getPosition(bodyId0));

        m_contactNormal = contactNormal.normalized();
        m_contactPt2    = contactPtOnBody2;

        // Infinite stiffness/completely rigid
        m_stiffness  = 0.0;
        m_compliance = 0.0;
    }

    bool computeValueAndGradient(PbdState&           bodies,
                                 double&             c,
                                 std::vector<Vec3d>& n,
                                 std::vector<Vec3d>& r) const override
    {
        const Vec3d& bodyPos0 = bodies.getPosition(m_particles[0]);
        const Quatd& bodyOrientation0 = bodies.getOrientation(m_particles[0]);
        const Vec3d  p0 = bodyPos0 + bodyOrientation0._transformVector(r[0]);

        const Vec3d diff = m_contactPt2 - p0;

        c = diff.dot(m_contactNormal);

        if (c > 0.0)
        {
            c = 0.0;
            return false;
        }

        // A
        n[0] = -m_contactNormal;
        // B
        n[1] = m_contactNormal;

        return true;
    }

protected:
    Vec3d m_contactNormal = Vec3d::Zero();
    Vec3d m_contactPt2    = Vec3d::Zero();
};

///
/// \class PbdBodyToBodyDistanceConstraint
///
/// \brief Constrain two points locally defined on two bodies to maintain
/// a provided distance from each other
///
class PbdBodyToBodyDistanceConstraint : PbdConstraint
{
public:
    PbdBodyToBodyDistanceConstraint() : PbdConstraint()
    {
    }

    ///
    /// \brief ptOnBody are globally defined
    ///
    void initConstraint(
        const PbdState&      state,
        const PbdParticleId& bodyId0,
        const Vec3d          ptOnBody0,
        const PbdParticleId& bodyId1,
        const Vec3d          ptOnBody1,
        const double         restLength,
        const double         stiffness)
    {
        m_particles[0] = bodyId0;
        m_r[0] = state.getOrientation(bodyId0).inverse()._transformVector(ptOnBody0 - state.getPosition(bodyId0));

        m_particles[1] = bodyId1;
        m_r[1] = state.getOrientation(bodyId1).inverse()._transformVector(ptOnBody1 - state.getPosition(bodyId1));

        setStiffness(stiffness);
    }

    bool computeValueAndGradient(PbdState& bodies,
                                 double& c, std::vector<Vec3d>& dcdx)
    {
        // Transform local position to acquire new global (for constraint reprojection)
        const Vec3d& bodyPos0 = bodies.getPosition(m_particles[0]);
        const Quatd& bodyOrientation0 = bodies.getOrientation(m_particles[0]);
        const Vec3d  p0 = bodyPos0 + bodyOrientation0._transformVector(m_r[0]);

        const Vec3d& bodyPos1 = bodies.getPosition(m_particles[1]);
        const Quatd& bodyOrientation1 = bodies.getOrientation(m_particles[1]);
        const Vec3d  p1 = bodyPos1 + bodyOrientation1._transformVector(m_r[0]);

        Vec3d        diff   = p1 - p0;
        const double length = diff.norm();
        if (std::abs(length) < IMSTK_DOUBLE_EPS)
        {
            return false;
        }
        diff /= length;

        // A
        dcdx[0] = diff.normalized();
        // B
        dcdx[1] = -dcdx[0];

        c = length;

        return true;
    }

protected:
    Vec3d m_r[2] = { Vec3d::Zero(), Vec3d::Zero() };
};

// Edge to edge on body?
} // namespace imstk