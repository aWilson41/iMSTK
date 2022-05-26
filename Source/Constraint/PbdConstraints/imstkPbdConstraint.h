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

#include "imstkPbdBody.h"

namespace imstk
{
using PbdParticleId = std::pair<int, int>;

///
/// \class PbdConstraint
///
/// \brief Base Constraint class for Position based dynamics constraints
///
class PbdConstraint
{
public:
    ///
    /// \brief Type of solvers
    ///
    enum class SolverType
    {
        xPBD = 0,
        PBD
    };

    PbdConstraint() = default;
    virtual ~PbdConstraint() = default;

    ///
    /// \brief Compute value and gradient of the constraint
    /// \param PbdState provides all the bodies
    /// \param Constraint value
    /// \param Normalized constraint gradients (per vertex)
    ///
    virtual bool computeValueAndGradient(PbdState& bodies,
                                         double& c, std::vector<Vec3d>& dcdx) const = 0;

    ///
    /// \brief Get the vertex indices of the constraint
    ///
    std::vector<PbdParticleId>& getParticles() { return m_particles; }

    ///
    /// \brief Set the tolerance used for pbd constraints
    ///@{
    double getTolerance() const { return m_epsilon; }
    void setTolerance(const double eps) { m_epsilon = eps; }
    ///@}

    ///
    /// \brief Get/Set resitution
    ///@{
    double getRestitution() const { return m_restitution; }
    void setRestitution(const double restitution) { m_restitution = restitution; }
    ///@}

    ///
    /// \brief Get/Set friction
    ///@{
    double getFriction() const { return m_friction; }
    void setFriction(const double friction) { m_friction = friction; }
    ///@}

    ///
    /// \brief Set the stiffness
    ///
    void setStiffness(const double stiffness)
    {
        m_stiffness  = stiffness;
        m_compliance = 1.0 / stiffness;
    }

    ///
    /// \brief  Get the stiffness
    ///
    double getStiffness() const { return m_stiffness; }

    ///
    /// \brief Get the force magnitude, valid after solving lambda
    ///
    double getForce(const double dt) const { return m_lambda / (dt * dt); }

    ///
    /// \brief Use PBD
    ///
    void zeroOutLambda() { m_lambda = 0.0; }

    ///
    /// \brief Update positions by projecting constraints.
    ///
    virtual void projectConstraint(PbdState& bodies, const double dt, const SolverType& type);

    ///
    /// \brief Solve the velocities given to the constraint
    ///
    virtual void correctVelocity(PbdState& bodies);

protected:
    PbdConstraint(const size_t numParticles)
    {
        m_particles.resize(numParticles);
        m_dcdx.resize(numParticles);
    }

    std::vector<PbdParticleId> m_particles; ///< body and index ids per particle

    double m_epsilon        = 1.0e-16;      ///< Tolerance used for the costraints
    double m_stiffness      = 1.0;          ///< used in PBD, [0, 1]
    double m_compliance     = 1e-7;         ///< used in xPBD, inverse of Young's Modulus
    mutable double m_lambda = 0.0;          ///< Lagrange multiplier

    std::vector<Vec3d> m_dcdx;              ///< Normalized constraint gradients (per particle)

    double m_friction    = 0.0;
    double m_restitution = 0.0;
};
} // namespace imstk