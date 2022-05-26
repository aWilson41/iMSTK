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

#include "imstkPbdDistanceConstraint.h"
#include "imstkPbdConstraintFunctor.h"

namespace imstk
{
///
/// \class PbdInflatableDistanceConstraint
///
class PbdInflatableDistanceConstraint : public PbdDistanceConstraint
{
public:
    PbdInflatableDistanceConstraint() : PbdDistanceConstraint() {}

    ///
    /// \brief Initializes the inflatable distance constraint
    ///
    void initConstraint(const Vec3d& p0, const Vec3d& p1,
                        const PbdParticleId& pIdx0,
                        const PbdParticleId& pIdx1,
                        const double k = 1e5)
    {
        PbdDistanceConstraint::initConstraint(p0, p1, pIdx0, pIdx1, k);
        m_initialRestLength = m_restLength;
        m_initialStiffness  = m_stiffness;
    }

    ///
    /// \brief Apply diffusion and update positions
    ///
    void projectConstraint(PbdState&         bodies,
                           const double      dt,
                           const SolverType& type) override
    {
        if (isInflated())
        {
            m_restLength =
                m_restLength + m_diffusionRate * dt * (m_initialRestLength - m_restLength);
            m_stiffness = m_stiffness + m_diffusionRate * dt * (m_initialStiffness - m_stiffness);
        }
        PbdConstraint::projectConstraint(bodies, dt, type);
    }

    void setRestLength(const double restLength) { m_restLength = restLength; }
    const double getRestLength() const { return m_restLength; }

    void setInitRestLength(const double initRestLength) { m_initialRestLength = initRestLength; }
    const double getInitRestLength() const { return m_initialRestLength; }

    ///
    /// \brief Reset constraint rest measurement
    ///
    void resetRestLength() { m_restLength = m_initialRestLength; }
    void resetStiffness() { m_stiffness = m_initialStiffness; }

    bool isInflated() const { return m_inflated; }
    void setInflated() { m_inflated = true; }

public:
    double m_initialRestLength = 0.0;
    double m_initialStiffness  = 1.0;
    double m_diffusionRate     = 0.00001;

    bool m_inflated = false;
};

struct PbdInflatableDistanceConstraintFunctor : public PbdDistanceConstraintFunctor
{
    public:
        PbdInflatableDistanceConstraintFunctor() : PbdDistanceConstraintFunctor() {}
        ~PbdInflatableDistanceConstraintFunctor() override = default;

        ///
        /// \brief Create the distance constraint
        ///
        std::shared_ptr<PbdDistanceConstraint> makeDistConstraint(
            const VecDataArray<double, 3>& vertices,
            int i1, int i2) override
        {
            auto constraint = std::make_shared<PbdInflatableDistanceConstraint>();
            constraint->initConstraint(vertices[i1], vertices[i2],
                { m_bodyIndex, i1 }, { m_bodyIndex, i2 }, m_stiffness);
            return constraint;
        }
};
} // imstk
