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

#include "imstkDynamicalModel.h"
#include "imstkPbdCollisionConstraint.h"
#include "imstkPbdFemConstraint.h"
#include "imstkPbdState.h"
#include "imstkPbdConstraintFunctor.h"

#include <unordered_map>
#include <unordered_set>

namespace imstk
{
class PointSet;
class PbdConstraintContainer;
class PbdSolver;

///
/// \struct PbdModelConfig
///
/// \brief Gives parameters for PBD simulation
///
struct PbdModelConfig
{
    public:
        ///
        /// \brief Gives the set of standard pbd constraint generation schemes/functors
        /// provided by iMSTK. Note, these do not correspond to constraint types
        /// as there may be multiple schemes for one constraint or even multiple
        /// constraints per scheme
        ///
        enum class ConstraintGenType
        {
            Custom,
            Distance,
            FemTet,
            Volume,
            Area,
            Bend,
            Dihedral,
            ConstantDensity
        };

    public:
        ///
        /// \brief Enables a regular constraint (not FEM constraint) with given stiffness
        /// If constraint of that type already exists, sets the stiffness on it
        ///
        void enableConstraint(ConstraintGenType type, const double stiffness, const int bodyId = 0);

        ///
        /// \brief Enables a bend constraint with given stiffness, stride, and flag for 0 rest length
        /// You may enable multiple with differing strides
        /// If constraint with same stride already exists, updates the stiffness and restLength0 on it
        /// \param Stiffness, how much bend is enforced
        /// \param Stride, distance between vertex connections
        /// \param When true rest length (and angle) are constrained to 0, useful when mesh initial/resting state
        /// is not 0 angled
        ///
        void enableBendConstraint(const double stiffness, const int stride, const bool restLength0 = true, const int bodyId = 0);

        ///
        /// \brief Enable a Fem constraint with the material provided
        ///
        void enableFemConstraint(PbdFemConstraint::MaterialType material, const int bodyId = 0);

        ///
        /// \brief If lame parameters (mu+lambda) are given in femParams, then youngs modulus and poissons ratio are computed
        /// Conversly if youngs and poissons are given, lame parameters are computed
        ///
        void computeElasticConstants();

        ///
        /// \brief Adds a functor to generate constraints
        /// \tparam Must contain operator(PbdConstraintContainer&), could be a PbdConstraintFunctor
        /// or std::function<void(PbdConstraintContainer&)>
        void addPbdConstraintFunctor(std::shared_ptr<PbdConstraintFunctor> functor)
        {
            m_functors[ConstraintGenType::Custom].push_back(functor);
        }

        void addPbdConstraintFunctor(std::function<void(PbdConstraintContainer&)> functor)
        {
            m_functors[ConstraintGenType::Custom].push_back(
                std::make_shared<PbdConstraintFunctorLambda>(functor));
        }

    public:
        double m_viscousDampingCoeff = 0.01;      ///< Viscous damping coefficient [0, 1]
        double m_contactStiffness    = 1.0;       ///< Stiffness for contact
        unsigned int m_iterations    = 10;        ///< Internal constraints pbd solver iterations
        double m_dt = 0.0;                        ///< Time step size
        bool m_doPartitioning = true;             ///< Does graph coloring to solve in parallel

        Vec3d m_gravity = Vec3d(0.0, -9.81, 0.0); ///< Gravity acceleration

        std::shared_ptr<PbdFemConstraintConfig> m_femParams =
            std::make_shared<PbdFemConstraintConfig>(PbdFemConstraintConfig
        {
            0.0,            // Lame constant, if constraint type is FEM
            0.0,            // Lame constant, if constraint type is FEM
            1000.0,         // FEM parameter, if constraint type is FEM
            0.2             // FEM parameter, if constraint type is FEM
            });

        PbdConstraint::SolverType m_solverType = PbdConstraint::SolverType::xPBD;

    protected:
        friend class PbdModel;

        std::unordered_map<ConstraintGenType, std::vector<std::shared_ptr<PbdConstraintFunctor>>> m_functors;
};

///
/// \class PbdModel
///
/// \brief This class implements the position based dynamics model. The
/// PbdModel is a constraint based model that iteratively solves constraints
/// to simulate the dynamics of a body. PbdModel supports SurfaceMesh,
/// LineMesh, or TetrahedralMesh. PointSet is also supported for PBD fluids.
///
/// One of the distinct properties of the PbdModel is that it is first order.
/// This means it simulates dynamics by modifying positions directly. Velocities
/// of the model are computed after positions are solved. Velocities from the
/// previous iteration are applied at the start of the update.
///
/// The PbdModel only takes care of internal body simulation. Collisions
/// are solved in separate systems afterwards to ensure non-penetration.
///
/// References:
/// Matthias Muller, Bruno Heidelberger, Marcus Hennix, and John Ratcliff. 2007. Position based dynamics. J. Vis. Comun. Image Represent. 18, 2 (April 2007), 109-118.
/// Miles Macklin, Matthias Muller, and Nuttapong Chentanez 1. XPBD: position-based simulation of compliant constrained dynamics. In Proc. of Motion in Games. 49-54
///
class PbdModel : public DynamicalModel<PbdState>
{
public:
    PbdModel();
    ~PbdModel() override = default;

    ///
    /// \brief Set simulation parameters
    ///
    void configure(std::shared_ptr<PbdModelConfig> params);

    ///
    /// \brief Add a new PbdBody, return handle to it
    ///
    std::shared_ptr<PbdBody> addPbdBody() { return getCurrentState()->addBody(); }

    void removePbdBody(std::shared_ptr<PbdBody> body) { getCurrentState()->removeBody(body); }

    ///
    /// \brief Get the simulation parameters
    ///
    std::shared_ptr<PbdModelConfig> getConfig() const
    {
        CHECK(m_config != nullptr) << "Cannot PbdModel::getConfig, config is nullptr";
        return m_config;
    }

    ///
    /// \brief Add constraints related to a set of vertices.
    /// \brief Does not check for duplicating pre-existed constraints.
    /// \todo: Move to containers and functors
    ///
    void addConstraints(std::shared_ptr<std::unordered_set<size_t>> vertices);

    void setTimeStep(const double timeStep) override { m_config->m_dt = timeStep; }
    double getTimeStep() const override { return m_config->m_dt; }

    ///
    /// \brief Return all constraints that are solved sequentially
    ///
    std::shared_ptr<PbdConstraintContainer> getConstraints() { return m_constraints; }

    ///
    /// \brief Time integrate the position
    ///@{
    void integratePosition();
    void integratePosition(const PbdBody& body);
    ///@}

    ///
    /// \brief Time integrate the velocity
    ///@{
    void updateVelocity();
    void updateVelocity(const PbdBody& body);
    ///@]

    ///
    /// \brief Solve the internal constraints
    ///
    void solveConstraints();

    ///
    /// \brief Initialize the PBD model
    ///
    bool initialize() override;

    ///
    /// \brief Set the threshold for constraint partitioning
    ///
    void setConstraintPartitionThreshold(size_t threshold) { m_partitionThreshold = threshold; }

    ///
    /// \brief Returns the solver used for internal constraints
    ///
    std::shared_ptr<PbdSolver> getSolver() const { return m_pbdSolver; }

    ///
    /// \brief Sets the solver used for internal constraints
    ///
    void setSolver(std::shared_ptr<PbdSolver> solver) { this->m_pbdSolver = solver; }

    std::shared_ptr<TaskNode> getIntegratePositionNode() const { return m_integrationPositionNode; }

    std::shared_ptr<TaskNode> getSolveNode() const { return m_solveConstraintsNode; }

    std::shared_ptr<TaskNode> getUpdateVelocityNode() const { return m_updateVelocityNode; }

protected:
    ///
    /// \brief Setup the computational graph of Pbd
    ///
    void initGraphEdges(std::shared_ptr<TaskNode> source, std::shared_ptr<TaskNode> sink) override;

    size_t m_partitionThreshold = 16;                      ///< Threshold for constraint partitioning

    std::shared_ptr<PbdSolver> m_pbdSolver = nullptr;      ///< PBD solver

    std::shared_ptr<PbdModelConfig> m_config = nullptr;    ///< Model parameters, must be set before simulation

    std::shared_ptr<PbdConstraintContainer> m_constraints; ///< The set of constraints to update/use

    ///< Computational Nodes
    ///@{
    std::shared_ptr<TaskNode> m_integrationPositionNode = nullptr;
    std::shared_ptr<TaskNode> m_solveConstraintsNode    = nullptr;
    std::shared_ptr<TaskNode> m_updateVelocityNode      = nullptr;
    ///@}
};
} // namespace imstk