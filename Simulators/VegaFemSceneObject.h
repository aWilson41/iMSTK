// This file is part of the SimMedTK project.
// Copyright (c) Center for Modeling, Simulation, and Imaging in Medicine,
//                        Rensselaer Polytechnic Institute
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
//---------------------------------------------------------------------------
//
// Authors:
//
// Contact:
//---------------------------------------------------------------------------

#ifndef SMVEGAFEMSCENEOBJECT_H
#define SMVEGAFEMSCENEOBJECT_H

// SimMedTK includes
#include "Simulators/SceneObjectDeformable.h"

// VEGA includes
#include "centralDifferencesSparse.h"
#include "configFile.h"
#include "corotationalLinearFEMForceModel.h"
#include "corotationalLinearFEM.h"
#include "corotationalLinearFEMMT.h"
#include "eulerSparse.h"
#include "generateMeshGraph.h"
#include "getIntegratorSolver.h"
#include "getopts.h"
#include "graph.h"
#include "implicitBackwardEulerSparse.h"
#include "isotropicHyperelasticFEMForceModel.h"
#include "isotropicHyperelasticFEM.h"
#include "isotropicHyperelasticFEMMT.h"
#include "isotropicMaterial.h"
#include "linearFEMForceModel.h"
#include "loadList.h"
#include "matrixIO.h"
#include "MooneyRivlinIsotropicMaterial.h"
#include "neoHookeanIsotropicMaterial.h"
#include "objMesh.h"
#include "StVKCubeABCD.h"
#include "StVKElementABCDLoader.h"
#include "StVKForceModel.h"
#include "StVKInternalForces.h"
#include "StVKInternalForcesMT.h"
#include "StVKIsotropicMaterial.h"
#include "StVKStiffnessMatrix.h"
#include "StVKStiffnessMatrixMT.h"
#include "StVKTetABCD.h"
#include "StVKTetHighMemoryABCD.h"
#include "tetMesh.h"
#include "volumetricMesh.h"
#include "volumetricMeshLoader.h"
#include "generateSurfaceMesh.h"
#include "generateMassMatrix.h"

class VegaVolumetricMesh;

/// \brief Base class for any scene object that is defmormable
/// and uses FE formulation to compute the evolution of configuration
/// in time.
class VegaFemSceneObject: public SceneObjectDeformable
{
public:

    /// \brief Constructor
    VegaFemSceneObject();

    /// \brief Constructor
    VegaFemSceneObject(const std::shared_ptr<ErrorLog> p_log, const std::string ConfigFile);

    /// \brief Destructor
    ~VegaFemSceneObject();

    /// \brief Initialize the parameters and properties of the simulation object
    void initialize() override;

    /// \brief configure the vega fem scene object using external config file
    bool configure(const std::string ConfigFile) override;

    /// \brief load initial displacements and velocities of the nodes
    void loadInitialStates() override;

    /// \brief reads the fixed nodes from .bou file
    int readBcFromFile(const char* filename, const int offset);

    /// \brief Load the data related to the vertices that will be fixed
    void loadFixedBC() override;

    /// \brief Load volume meshes
    void loadVolumeMesh() override;

    /// \brief Load the surface mesh
    void loadSurfaceMesh() override;

    /// \brief Forces as a result of user interaction
    /// (through an interface such as mouse or haptic device)
    /// with the scene during runtime are added here
    void applyUserInteractionForces() override;

    /// \brief Use the computed displacemetnt update
    /// to interpolate to the secondary display mesh
    void updateSecondaryRenderingMesh() override;

    /// \brief print object specific data
    void printInfo() const override;

    /// \brief Update the deformations by time stepping
    void advanceDynamics() override;

    /// \brief Advance in time by a specificed amount and a chosen time stepping scheme
    inline void advanceOneTimeStep();

    /// \brief rest the object to inital configuration and reset initial states
    void resetToInitialState() override;

    /// \brief Set the type of formulation used to model the deformation
    void setDeformableModelType();

    /// \brief load the scripted externalloads
    void loadScriptedExternalForces();

    /// \brief Create the force model (underlying formulation)
    void createForceModel();

    /// \brief Inititialize the time integrator
    void initializeTimeIntegrator();

    /// \brief Forces that are defined by the user
    ///  before the start of the simulation
    ///  is added to the external force vector here
    inline void applyScriptedExternalForces();

    /// \brief Updates the stats related to timing, fps etc.
    /// Also updates window title with real-time information
    inline void updatePerformanceMetrics();

    /// \brief check all the surface nodes for the closest node within
    /// certain threshold and set it to be the pulled vertex
    void setPulledVertex(const core::Vec3d &userPos);

    /// \brief returns velocity given the
    /// localtion in the global velocity vector
    core::Vec3d getVelocity(const int dofID) const;

    /// \brief returns displacement given the
    /// localtion in the global displacement vector
    core::Vec3d getDisplacementOfNodeWithDofID(const int dofID) const;

    /// \brief returns acceleration given the
    /// localtion in the global acceleration vector
    core::Vec3d getAccelerationOfNodeWithDofID(const int dofID) const;

    /// \brief returns the number of nodes
    int getNumNodes() const;

    /// \brief returns total degree of freeedom (including fixed)
    int getNumTotalDof() const;

    /// \brief returns unknown degree of freeedom
    int getNumDof() const;

    /// \brief returns number of nodes that are fixed in space
    int getNumFixedNodes() const;

    /// \brief returns the degree od freedom that are known or fixed
    int getNumFixedDof() const;

    /// \brief serialize function explicity writes the object to the memory block
    ///each scene object should know how to write itself to a memory block
    virtual void serialize(void */*p_memoryBlock*/) override {};

    /// \brief Unserialize function can recover the object from the memory location
    virtual void unSerialize(void */*p_memoryBlock*/) override {};

    /// \brief this function may not be used
    ///every Scene Object should know how to clone itself.
    /// Since the data structures will be
    ///in the beginning of the modules(such as simulator, viewer, collision etc.)
    //virtual std::shared_ptr<SceneObject> clone() override { return nullptr; };

    /// \brief not implemented yet.
    std::shared_ptr<SceneObject> clone() override;

    std::shared_ptr<VegaVolumetricMesh> getVolumetricMesh();
    void setVolumetricMesh(std::shared_ptr<VegaVolumetricMesh> mesh)
    {
        this->volumetricMesh = mesh;
    }

private:

    int staticSolver;
    int graphicFrame;
    int explosionFlag; ///< 1 if the simulation goes unstable
    int positiveDefinite; ///< 1 if the effective matrix is positive definite

    VegaPerformanceCounter performaceTracker;

    std::shared_ptr<VegaObjectConfig> femConfig;

    // Time integrators
    std::shared_ptr<IntegratorBase> integratorBase; ///< integrator
    std::shared_ptr<ImplicitNewmarkSparse> implicitNewmarkSparse;
    std::shared_ptr<IntegratorBaseSparse> integratorBaseSparse;

    // Force models
    std::shared_ptr<ForceModel> forceModel; ///< Type of formulation driving the FEM simulation
    std::shared_ptr<StVKInternalForces> stVKInternalForces;
    std::shared_ptr<StVKStiffnessMatrix> stVKStiffnessMatrix;

    // Volume meshes and related graphs
    std::shared_ptr<VegaVolumetricMesh> volumetricMesh; ///< volume mesh

    // Sparse matrices
    std::shared_ptr<SparseMatrix> massMatrix; ///< sparse mass matrix need for FEM simulation
    std::shared_ptr<SparseMatrix> LaplacianDampingMatrix; ///< sparse damping matrix need for FEM simulation

    std::shared_ptr<LinearSolver> linearSolver;
};

#endif //SMVEGAFEMSCENEOBJECT_H
