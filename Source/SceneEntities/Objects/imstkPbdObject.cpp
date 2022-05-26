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

#include "imstkPbdObject.h"
#include "imstkLogger.h"
#include "imstkPbdModel.h"
#include "imstkPointSet.h"

namespace imstk
{
std::shared_ptr<PbdModel>
PbdObject::getPbdModel()
{
    m_pbdModel = std::dynamic_pointer_cast<PbdModel>(m_dynamicalModel);
    return m_pbdModel;
}

void
PbdObject::setBodyFromGeometry()
{
    std::shared_ptr<PbdBody> body = getPbdBody();
    if (body->bodyType == PbdBody::Type::RIGID)
    {
        setRigidBodyFromValues(*body, body->initPosTest, body->initOrientationTest, Vec3d::Zero(), Vec3d::Zero());
    }
    else
    {
        if (m_physicsGeometry != nullptr)
        {
            auto pointSet = std::dynamic_pointer_cast<PointSet>(m_physicsGeometry);
            CHECK(pointSet != nullptr) << "PbdObject " << m_name << " only supports PointSet geometries";
            setDeformBodyFromGeometry(*body, pointSet);
        }
    }

    // Set the geometry to all the functors that need it
    // Not a good solution
    auto functors = m_pbdModel->getConfig()->getFunctors();
    for (auto& functorArray : functors)
    {
        for (auto& functor : functorArray.second)
        {
            if (auto bodyFunctor = std::dynamic_pointer_cast<PbdBodyConstraintFunctor>(functor))
            {
                if (bodyFunctor->m_bodyIndex == body->bodyHandle)
                {
                    auto pointSet = std::dynamic_pointer_cast<PointSet>(m_physicsGeometry);
                    CHECK(pointSet != nullptr) << "PbdObject " << m_name << " only supports PointSet geometries";
                    bodyFunctor->setGeometry(pointSet);
                }
            }
        }
    }
}

bool
PbdObject::initialize()
{
    m_pbdModel = std::dynamic_pointer_cast<PbdModel>(m_dynamicalModel);
    if (m_pbdModel == nullptr)
    {
        LOG(FATAL) << "PbdObject " << m_name << " was not given a PbdModel. Please PbdObject::setDynamicalModel";
        return false;
    }

    setBodyFromGeometry();

    updateGeometries();

    // Sets up maps
    DynamicObject::initialize();

    return true;
}

void
PbdObject::setDynamicalModel(std::shared_ptr<AbstractDynamicalModel> dynaModel)
{
    // todo: If already has another model, should remove the corresponding body?
    m_pbdModel       = std::dynamic_pointer_cast<PbdModel>(dynaModel);
    m_dynamicalModel = dynaModel;

    // If the model already has a pbd body for this PbdObject remove it from
    // that prior model
    if (m_pbdBody != nullptr)
    {
        CHECK(m_pbdModel != nullptr) <<
            "PbdObject has a PbdBody but cannot find associated PbdModel?";
        m_pbdModel->removeBody(m_pbdBody);
    }
    m_pbdBody = m_pbdModel->addBody();
}

void
PbdObject::updatePhysicsGeometry()
{
    DynamicObject::updatePhysicsGeometry();

    if (m_pbdBody->bodyType == PbdBody::Type::RIGID)
    {
        // Pbd has a special case here for rigid body as it cannot apply the
        // Apply the transform back to the geometry
        // If called before body is init'd use initial pose
        if (m_pbdBody->vertices->size() > 0)
        {
            m_physicsGeometry->setTranslation((*m_pbdBody->vertices)[0]);
            m_physicsGeometry->setRotation((*m_pbdBody->orientations)[0]);
        }
        else
        {
            m_physicsGeometry->setTranslation(m_pbdBody->initPosTest);
            m_physicsGeometry->setRotation(m_pbdBody->initOrientationTest);
        }
        m_physicsGeometry->updatePostTransformData();
    }
}

///
/// \brief Convience function. Set the bodyArrPtr to the specified attribute
/// or allocate it with initValue set
///
template<typename T>
static void
setOrAllocate(std::shared_ptr<T>& bodyArrPtr,
              std::shared_ptr<PointSet> pointSet, const std::string& attributeName,
              typename T::ValueType initValue)
{
    // If the input mesh has the attribute already use those
    std::shared_ptr<AbstractDataArray> velocities = pointSet->getVertexAttribute(attributeName);
    if (velocities != nullptr
        && velocities->getNumberOfComponents() == T::NumComponents
        && velocities->getScalarType() == IMSTK_DOUBLE
        && std::dynamic_pointer_cast<T>(velocities)->size() == pointSet->getNumVertices())
    {
        bodyArrPtr = std::dynamic_pointer_cast<T>(velocities);
    }
    // If not, put existing (0 initialized velocities) on mesh
    else
    {
        bodyArrPtr = std::make_shared<T>(pointSet->getNumVertices());
        bodyArrPtr->fill(initValue);
        pointSet->setVertexAttribute(attributeName, bodyArrPtr);
    }
}

void
PbdObject::setDeformBodyFromGeometry(PbdBody& body, std::shared_ptr<PointSet> geom)
{
    body.vertices     = geom->getVertexPositions();
    body.prevVertices = std::make_shared<VecDataArray<double, 3>>(*body.vertices);

    const int numParticles = body.vertices->size();

    // Initialize Mass+InvMass
    {
        // If the input mesh has per vertex masses, use those
        std::shared_ptr<AbstractDataArray> masses = geom->getVertexAttribute("Mass");
        if (masses != nullptr && masses->getNumberOfComponents() == 1
            && masses->getScalarType() == IMSTK_DOUBLE && masses->size() == numParticles)
        {
            body.masses = std::dynamic_pointer_cast<DataArray<double>>(masses);
            body.invMasses->resize(body.masses->size());
            for (int i = 0; i < body.masses->size(); i++)
            {
                (*body.invMasses)[i] = ((*body.masses)[i] == 0.0) ? 0.0 : 1.0 / (*body.masses)[i];
            }
        }
        // If not, initialize as uniform and put on mesh
        else
        {
            // Initialize as uniform
            body.masses    = std::make_shared<DataArray<double>>(numParticles);
            body.invMasses = std::make_shared<DataArray<double>>(numParticles);

            body.masses->fill(body.uniformMassValue);
            body.invMasses->fill((body.uniformMassValue != 0.0) ?
                1.0 / body.uniformMassValue : 0.0);

            geom->setVertexAttribute("Mass", body.masses);
        }
        geom->setVertexAttribute("InvMass", body.invMasses);
    }

    setOrAllocate(body.velocities, geom, "Velocities", Vec3d::Zero());
    setOrAllocate(body.accelerations, geom, "Accelerations", Vec3d::Zero());

    if (body.getOriented())
    {
        // Initialize Inertia + Inv tensors
        {
            // Initialize as uniform
            body.inertias    = std::make_shared<StdVectorOfMat3d>(numParticles);
            body.invInertias = std::make_shared<StdVectorOfMat3d>(numParticles);

            std::fill(body.inertias->begin(), body.inertias->end(), Mat3d::Identity());
            std::fill(body.invInertias->begin(), body.invInertias->end(), Mat3d::Identity());
        }

        // Initialize orientations
        {
            std::shared_ptr<AbstractDataArray> orientations = geom->getVertexAttribute("Orientations");
            if (orientations != nullptr && orientations->getNumberOfComponents() == 4 && orientations->getScalarType() == IMSTK_DOUBLE
                && std::dynamic_pointer_cast<VecDataArray<double, 3>>(orientations)->size() == numParticles)
            {
                auto vec = std::dynamic_pointer_cast<VecDataArray<double, 4>>(orientations);
                body.orientations = std::make_shared<StdVectorOfQuatd>(numParticles);
                for (int i = 0; i < orientations->size(); i++)
                {
                    (*body.orientations)[i] = Quatd((*vec)[i][3], (*vec)[i][1], (*vec)[i][2], (*vec)[i][0]);
                }
            }
            else
            {
                body.orientations = std::make_shared<StdVectorOfQuatd>(numParticles);
                std::fill(body.orientations->begin(), body.orientations->end(), Quatd::Identity());
            }
        }
        body.prevOrientations = std::make_shared<StdVectorOfQuatd>(*body.orientations);

        setOrAllocate(body.angularVelocities, geom, "AngularVelocities", Vec3d::Zero());
        setOrAllocate(body.angularAccel, geom, "AngularAccelerations", Vec3d::Zero());
    }

    // Overwrite some masses for specified fixed points
    body.fixedNodeInvMass = std::unordered_map<int, double>();
    for (auto i : body.fixedNodeIds)
    {
        DataArray<double>& invMasses = *body.invMasses;
        CHECK(i < numParticles && i >= 0) << "Tried to fix particle " << i
                                          << " but there only exist " << numParticles << " particles";
        body.fixedNodeInvMass[i] = invMasses[i];
        invMasses[i] = 0.0;
    }
}

void
PbdObject::setRigidBodyFromValues(PbdBody& body,
                                  const Vec3d& pos, const Quatd& orientation,
                                  const Vec3d& velocity, const Vec3d& angularVelocity)
{
    // Basically a PbdBody with a single particle

    // Currently doesn't support an initial state, zeros it out
    body.vertices     = std::make_shared<VecDataArray<double, 3>>();
    *body.vertices    = { pos };
    body.prevVertices = std::make_shared<VecDataArray<double, 3>>(*body.vertices);

    // Initialize Mass+InvMass
    body.masses    = std::make_shared<DataArray<double>>();
    body.invMasses = std::make_shared<DataArray<double>>();

    *body.masses    = { body.uniformMassValue };
    *body.invMasses = { (body.uniformMassValue != 0.0) ? 1.0 / body.uniformMassValue : 0.0 };

    // Initialize Velocities
    body.velocities  = std::make_shared<VecDataArray<double, 3>>();
    *body.velocities = { velocity };

    // Initialize Accelerations
    body.accelerations  = std::make_shared<VecDataArray<double, 3>>();
    *body.accelerations = { Vec3d::Zero() };

    // Initialize Inertia + Inv tensors
    body.inertias     = std::make_shared<StdVectorOfMat3d>();
    *body.inertias    = { body.initInertiaTest };
    body.invInertias  = std::make_shared<StdVectorOfMat3d>();
    *body.invInertias = { body.initInertiaTest.inverse() };

    // Initialize orientations
    body.orientations     = std::make_shared<StdVectorOfQuatd>();
    *body.orientations    = { orientation };
    body.prevOrientations = std::make_shared<StdVectorOfQuatd>(*body.orientations);

    // Initialize Angular Velocities
    body.angularVelocities  = std::make_shared<VecDataArray<double, 3>>();
    *body.angularVelocities = { angularVelocity };

    // Initialize Angular Accelerations
    body.angularAccel  = std::make_shared<VecDataArray<double, 3>>();
    *body.angularAccel = { Vec3d::Zero() };

    // Overwrite some masses for specified fixed points
    body.fixedNodeInvMass = std::unordered_map<int, double>();
}
} // namespace imstk