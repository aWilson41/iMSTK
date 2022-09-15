/*
** This file is part of the Interactive Medical Simulation Toolkit (iMSTK)
** iMSTK is distributed under the Apache License, Version 2.0.
** See accompanying NOTICE for details.
*/

#include "imstkCamera.h"
#include "imstkCapsule.h"
#include "imstkDeviceManager.h"
#include "imstkDeviceManagerFactory.h"
#include "imstkDirectionalLight.h"
#include "imstkIsometricMap.h"
#include "imstkKeyboardDeviceClient.h"
#include "imstkMeshIO.h"
#include "imstkObjectControllerGhost.h"
#include "imstkPbdContactConstraint.h"
#include "imstkPbdModel.h"
#include "imstkPbdModelConfig.h"
#include "imstkPbdObject.h"
#include "imstkPbdObjectCollision.h"
#include "imstkPbdObjectController.h"
#include "imstkPbdObjectGrasping.h"
#include "imstkPointwiseMap.h"
#include "imstkPortHoleInteraction.h"
#include "imstkRenderMaterial.h"
#include "imstkScene.h"
#include "imstkSceneManager.h"
#include "imstkSignedDistanceField.h"
#include "imstkSimulationManager.h"
#include "imstkSimulationUtils.h"
#include "imstkSphere.h"
#include "imstkTetrahedralize.h"
#include "imstkVisualModel.h"
#include "imstkVTKViewer.h"
#include "imstkDummyClient.h"

using namespace imstk;

//#define USE_TWO_HAPTIC_DEVICES

const double stomachCapsuleCompliance = 0.00001;

const double capsuleMass = 1.0;
const double capsuleInertiaScale = 0.01;
const double capsuleDistConstraintStiffness      = 100.0;
const double capsuleHingeConstraintCompliance    = 0.000001;
const double capsuleHingeConstraintEndCompliance = 0.001; // The two end points of the capsule chain

const double tissueParticleMass      = 0.1;
const double tissueDistStiffness     = 10000.0;
const double tissueDihedralStiffness = 1.0;

const double globalLinearDamping  = 0.1;
const double globalAngularDamping = 0.03;

//const double globalLinearDamping = 0.075;
//const double globalAngularDamping = 0.01;

const int    iterations        = 1;     // Prefer lower timestep with XPBD
const double simTimestep       = 0.005; // 0.0025
const double executionTimestep = 0.005;

const double graspStiffness = 0.01;

const double tetYoungsModulus = 4000.0;
const double tetPossionsRatio = 0.47;

const Vec3d leftPortHolePos   = Vec3d(-0.065, 0.59, -0.951);
const Vec3d rightPortHolePos  = Vec3d(0.016, 0.597, -0.939);
const Vec3d cameraPortHolePos = Vec3d(-0.003, 0.589, -0.923);

const Vec3d focalPoint = Vec3d(0.03, 0.569, -1.1); // Can be used to quickly find the esophagus/diaphgram junction

const double lapToolCapsuleLength = 0.3;

// There a couple major things happening in this example
// - We have a deformable stomach
// - We a have a diaphgram SDF which helps hold the stomach in place
// - We have 3 ligaments attached to the diaphgram which help hold it in place.
//  In particular the gastrosplenic ligament helps hold the stomach back when it
//  is pulled to the left.
// - We have 3 tools. 2 controlled by devices. 'f' can switch the camera tool. These
// can grasp the stomach to wrap around the esophagus.
// - We have a string of connected capsule running down the esophagus and right side of stomach
// to avoid mesh-mesh collision and keep the shape of the esophagus.
// - We have 2 tools, though 3 is desirable, this allows wrapping of the stomach around
//  the esophagus. The most critical component in wrapping.
//

///
/// \brief Defines a custom way to add Fem Constraint so we add ones of high stiffness
/// around the esophagus, but ones of low stiffness on the stomach
///
struct PbdFemTetConstraintFunctorCustom : public PbdFemTetConstraintFunctor
{
    public:
        PbdFemTetConstraintFunctorCustom() = default;
        ~PbdFemTetConstraintFunctorCustom() override = default;

        void operator()(PbdConstraintContainer& constraints) override
        {
            // Check for correct mesh type
            CHECK(std::dynamic_pointer_cast<TetrahedralMesh>(m_geom) != nullptr)
                << "PbdFemTetConstraint can only be generated with a TetrahedralMesh";

            // Create constraints
            auto                                     tetMesh     = std::dynamic_pointer_cast<TetrahedralMesh>(m_geom);
            std::shared_ptr<VecDataArray<double, 3>> verticesPtr = m_geom->getVertexPositions();
            const VecDataArray<double, 3>&           vertices    = *verticesPtr;
            std::shared_ptr<VecDataArray<int, 4>>    elementsPtr = tetMesh->getCells();
            const VecDataArray<int, 4>&              elements    = *elementsPtr;
            const DataArray<unsigned char>&          mask = *std::dynamic_pointer_cast<DataArray<unsigned char>>(tetMesh->getVertexAttribute("mask"));

            ParallelUtils::parallelFor(elements.size(),
                [&](const size_t k)
                {
                    const Vec4i& tet = elements[k];

                    if (mask[k] == 0)
                    {
                        auto c = std::make_shared<PbdFemTetConstraint>(m_matTypeStomach);
                        c->initConstraint(
                        vertices[tet[0]], vertices[tet[1]], vertices[tet[2]], vertices[tet[3]],
                            { m_bodyIndex, tet[0] }, { m_bodyIndex, tet[1] },
                            { m_bodyIndex, tet[2] }, { m_bodyIndex, tet[3] },
                        m_femConfigStomach);
                        constraints.addConstraint(c);
                    }
                    else
                    {
                        auto c = std::make_shared<PbdFemTetConstraint>(m_matTypeEsophagus);
                        c->initConstraint(
                        vertices[tet[0]], vertices[tet[1]], vertices[tet[2]], vertices[tet[3]],
                            { m_bodyIndex, tet[0] }, { m_bodyIndex, tet[1] },
                            { m_bodyIndex, tet[2] }, { m_bodyIndex, tet[3] },
                        m_femConfigEsophagus);
                        constraints.addConstraint(c);
                    }
            }, elements.size() > 100);
        }

        void setMaterialTypeEsophagus(const PbdFemTetConstraint::MaterialType materialType) { m_matTypeEsophagus = materialType; }
        PbdFemTetConstraint::MaterialType getMaterialTypeEsophagus() const { return m_matTypeEsophagus; }

        void setFemConfigEsophagus(std::shared_ptr<PbdFemConstraintConfig> femConfig) { m_femConfigEsophagus = femConfig; }
        std::shared_ptr<PbdFemConstraintConfig> getFemConfigEsophagus() const { return m_femConfigEsophagus; }

        void setMaterialTypeStomach(const PbdFemTetConstraint::MaterialType materialType) { m_matTypeStomach = materialType; }
        PbdFemTetConstraint::MaterialType getMaterialTypeStomach() const { return m_matTypeStomach; }

        void setFemConfigStomach(std::shared_ptr<PbdFemConstraintConfig> femConfig) { m_femConfigStomach = femConfig; }
        std::shared_ptr<PbdFemConstraintConfig> getFemConfigStomach() const { return m_femConfigStomach; }

    protected:
        PbdFemTetConstraint::MaterialType m_matTypeEsophagus = PbdFemTetConstraint::MaterialType::StVK;
        std::shared_ptr<PbdFemConstraintConfig> m_femConfigEsophagus = nullptr;

        PbdFemTetConstraint::MaterialType m_matTypeStomach = PbdFemTetConstraint::MaterialType::StVK;
        std::shared_ptr<PbdFemConstraintConfig> m_femConfigStomach = nullptr;
};

///
/// \brief Compute the a vector of points indices that are coincident with the on the parent
///
static std::vector<int>
computeFixedPtsViaMap(std::shared_ptr<PointSet> parent,
                      std::shared_ptr<PointSet> child,
                      const double              tolerance = 0.00001)
{
    std::vector<int> fixedPts;

    auto map = std::make_shared<PointwiseMap>();
    map->setParentGeometry(parent);
    map->setChildGeometry(child);
    map->setTolerance(tolerance);
    map->compute();
    fixedPts.reserve(child->getNumVertices());
    for (int i = 0; i < child->getNumVertices(); i++)
    {
        fixedPts.push_back(map->getParentVertexId(i));
    }
    return fixedPts;
}

///
/// \brief Given two tet meshes, compute a mask array for the parent that tells
/// which tets are from the child. Assigned specified value
///
static std::shared_ptr<DataArray<unsigned char>>
computeTetMask(std::shared_ptr<TetrahedralMesh> parentMesh,
               std::shared_ptr<TetrahedralMesh> childMesh,
               std::shared_ptr<DataArray<unsigned char>> maskPtr, const unsigned char value)
{
    int foundCount = 0;
    for (int i = 0; i < childMesh->getNumCells(); i++)
    {
        const Vec4i& esoTet  = (*childMesh->getCells())[i];
        const Vec3d& a0      = (*childMesh->getVertexPositions())[esoTet[0]];
        const Vec3d& a1      = (*childMesh->getVertexPositions())[esoTet[1]];
        const Vec3d& a2      = (*childMesh->getVertexPositions())[esoTet[2]];
        const Vec3d& a3      = (*childMesh->getVertexPositions())[esoTet[3]];
        const Vec3d  centerA = (a0 + a1 + a2 + a3) * 0.25;

        // For every tet of the esophagus find corresponding tet in tissue
        for (int j = 0; j < parentMesh->getNumCells(); j++)
        {
            const Vec4i& tissueTet = (*parentMesh->getCells())[j];
            const Vec3d& b0      = (*parentMesh->getVertexPositions())[tissueTet[0]];
            const Vec3d& b1      = (*parentMesh->getVertexPositions())[tissueTet[1]];
            const Vec3d& b2      = (*parentMesh->getVertexPositions())[tissueTet[2]];
            const Vec3d& b3      = (*parentMesh->getVertexPositions())[tissueTet[3]];
            const Vec3d  centerB = (b0 + b1 + b2 + b3) * 0.25;

            // Compare all arrangements of vertices. Can we assume widing is the same?
            if (centerA.isApprox(centerB, 0.00001))
            {
                (*maskPtr)[j] = 1;
                foundCount++;
            }
        }
    }
    LOG(INFO) << "computeTetMask found " << foundCount << " / " << childMesh->getNumCells() << " corresponding tets, there are ";
    return maskPtr;
}

///
/// \brief Creates tissue object
/// \param name
/// \param model dynamical model with which this is simulated
///
static std::shared_ptr<PbdObject>
makeTissueObj(const std::string&        name,
              std::shared_ptr<PbdModel> model)
{
    std::shared_ptr<TetrahedralMesh> tissueMesh =
        MeshIO::read<TetrahedralMesh>(VLAHHS_DATA_PATH "/hiatal_rotated/stomach_tet.msh");
    std::shared_ptr<TetrahedralMesh> esophagusMesh =
        MeshIO::read<TetrahedralMesh>(VLAHHS_DATA_PATH "/hiatal_rotated/stomach_tet_esophagus.vtk");
    std::shared_ptr<SurfaceMesh> tissueSurfMesh = tissueMesh->extractSurfaceMesh();

    // Compute which tets are part of the esophagus, notate this in the mask with 1
    auto maskPtr = std::make_shared<DataArray<unsigned char>>(tissueMesh->getNumCells());
    maskPtr->fill(0);
    computeTetMask(tissueMesh, esophagusMesh, maskPtr, 1);
    tissueMesh->setVertexAttribute("mask", maskPtr);

    // Setup the material
    auto material = std::make_shared<RenderMaterial>();
    material->setColor(Color(93.0 / 255.0, 38.0 / 255.0, 37.0 / 255.0));
    material->setDisplayMode(RenderMaterial::DisplayMode::WireframeSurface);
    material->setBackFaceCulling(false);
    //material->setOpacity(0.5);

    // Setup the Object
    auto tissueObj = std::make_shared<PbdObject>(name);
    tissueObj->setVisualGeometry(tissueSurfMesh);
    auto colAndVisualMap = std::make_shared<PointwiseMap>(tissueMesh, tissueSurfMesh);
    colAndVisualMap->setTolerance(0.0001);
    tissueObj->setPhysicsToCollidingMap(colAndVisualMap);
    tissueObj->setCollidingGeometry(tissueSurfMesh);
    tissueObj->getVisualModel(0)->setRenderMaterial(material);
    tissueObj->setPhysicsGeometry(tissueMesh);

    tissueObj->setDynamicalModel(model);
    tissueObj->getPbdBody()->uniformMassValue = tissueParticleMass;

    auto         fixedVerts = MeshIO::read<PointSet>(VLAHHS_DATA_PATH "/hiatal_rotated/stomach_fixed_verts.obj");
    PointwiseMap fixedMapper;
    fixedMapper.setTolerance(0.00001);
    fixedMapper.setParentGeometry(tissueObj->getPhysicsGeometry());
    fixedMapper.setChildGeometry(fixedVerts);
    fixedMapper.compute();
    const std::unordered_map<int, int>& fixedMap = fixedMapper.getMap();
    for (auto i : fixedMap)
    {
        tissueObj->getPbdBody()->fixedNodeIds.push_back(i.second);
    }

    auto customFunctor = std::make_shared<PbdFemTetConstraintFunctorCustom>();
    customFunctor->setBodyIndex(tissueObj->getPbdBody()->bodyHandle);
    const double possionsRatio = tetPossionsRatio;
    const double stomachYoungsModulus   = tetYoungsModulus * 0.75;
    const double esophagusYoungsModulus = tetYoungsModulus * 3.0;
    auto         femConfigStomach       = std::make_shared<PbdFemConstraintConfig>(
        stomachYoungsModulus / 2.0 / (1.0 + possionsRatio), stomachYoungsModulus * possionsRatio / ((1.0 + possionsRatio) * (1.0 - 2.0 * possionsRatio)),
        stomachYoungsModulus, possionsRatio);
    customFunctor->setFemConfigStomach(femConfigStomach);
    auto femConfigEsophagus = std::make_shared<PbdFemConstraintConfig>(
        esophagusYoungsModulus / 2.0 / (1.0 + possionsRatio), esophagusYoungsModulus * possionsRatio / ((1.0 + possionsRatio) * (1.0 - 2.0 * possionsRatio)),
        esophagusYoungsModulus, possionsRatio);
    customFunctor->setFemConfigEsophagus(femConfigEsophagus);
    customFunctor->setMaterialTypeEsophagus(PbdFemConstraint::MaterialType::StVK);     // Cheaper less large deformation
    customFunctor->setMaterialTypeStomach(PbdFemConstraint::MaterialType::NeoHookean); // Large deformation
    model->getConfig()->addPbdConstraintFunctor(customFunctor);

    return tissueObj;
}

///
/// \brief Creates and PbdObject that links/connects to the tissueObj provided
/// via computation of overlapping points
///
static std::shared_ptr<PbdObject>
makeLigamentObj(std::string                name,
                std::string                ligamentGeomFileName,
                std::string                fixedVertexFileName,
                const double               stiffness,
                std::shared_ptr<PbdObject> tissueObj)
{
    auto ligamentObj = std::make_shared<PbdObject>(name);

    auto surfMesh = MeshIO::read<SurfaceMesh>(ligamentGeomFileName);

    std::shared_ptr<PbdModel> model = tissueObj->getPbdModel();
    ligamentObj->setVisualGeometry(surfMesh);
    ligamentObj->getVisualModel(0)->getRenderMaterial()->setColor(
        Color(183.0 / 255.0, 86.0 / 255.0, 49.0 / 255.0));
    ligamentObj->getVisualModel(0)->getRenderMaterial()->setBackFaceCulling(false);
    //ligamentObj->getVisualModel(0)->getRenderMaterial()->setOpacity(0.5);
    //ligamentObj->getVisualModel(0)->getRenderMaterial()->setDisplayMode(RenderMaterial::DisplayMode::WireframeSurface);
    ligamentObj->setPhysicsGeometry(surfMesh);
    ligamentObj->setCollidingGeometry(surfMesh);
    ligamentObj->setDynamicalModel(model);
    ligamentObj->getPbdBody()->uniformMassValue = tissueParticleMass;
    auto fixedPtMesh = MeshIO::read<LineMesh>(fixedVertexFileName);
    ligamentObj->getPbdBody()->fixedNodeIds = computeFixedPtsViaMap(surfMesh, fixedPtMesh);

    model->getConfig()->enableConstraint(PbdModelConfig::ConstraintGenType::Distance, tissueDistStiffness,
        ligamentObj->getPbdBody()->bodyHandle);
    model->getConfig()->enableConstraint(PbdModelConfig::ConstraintGenType::Dihedral, tissueDihedralStiffness,
        ligamentObj->getPbdBody()->bodyHandle);
    model->getConfig()->addPbdConstraintFunctor([ = ](PbdConstraintContainer& container)
        {
            // Setup the map between the tissueObj and current object
            auto tissueToLigamentMap = std::make_shared<PointwiseMap>();
            tissueToLigamentMap->setParentGeometry(tissueObj->getPhysicsGeometry());
            tissueToLigamentMap->setChildGeometry(ligamentObj->getPhysicsGeometry());
            tissueToLigamentMap->setTolerance(0.001);
            std::unordered_map<int, int> tissueMap;
            tissueToLigamentMap->computeMap(tissueMap);

            for (auto i : tissueMap)
            {
                const PbdParticleId id1 = { tissueObj->getPbdBody()->bodyHandle, i.second };
                const PbdParticleId id2 = { ligamentObj->getPbdBody()->bodyHandle, i.first };

                auto constraint = std::make_shared<PbdDistanceConstraint>();
                constraint->initConstraint(
                    tissueObj->getPbdModel()->getBodies().getPosition(id1),
                    tissueObj->getPbdModel()->getBodies().getPosition(id2),
                    id1, id2, 100000.0);
                container.addConstraint(constraint);
            }
        });

    return ligamentObj;
}

void
makeCapsules(const std::string&                       name,
             std::vector<std::shared_ptr<PbdObject>>& capsuleObjs,
             std::shared_ptr<PbdModel>                model,
             const double                             capsuleRadius)
{
    capsuleObjs.clear();

    // Read in a medial line
    auto maLineMesh =
        MeshIO::read<LineMesh>(VLAHHS_DATA_PATH "/hiatal_rotated/capsules.obj");
    std::shared_ptr<VecDataArray<double, 3>> maVerticesPtr = maLineMesh->getVertexPositions();
    const VecDataArray<double, 3>&           maVertices    = *maVerticesPtr;
    std::shared_ptr<VecDataArray<int, 2>>    maIndicesPtr  = maLineMesh->getCells();
    const VecDataArray<int, 2>&              maIndices     = *maIndicesPtr;

    // For every line segment in the line mesh create a rigid body capsule
    for (int i = 0; i < maIndices.size(); i++)
    {
        const Vec2i& cell   = maIndices[i];
        const Vec3d  p      = maVertices[cell[0]];
        const Vec3d  q      = maVertices[cell[1]];
        const Vec3d  center = (q + p) * 0.5;
        const Vec3d  diff   = q - p;
        const Vec3d  dir    = diff.normalized();

        auto capsule = std::make_shared<Capsule>(Vec3d::Zero(), capsuleRadius, diff.norm());
        //auto capsule = std::make_shared<Sphere>(Vec3d::Zero(), capsuleRadius);

        auto capsuleObj = std::make_shared<PbdObject>(name);
        capsuleObj->setDynamicalModel(model);
        capsuleObj->setVisualGeometry(capsule);
        capsuleObj->setPhysicsGeometry(capsule);
        capsuleObj->setCollidingGeometry(capsule);
        //capsuleObj->getVisualModel(0)->getRenderMaterial()->setColor(Color::Red);
        //capsuleObj->getVisualModel(0)->getRenderMaterial()->setOpacity(0.5);

        // Setup body
        const Quatd orientation = Quatd::FromTwoVectors(Vec3d(0.0, 1.0, 0.0), dir);
        capsuleObj->getPbdBody()->setRigid(center,    // Position
            capsuleMass,                              // Mass
            orientation,                              // Orientation
            Mat3d::Identity() * capsuleInertiaScale); // Inertia
        model->getConfig()->setBodyDamping(capsuleObj->getPbdBody()->bodyHandle, 0.1, 0.0);

        capsuleObjs.push_back(capsuleObj);
    }

    // This one just connects every single capsule to be fixed in place in positin
    // and orientation
    model->getConfig()->addPbdConstraintFunctor(
        [ = ](PbdConstraintContainer& container)
        {
            // Also add capsule->capsule joint constraints to keep them connected
            for (int i = 0; i < maIndices.size() - 1; i++)
            {
                // Assumes consistent directions of the capsules
                auto capsule0 = std::dynamic_pointer_cast<Capsule>(capsuleObjs[i]->getPhysicsGeometry());
                Vec3d jp0     = Vec3d::Zero();
                {
                    const Vec3d& capsulePos         = capsule0->getPosition();
                    const double capsuleLength      = capsule0->getLength();
                    const Quatd& capsuleOrientation = capsule0->getOrientation();
                    const Vec3d& capsulePosA        = capsulePos - 0.5 * capsuleLength * capsuleOrientation.toRotationMatrix().col(1);
                    const Vec3d& capsulePosB        = capsulePos + (capsulePos - capsulePosA);
                    jp0 = capsulePosB;
                }
                auto capsule1 = std::dynamic_pointer_cast<Capsule>(capsuleObjs[i + 1]->getPhysicsGeometry());
                Vec3d jp1     = Vec3d::Zero();
                {
                    const Vec3d& capsulePos         = capsule1->getPosition();
                    const double capsuleLength      = capsule1->getLength();
                    const Quatd& capsuleOrientation = capsule1->getOrientation();
                    const Vec3d& capsulePosA        = capsulePos - 0.5 * capsuleLength * capsuleOrientation.toRotationMatrix().col(1);
                    const Vec3d& capsulePosB        = capsulePos + (capsulePos - capsulePosA);
                    jp1 = capsulePosA;
                }
                const Vec3d jointPos = (jp0 + jp1) * 0.5;

                auto bodyToBodyConstraint = std::make_shared<PbdBodyToBodyDistanceConstraint>();
                bodyToBodyConstraint->initConstraint(model->getBodies(),
                    { capsuleObjs[i]->getPbdBody()->bodyHandle, 0 },
                    jointPos,
                    { capsuleObjs[i + 1]->getPbdBody()->bodyHandle, 0 },
                    jointPos, capsuleHingeConstraintCompliance);
                container.addConstraint(bodyToBodyConstraint);
            }
            // Fix the two ends
            {
                // Assumes consistent directions of the capsules
                auto capsule0 = std::dynamic_pointer_cast<Capsule>(capsuleObjs[0]->getPhysicsGeometry());
                Vec3d jp0     = Vec3d::Zero();
                {
                    const Vec3d& capsulePos         = capsule0->getPosition();
                    const double capsuleLength      = capsule0->getLength();
                    const Quatd& capsuleOrientation = capsule0->getOrientation();
                    const Vec3d& capsulePosA        = capsulePos - 0.5 * capsuleLength * capsuleOrientation.toRotationMatrix().col(1);
                    const Vec3d& capsulePosB        = capsulePos + (capsulePos - capsulePosA);
                    jp0 = capsulePosA;
                }
                // Assumes consistent directions of the capsules
                auto capsule1 = std::dynamic_pointer_cast<Capsule>(capsuleObjs.back()->getPhysicsGeometry());
                Vec3d jp1     = Vec3d::Zero();
                {
                    const Vec3d& capsulePos         = capsule1->getPosition();
                    const double capsuleLength      = capsule1->getLength();
                    const Quatd& capsuleOrientation = capsule1->getOrientation();
                    const Vec3d& capsulePosA        = capsulePos - 0.5 * capsuleLength * capsuleOrientation.toRotationMatrix().col(1);
                    const Vec3d& capsulePosB        = capsulePos + (capsulePos - capsulePosA);
                    jp1 = capsulePosB;
                }

                const PbdParticleId& fixedVid0 = model->addVirtualParticle(jp0, 0.0, Vec3d::Zero(), true);
                auto bodyToBodyConstraint0     = std::make_shared<PbdBodyToBodyDistanceConstraint>();
                bodyToBodyConstraint0->initConstraint(model->getBodies(),
                    { capsuleObjs[0]->getPbdBody()->bodyHandle, 0 },
                    jp0,
                    fixedVid0,
                    0.0, capsuleHingeConstraintEndCompliance);
                container.addConstraint(bodyToBodyConstraint0);

                const PbdParticleId& fixedVid1 = model->addVirtualParticle(jp1, 0.0, Vec3d::Zero(), true);
                auto bodyToBodyConstraint1     = std::make_shared<PbdBodyToBodyDistanceConstraint>();
                bodyToBodyConstraint1->initConstraint(model->getBodies(),
                    { capsuleObjs.back()->getPbdBody()->bodyHandle, 0 },
                    jp1,
                    fixedVid1,
                    0.0, capsuleHingeConstraintEndCompliance);
                container.addConstraint(bodyToBodyConstraint1);
            }
        });
}

std::shared_ptr<PbdObject>
makeLapToolObj(const std::string&        name,
               std::shared_ptr<PbdModel> model)
{
    auto lapTool = std::make_shared<PbdObject>(name);

    auto toolGeom = std::make_shared<Capsule>(
        Vec3d(0.0, 0.0, lapToolCapsuleLength * 0.5 - 0.005), // Position
        0.002,                                               // Radius
        lapToolCapsuleLength,                                // Length
        Quatd(Rotd(PI_2, Vec3d(1.0, 0.0, 0.0))));            // Orientation

    const double lapToolHeadLength = 0.005;
    auto         graspCapsule      = std::make_shared<Capsule>(
        Vec3d(0.0, 0.0, lapToolHeadLength * 0.5),                           // Position
        0.0035,                                                             // Radius
        lapToolHeadLength,                                                  // Length
        Quatd::FromTwoVectors(Vec3d(0.0, 1.0, 0.0), Vec3d(0.0, 0.0, 1.0))); // Orientation

    auto lapToolVisualGeom = MeshIO::read<SurfaceMesh>(iMSTK_DATA_ROOT "/Surgical Instruments/LapTool/laptool_all_in_one.obj");

    lapTool->setDynamicalModel(model);
    lapTool->setPhysicsGeometry(toolGeom);
    lapTool->setCollidingGeometry(toolGeom);
    lapTool->setVisualGeometry(lapToolVisualGeom);
    lapTool->setPhysicsToVisualMap(std::make_shared<IsometricMap>(toolGeom, lapToolVisualGeom));

    // Add grasp capsule for visualization
    auto graspVisualModel = std::make_shared<VisualModel>();
    graspVisualModel->setGeometry(graspCapsule);
    graspVisualModel->getRenderMaterial()->setIsDynamicMesh(false);
    graspVisualModel->setIsVisible(true);
    lapTool->addVisualModel(graspVisualModel);

    std::shared_ptr<RenderMaterial> material = lapTool->getVisualModel(0)->getRenderMaterial();
    material->setIsDynamicMesh(false);
    material->setMetalness(1.0);
    material->setRoughness(0.2);
    material->setShadingModel(RenderMaterial::ShadingModel::PBR);

    lapTool->getPbdBody()->setRigid(
        Vec3d(0.0, 0.0, lapToolCapsuleLength * 0.5) + Vec3d(0.0, 0.1, -1.0),
        5.0,
        Quatd::Identity(),
        Mat3d::Identity() * 0.08);
    lapTool->getPbdBody()->hasGravity = false;

    auto controller = lapTool->addComponent<PbdObjectController>();
    controller->setControlledObject(lapTool);
    controller->setLinearKs(5000.0);
    controller->setAngularKs(1.0);
    //controller->setForceScaling(0.00001);
    controller->setForceScaling(0.0);
    controller->setSmoothingKernelSize(15);
    controller->setUseForceSmoothening(true);

    auto ghostController = lapTool->addComponent<ObjectControllerGhost>();
    ghostController->setController(controller);

    // The center of mass of the object is at the tip this allows most force applied
    // to the tool at the tip upon touch to be translated into linear force. Suitable
    // for 3dof devices.
    //
    // However, the point at which you actually apply force is on the back of the tool,
    // this is important for the inversion of control in lap tools (right movement at the
    // back should move the tip left).
    controller->setHapticOffset(Vec3d(0.0, 0.0, lapToolCapsuleLength));

    // \todo: The grasp capsule and its map can't be placed as components yet.
    // For now grasp capsule is placed as a VisualModel and the map updated in this lambda
    auto graspCapsuleMap    = std::make_shared<IsometricMap>(toolGeom, graspCapsule);
    auto graspCapsuleUpdate = lapTool->addComponent<LambdaBehaviour>("graspCapsuleUpdate");
    graspCapsuleUpdate->setUpdate([ = ](const double&)
        {
            graspCapsuleMap->update();
        });

    return lapTool;
}

///
/// \brief Visual represents the diaphgram but also serves as collision container
/// for stomach
///
std::shared_ptr<CollidingObject>
makeDiaphgramObj(const std::string& name)
{
    auto diaphgramObj = std::make_shared<CollidingObject>(name);

    // \todo: Generate this from geometry and cache it, use vtk for hashing
    auto diaphgramSdfImage =
        MeshIO::read<ImageData>(VLAHHS_DATA_PATH "/hiatal_rotated/stomach_container_SDF.mhd")
        ->cast(IMSTK_DOUBLE);
    diaphgramObj->setCollidingGeometry(std::make_shared<SignedDistanceField>(diaphgramSdfImage));

    /*imstkNew<SurfaceMeshFlyingEdges> isoExtract;
    isoExtract->setInputImage(diaphgramSdfImage);
    isoExtract->update();
    diaphgramObj->setVisualGeometry(isoExtract->getOutputMesh());
    diaphgramObj->getVisualModel(0)->getRenderMaterial()->setOpacity(0.5);*/

    auto surfMesh = MeshIO::read<SurfaceMesh>(VLAHHS_DATA_PATH "/hiatal_rotated/diaphragm.obj");
    diaphgramObj->setVisualGeometry(surfMesh);

    auto material = std::make_shared<RenderMaterial>();
    material->setMetalness(0.0);
    material->setRoughness(0.26);
    material->setNormalStrength(2.0);
    material->setOcclusionStrength(0.0);
    material->setCoatRoughness(0.1);
    material->setCoatStrength(1.0);
    material->setCoatColor(Color::White);
    material->setCoatIOR(5.0);
    material->setBaseIOR(5.0);
    material->setCoatNormalScale(0.5);
    material->setEdgeTint(Color::White);
    material->setIsDynamicMesh(false);
    material->addTexture(std::make_shared<Texture>(VLAHHS_DATA_PATH "/textures/Diaphragm_SG1.png",
        Texture::Type::Diffuse));
    material->addTexture(std::make_shared<Texture>(
        VLAHHS_DATA_PATH "/textures/Diaphragm_SG1_normal.png", Texture::Type::Normal));
    material->addTexture(std::make_shared<Texture>(
        VLAHHS_DATA_PATH "/textures/Diaphragm_SG1_normal.png", Texture::Type::CoatNormal));
    material->setShadingModel(RenderMaterial::ShadingModel::PBR);
    diaphgramObj->getVisualModel(0)->setRenderMaterial(material);

    return diaphgramObj;
}

int
main()
{
    // Write log to stdout and file
    Logger::startLogger();

    // Setup a scene
    auto scene = std::make_shared<Scene>("Hernia");
    // Camera setup and copied for debug, camera is overwritten later by tool
    scene->getActiveCamera()->setFocalPoint(focalPoint);
    scene->getActiveCamera()->setPosition((cameraPortHolePos + focalPoint) * 0.5);
    const Vec3d camDiff = focalPoint - cameraPortHolePos;
    scene->getActiveCamera()->setViewUp(-camDiff.cross(Vec3d(1.0, 0.0, 0.0)).normalized());
    scene->getActiveCamera()->setFieldOfView(45.0);  // 45deg fov for this lap tool camera
    scene->getConfig()->debugCamBoundingBox = false; // Don't init debug camera to bounding box of scene
    *scene->getCamera("debug") = *scene->getActiveCamera();

    auto bodyObject = std::make_shared<SceneObject>("body");
    {
        auto surfMesh = MeshIO::read<SurfaceMesh>(VLAHHS_DATA_PATH "/hiatal_rotated/body.obj");
        bodyObject->setVisualGeometry(surfMesh);
        bodyObject->getVisualModel(0)->getRenderMaterial()->setShadingModel(
            RenderMaterial::ShadingModel::PBR);
        std::shared_ptr<RenderMaterial> material =
            bodyObject->getVisualModel(0)->getRenderMaterial();
        material->setRoughness(0.8);
        material->setMetalness(0.1);
        material->setOpacity(0.5);
    }
    scene->addSceneObject(bodyObject);

    auto                            pbdModel  = std::make_shared<PbdModel>();
    std::shared_ptr<PbdModelConfig> pbdConfig = pbdModel->getConfig();
    //pbdConfig->m_gravity = Vec3d(0.0, 0.0, 0.0);
    pbdConfig->m_gravity    = Vec3d(0.0, -0.1, 0.0);
    pbdConfig->m_dt         = simTimestep;
    pbdConfig->m_iterations = iterations;
    pbdConfig->m_linearDampingCoeff  = globalLinearDamping;
    pbdConfig->m_angularDampingCoeff = globalAngularDamping;
    pbdConfig->m_doPartitioning      = false;

    std::shared_ptr<PbdObject> stomachObj = makeTissueObj("stomach", pbdModel);
    scene->addSceneObject(stomachObj);

    std::shared_ptr<PbdObject> leftToolObj = makeLapToolObj("leftToolObj", pbdModel);
    scene->addSceneObject(leftToolObj);
    std::shared_ptr<PbdObject> rightToolObj = makeLapToolObj("rightToolObj", pbdModel);
    scene->addSceneObject(rightToolObj);
    std::shared_ptr<PbdObject> cameraToolObj = makeLapToolObj("cameraToolObj", pbdModel);
    (*cameraToolObj->getPbdBody()->vertices)[0]     = Vec3d(0.001, 0.078, -1.206);
    (*cameraToolObj->getPbdBody()->orientations)[0] = Quatd(0.826, -0.54, -0.149, 0.061);
    scene->addSceneObject(cameraToolObj);

    std::vector<std::shared_ptr<PbdObject>> esophagusCapsuleObjs;
    makeCapsules("capsuleObj", esophagusCapsuleObjs, pbdModel, 0.004);
    for (int i = 0; i < esophagusCapsuleObjs.size(); i++)
    {
        scene->addSceneObject(esophagusCapsuleObjs[i]);
        auto capsCollision = std::make_shared<PbdObjectCollision>(
            stomachObj, esophagusCapsuleObjs[i]);
        capsCollision->setRestitution(0.0);
        capsCollision->setFriction(0.0);
        capsCollision->setUseCorrectVelocity(false);
        capsCollision->setRigidBodyCompliance(stomachCapsuleCompliance);
        scene->addInteraction(capsCollision);
    }

    std::shared_ptr<CollidingObject> diaphgramObj = makeDiaphgramObj("diaphgram");
    scene->addSceneObject(diaphgramObj);

    #pragma region Ligaments
    /* std::shared_ptr<PbdObject> fatObj = makeLigamentObj("fat_obj",
         "C:/Users/Andx_/Desktop/NewHernia/fat.obj",
         "C:/Users/Andx_/Desktop/NewHernia/fat_fixed.obj",
         1000.0,
         tissueObj);
     scene->addSceneObject(fatObj);*/

    /*std::shared_ptr<PbdObject> lesserOmentumObj = makeLigamentObj("lesser_omentum",
        HERNIA_RESOURCE_DIR "/lesser_omentum_sim.obj",
        HERNIA_RESOURCE_DIR "/lesser_omentum_sim_fixed.obj",
        1000.0,
        tissueObj);
    scene->addSceneObject(lesserOmentumObj);

    std::shared_ptr<PbdObject> greaterOmentumObj =
        makeLigamentObj("greater_omentum",
            HERNIA_RESOURCE_DIR "/greater_omentum_sim.obj",
            HERNIA_RESOURCE_DIR "/greater_omentum_sim_fixed.obj",
            1000.0,
            tissueObj);
    scene->addSceneObject(greaterOmentumObj);

    std::shared_ptr<PbdObject> gastrosplenicLigamentObj =
        makeLigamentObj("gastrosplenic_ligament",
            HERNIA_RESOURCE_DIR "/FundoTest1/gastrosplenic_ligament_sim.obj",
            HERNIA_RESOURCE_DIR "/FundoTest1/gastrosplenic_ligament_sim_fixed.obj",
            1000.0,
            tissueObj);
    scene->addSceneObject(gastrosplenicLigamentObj);*/
    #pragma endregion

    // Add picking interaction for both jaws of the tool
    auto leftGrasping = std::make_shared<PbdObjectGrasping>(stomachObj);
    leftGrasping->setGeometryToPick(stomachObj->getCollidingGeometry(),
        std::dynamic_pointer_cast<PointwiseMap>(stomachObj->getPhysicsToCollidingMap()));
    leftGrasping->setStiffness(graspStiffness);
    scene->addInteraction(leftGrasping);
    auto rightGrasping = std::make_shared<PbdObjectGrasping>(stomachObj);
    rightGrasping->setGeometryToPick(stomachObj->getCollidingGeometry(),
        std::dynamic_pointer_cast<PointwiseMap>(stomachObj->getPhysicsToCollidingMap()));
    rightGrasping->setStiffness(graspStiffness);
    scene->addInteraction(rightGrasping);

    auto leftToolCollision = std::make_shared<PbdObjectCollision>(stomachObj, leftToolObj);
    leftToolCollision->setFriction(0.0);
    leftToolCollision->setRestitution(0.0);
    leftToolCollision->setRigidBodyCompliance(0.0000001);
    leftToolCollision->setUseCorrectVelocity(false);
    scene->addInteraction(leftToolCollision);
    auto rightToolCollision = std::make_shared<PbdObjectCollision>(stomachObj, rightToolObj);
    rightToolCollision->setFriction(0.0);
    rightToolCollision->setRestitution(0.0);
    rightToolCollision->setRigidBodyCompliance(0.0000001);
    rightToolCollision->setUseCorrectVelocity(false);
    scene->addInteraction(rightToolCollision);

    auto diaphgramCollision = std::make_shared<PbdObjectCollision>(stomachObj, diaphgramObj);
    diaphgramCollision->setFriction(0.0);
    diaphgramCollision->setRestitution(0.0);
    diaphgramCollision->setDeformableStiffnessA(0.3);
    scene->addInteraction(diaphgramCollision);

    // Light
    auto light = std::make_shared<DirectionalLight>();
    light->setFocalPoint(Vec3d(5.0, -8.0, -5.0));
    light->setIntensity(1.0);
    scene->addLight("Light", light);

    // Run the simulation
    {
        // Setup a viewer to render
        auto viewer = std::make_shared<VTKViewer>();
        viewer->setVtkLoggerMode(VTKViewer::VTKLoggerMode::MUTE);
        viewer->setActiveScene(scene);
        viewer->setDebugAxesLength(0.1, 0.1, 0.1);

        // Setup a scene manager to advance the scene
        auto sceneManager = std::make_shared<SceneManager>();
        sceneManager->setActiveScene(scene);
        sceneManager->pause(); // Start simulation paused

        // Setup a haptic manager to render haptics on a separate thread
        auto hapticManager = DeviceManagerFactory::makeDeviceManager();

        auto simManager = std::make_shared<SimulationManager>();
        simManager->addModule(viewer);
        simManager->addModule(sceneManager);
        simManager->addModule(hapticManager);
        simManager->setDesiredDt(executionTimestep);

        // Setup left haptic device
        std::shared_ptr<DeviceClient> leftHapticDevice = hapticManager->makeDeviceClient("LeftDevice");
        auto                          leftController   = leftToolObj->getComponent<PbdObjectController>();
        leftController->setDevice(leftHapticDevice);
        leftController->setTranslationOffset(leftPortHolePos -
            (focalPoint - leftPortHolePos).normalized() * lapToolCapsuleLength * 0.5);

#ifdef USE_TWO_HAPTIC_DEVICES
        // Setup right haptic device
        std::shared_ptr<DeviceClient> rightHapticDevice = hapticManager->makeDeviceClient("RightDevice");
        auto                          rightController   = rightToolObj->getComponent<PbdObjectController>();
        rightController->setDevice(rightHapticDevice);
        rightController->setTranslationOffset(rightPortHolePos -
            (focalPoint - rightPortHolePos).normalized() * lapToolCapsuleLength * 0.5);
#else
        auto rightHapticDevice = std::make_shared<DummyClient>();
        auto                          rightController = rightToolObj->getComponent<PbdObjectController>();
        rightController->setDevice(rightHapticDevice);
#endif

        auto camController = cameraToolObj->getComponent<PbdObjectController>();
        camController->setTranslationOffset(cameraPortHolePos -
            (focalPoint - cameraPortHolePos).normalized() * lapToolCapsuleLength * 0.5);

        auto rightPortHole = rightToolObj->addComponent<PortHoleInteraction>();
        rightPortHole->setTool(rightToolObj);
        (*rightToolObj->getPbdBody()->vertices)[0]     = Vec3d(rightPortHolePos);
        (*rightToolObj->getPbdBody()->orientations)[0] =
            Quatd::FromTwoVectors(Vec3d(0.0, 0.0, -1.0), (focalPoint - rightPortHolePos).normalized());
        rightPortHole->setPortHoleLocation(rightPortHolePos);
        rightPortHole->setToolGeometry(rightToolObj->getCollidingGeometry());
        rightPortHole->setCompliance(0.00000001);
        auto rightPortVisuals = rightToolObj->addComponent<VisualModel>();
        auto rightPortSphere = std::make_shared<Sphere>(rightPortHolePos, 0.005);
        rightPortVisuals->setGeometry(rightPortSphere);

        auto leftPortHole = leftToolObj->addComponent<PortHoleInteraction>();
        leftPortHole->setTool(leftToolObj);
        (*leftToolObj->getPbdBody()->vertices)[0]     = Vec3d(leftPortHolePos);
        (*leftToolObj->getPbdBody()->orientations)[0] =
            Quatd::FromTwoVectors(Vec3d(0.0, 0.0, -1.0), (focalPoint - leftPortHolePos).normalized());
        leftPortHole->setPortHoleLocation(leftPortHolePos);
        leftPortHole->setToolGeometry(leftToolObj->getCollidingGeometry());
        leftPortHole->setCompliance(0.00000001);
        auto leftPortVisuals = leftToolObj->addComponent<VisualModel>();
        auto leftPortSphere = std::make_shared<Sphere>(leftPortHolePos, 0.005);
        leftPortVisuals->setGeometry(leftPortSphere);

        auto cameraPortHole = cameraToolObj->addComponent<PortHoleInteraction>();
        cameraPortHole->setTool(cameraToolObj);
        (*cameraToolObj->getPbdBody()->vertices)[0]     = Vec3d(cameraPortHolePos);
        (*cameraToolObj->getPbdBody()->orientations)[0] =
            Quatd::FromTwoVectors(Vec3d(0.0, 0.0, -1.0), (focalPoint - cameraPortHolePos).normalized());
        cameraPortHole->setPortHoleLocation(cameraPortHolePos);
        cameraPortHole->setToolGeometry(cameraToolObj->getCollidingGeometry());
        cameraPortHole->setCompliance(0.00000001);
        auto camPortVisuals = cameraToolObj->addComponent<VisualModel>();
        auto camPortSphere = std::make_shared<Sphere>(cameraPortHolePos, 0.005);
        camPortVisuals->setGeometry(camPortSphere);

        // Add mouse and keyboard controls to the viewer
        {
            // Couple extra key press controls
            connect<KeyEvent>(viewer->getKeyboardDevice(), &KeyboardDeviceClient::keyPress, [&](KeyEvent* e)
                {
                    if (e->m_key == 'u')
                    {
                        LOG(INFO) << "One sim step + render";
                        scene->advance(sceneManager->getDt());
                        viewer->update();
                    }
                    else if (e->m_key == 'f')
                    {
                        std::shared_ptr<DeviceClient> tmpDevice = leftController->getDevice();
                        leftController->setDevice(camController->getDevice());
                        camController->setDevice(tmpDevice);
                    }
                });
            connect<Event>(viewer, &VTKViewer::preUpdate, [&](Event*)
                {
                    auto capsule        = std::dynamic_pointer_cast<Capsule>(cameraToolObj->getCollidingGeometry());
                    const Vec3d dir     = -capsule->getRotation().col(2).normalized();
                    const Vec3d up      = capsule->getRotation().col(1).normalized();
                    const double r      = capsule->getRadius();
                    const double length = capsule->getLength();

                    std::shared_ptr<Camera> cam = scene->getCamera("default");
                    cam->setPosition(capsule->getPosition() + dir * (length * 0.48));
                    cam->setFocalPoint(capsule->getPosition() + dir * (length * 0.5));
                    cam->setViewUp(up);
                });
            connect<ButtonEvent>(leftHapticDevice, &DeviceClient::buttonStateChanged,
                [&](ButtonEvent* e)
                {
                    if (e->m_buttonState == BUTTON_PRESSED)
                    {
                        if (e->m_button == 0)
                        {
                            LOG(INFO) << "Left Grasp!";
                            //leftGrasping->beginVertexGrasp(std::dynamic_pointer_cast<AnalyticalGeometry>(leftToolObj->getCollidingGeometry()));
                            //leftToolCollision->setEnabled(false);
                            auto graspGeom = std::dynamic_pointer_cast<AnalyticalGeometry>(leftToolObj->getVisualModel(1)->getGeometry());
                            leftGrasping->beginCellGrasp(graspGeom);
                        }
                    }
                    else if (e->m_buttonState == BUTTON_RELEASED)
                    {
                        if (e->m_button == 0)
                        {
                            LOG(INFO) << "Left Release!";
                            leftGrasping->endGrasp();
                            //leftToolCollision->setEnabled(true);
                        }
                    }
                });
            connect<ButtonEvent>(rightHapticDevice, &DeviceClient::buttonStateChanged,
                [&](ButtonEvent* e)
                {
                    if (e->m_buttonState == BUTTON_PRESSED)
                    {
                        if (e->m_button == 0)
                        {
                            LOG(INFO) << "Right Grasp!";
                            auto graspGeom = std::dynamic_pointer_cast<AnalyticalGeometry>(rightToolObj->getVisualModel(1)->getGeometry());
                            rightGrasping->beginCellGrasp(graspGeom);
                        }
                    }
                    else if (e->m_buttonState == BUTTON_RELEASED)
                    {
                        if (e->m_button == 0)
                        {
                            LOG(INFO) << "Right Release!";
                            rightGrasping->endGrasp();
                            //rightToolCollision->setEnabled(true);
                        }
                    }
                });

            // Add default mouse and keyboard controls to the viewer
            std::shared_ptr<Entity> mouseAndKeyControls =
                SimulationUtils::createDefaultSceneControl(simManager);
            scene->addSceneObject(mouseAndKeyControls);
        }

        simManager->start();
    }
}