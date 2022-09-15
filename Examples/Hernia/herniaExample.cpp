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

using namespace imstk;

double tissueCapsuleCompliance = 0.0;

double capsuleMass = 1.0;
double capsuleInertiaScale = 0.01;
double capsuleDistConstraintStiffness      = 100.0;
double capsuleHingeConstraintCompliance    = 0.000001;
double capsuleHingeConstraintEndCompliance = 0.001; // The two end points of the capsule chain

double tissueParticleMass      = 0.1;
double tissueDistStiffness     = 10000.0;
double tissueDihedralStiffness = 1.0;

double globalLinearDamping  = 0.075;
double globalAngularDamping = 0.01;

int    iterations = 1; // Prefer lower timestep with XPBD
int    collisionIterations = 1;
double simTimestep       = 0.0025;
double executionTimestep = 0.0025;

double graspStiffness = 0.01;

double tetYoungsModulus = 4000.0;
double tetPossionsRatio = 0.47;

const Vec3d rightPortHolePos  = Vec3d(0.015, 0.092, -1.117);
const Vec3d leftPortHolePos   = Vec3d(-0.065, 0.078, -1.127);
const Vec3d cameraPortHolePos = Vec3d(-0.004, 0.093, -1.2);

// There a couple major things happening in this example
// - We have a deformable stomach
// - We a have a diaphgram SDF which helps hold the stomach in place
// - We have 3 ligaments attached to the diaphgram which help hold it in place.
//  In particular the gastrosplenic ligament helps hold the stomach back when it
//  is pulled to the left.
// - We have 2 tools, though 3 is desirable, this allows wrapping of the stomach around
//  the esophagus.
// - We have stuck 3 capsules down the esophagus and 3 along the edge of the lesser omentum.
//  this avoids mesh-mesh collision whilst keeping the shape of the esohpagus, it is the most
//  critical component that allows the wrapping.

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
    // Setup the Geometry
    /*std::shared_ptr<SurfaceMesh> tissueMesh =
        MeshIO::read<SurfaceMesh>("C:/Users/Andx_/Desktop/NewHernia/stomach_surface.obj");*/

    /* auto testSurfMesh = MeshIO::read<SurfaceMesh>(VLAHHS_DATA_PATH"/volume/stomach_surface_tmp.obj");
     Tetrahedralize tetrahedralize;
     tetrahedralize.setInput(testSurfMesh);
     tetrahedralize.update();
     std::shared_ptr<TetrahedralMesh> tetMesh = tetrahedralize.getOutputTetMesh();*/

    std::shared_ptr<TetrahedralMesh> tissueMesh =
        MeshIO::read<TetrahedralMesh>(VLAHHS_DATA_PATH "/volume/stomach_surface_tmp_.msh");
    std::shared_ptr<TetrahedralMesh> esophagusMesh =
        MeshIO::read<TetrahedralMesh>(VLAHHS_DATA_PATH "/volume/esophagus_tets.vtk");
    std::shared_ptr<SurfaceMesh> tissueSurfMesh = tissueMesh->extractSurfaceMesh();
    std::shared_ptr<SurfaceMesh> colMesh = MeshIO::read<SurfaceMesh>(
        VLAHHS_DATA_PATH "/volume/stomach_surface_collision.obj");

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
    auto visualMap = std::make_shared<PointwiseMap>(tissueMesh, tissueSurfMesh);
    visualMap->setTolerance(0.0001);
    tissueObj->setPhysicsToVisualMap(visualMap);
    tissueObj->setCollidingGeometry(colMesh);
    auto colMap = std::make_shared<PointwiseMap>(tissueMesh, colMesh);
    colMap->setTolerance(0.0001);
    tissueObj->setPhysicsToCollidingMap(colMap);
    tissueObj->getVisualModel(0)->setRenderMaterial(material);
    tissueObj->setPhysicsGeometry(tissueMesh);

    tissueObj->setDynamicalModel(model);
    tissueObj->getPbdBody()->uniformMassValue = tissueParticleMass;

    auto         fixedVerts = MeshIO::read<PointSet>(VLAHHS_DATA_PATH "/volume/fixed_verts.obj");
    PointwiseMap fixedMapper;
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
        MeshIO::read<LineMesh>(VLAHHS_DATA_PATH "/capsules.obj");
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
        //capsuleObj->setVisualGeometry(capsule);
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
            // For every capsule
            //for (int i = 0; i < maIndices.size(); i++)
            //{
            //    auto capsule = std::dynamic_pointer_cast<Capsule>(capsuleObjs[i]->getVisualGeometry());
            //    const Vec3d center = (*capsuleObjs[i]->getPbdBody()->vertices)[0];
            //    const Quatd orientation = (*capsuleObjs[i]->getPbdBody()->orientations)[0];
            //    const Vec3d dir = orientation._transformVector(Vec3d(0.0, 1.0, 0.0));
            //    const PbdParticleId fixedParticle =
            //        model->addVirtualParticle(
            //            center, // Position
            //            orientation, // Orientation
            //            0.0, // Mass
            //            Mat3d::Identity(), // Inertia
            //            Vec3d::Zero(), // Velocity
            //            Vec3d::Zero(), // Angular Velocity
            //            true); // Persist

            //    auto distConstraint = std::make_shared<PbdDistanceConstraint>();
            //    distConstraint->initConstraint(
            //        0.0,
            //        { capsuleObjs[i]->getPbdBody()->bodyHandle, 0 },
            //        fixedParticle, capsuleDistConstraintStiffness);
            //    container.addConstraint(distConstraint);

            //    /*auto angularDistConstraint = std::make_shared<PbdAngularDistanceConstraint>();
            //    angularDistConstraint->initConstraint(
            //        { capsuleObjs[i]->getPbdBody()->bodyHandle, 0 },
            //        fixedParticle,
            //        capsuleHingeConstraintCompliance);
            //    container.addConstraint(angularDistConstraint);*/
            //}

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

void
makeCollidingCapsules(const std::string&                             name,
                      std::vector<std::shared_ptr<CollidingObject>>& capsuleObjs,
                      std::shared_ptr<PbdModel>                      model,
                      const double                                   capsuleRadius)
{
    capsuleObjs.clear();

    // Read in a medial line
    auto maLineMesh =
        MeshIO::read<LineMesh>(VLAHHS_DATA_PATH "/capsules.obj");
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

        const Quatd orientation = Quatd::FromTwoVectors(Vec3d(0.0, 1.0, 0.0), dir);
        auto        capsule     = std::make_shared<Capsule>(center, capsuleRadius, diff.norm(), orientation);
        //auto capsule = std::make_shared<Sphere>(Vec3d::Zero(), capsuleRadius);

        auto capsuleObj = std::make_shared<CollidingObject>(name);
        capsuleObj->setVisualGeometry(capsule);
        capsuleObj->setCollidingGeometry(capsule);
        //capsuleObj->getVisualModel(0)->getRenderMaterial()->setColor(Color::Red);
        //capsuleObj->getVisualModel(0)->getRenderMaterial()->setOpacity(0.5);

        capsuleObjs.push_back(capsuleObj);
    }
}

std::shared_ptr<PbdObject>
makeLapToolObj(const std::string&        name,
               std::shared_ptr<PbdModel> model)
{
    auto lapTool = std::make_shared<PbdObject>(name);

    const double capsuleLength = 0.3;
    auto         toolGeom      = std::make_shared<Capsule>(Vec3d(0.0, 0.0, capsuleLength * 0.5),
        0.002, capsuleLength,
        Quatd::FromTwoVectors(Vec3d(0.0, 1.0, 0.0), Vec3d(0.0, 0.0, 1.0)));

    // Slightly large one for grasping
    auto toolGeomDilated = std::make_shared<Capsule>(Vec3d(0.0, 0.0, capsuleLength * 0.5),
        0.0022, capsuleLength,
        Quatd::FromTwoVectors(Vec3d(0.0, 1.0, 0.0), Vec3d(0.0, 0.0, 1.0)));

    lapTool->setDynamicalModel(model);
    lapTool->setPhysicsGeometry(toolGeom);
    lapTool->setCollidingGeometry(toolGeom);
    lapTool->setVisualGeometry(toolGeomDilated);
    lapTool->setPhysicsToVisualMap(std::make_shared<IsometricMap>(toolGeom, toolGeomDilated));

    std::shared_ptr<RenderMaterial> material = lapTool->getVisualModel(0)->getRenderMaterial();
    material->setIsDynamicMesh(false);
    material->setMetalness(1.0);
    material->setRoughness(0.2);
    material->setShadingModel(RenderMaterial::ShadingModel::PBR);

    lapTool->getPbdBody()->setRigid(
        Vec3d(0.0, 0.0, 0.0),
        10.0,
        Quatd::Identity(),
        Mat3d::Identity() * 0.1);

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
        MeshIO::read<ImageData>(VLAHHS_DATA_PATH "/stomach_container_SDF.mhd")
        ->cast(IMSTK_DOUBLE);
    diaphgramObj->setCollidingGeometry(std::make_shared<SignedDistanceField>(diaphgramSdfImage));

    /*imstkNew<SurfaceMeshFlyingEdges> isoExtract;
    isoExtract->setInputImage(diaphgramSdfImage);
    isoExtract->update();
    diaphgramObj->setVisualGeometry(isoExtract->getOutputMesh());
    diaphgramObj->getVisualModel(0)->getRenderMaterial()->setOpacity(0.5);*/

    auto surfMesh = MeshIO::read<SurfaceMesh>(VLAHHS_DATA_PATH "/diaphragm.obj");
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
    scene->getActiveCamera()->setFocalPoint(0.0158845, -0.00876985, -1.21677);
    scene->getActiveCamera()->setPosition(0.0784281, 0.333918, -1.33171);
    scene->getActiveCamera()->setViewUp(-0.00193505, -0.32117, -0.947019);
    scene->getConfig()->debugCamBoundingBox = false;
    *scene->getCamera("debug") = *scene->getActiveCamera();

    auto bodyObject = std::make_shared<SceneObject>("body");
    {
        auto surfMesh = MeshIO::read<SurfaceMesh>(iMSTK_DATA_ROOT "/human/full_body/body.obj");
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
    pbdConfig->m_gravity = Vec3d(0.0, 0.0, 0.0);
    //pbdConfig->m_gravity = Vec3d(0.0, -20.0, 0.0);
    //pbdConfig->m_gravity = Vec3d(9.8, 0.0, 0.0);
    pbdConfig->m_dt = simTimestep;
    pbdConfig->m_iterations = iterations;
    pbdConfig->m_linearDampingCoeff  = globalLinearDamping;
    pbdConfig->m_angularDampingCoeff = globalAngularDamping;
    pbdConfig->m_doPartitioning      = false;
    //pbdConfig->m_collisionSolverType = PbdConstraint::SolverType::PBD;

    std::shared_ptr<PbdObject> tissueObj = makeTissueObj("Tissue", pbdModel);
    scene->addSceneObject(tissueObj);

    std::shared_ptr<PbdObject> leftToolObj = makeLapToolObj("leftToolObj", pbdModel);
    scene->addSceneObject(leftToolObj);
    std::shared_ptr<PbdObject> rightToolObj = makeLapToolObj("rightToolObj", pbdModel);
    scene->addSceneObject(rightToolObj);
    std::shared_ptr<PbdObject> cameraToolObj = makeLapToolObj("cameraToolObj", pbdModel);
    (*cameraToolObj->getPbdBody()->vertices)[0]     = Vec3d(0.001, 0.078, -1.206);
    (*cameraToolObj->getPbdBody()->orientations)[0] = Quatd(0.826, -0.54, -0.149, 0.061);
    scene->addSceneObject(cameraToolObj);

    std::vector<std::shared_ptr<PbdObject>> esophagusCapsuleObjs;
    makeCapsules("CapsuleObjs", esophagusCapsuleObjs, pbdModel, 0.004);
    for (int i = 0; i < esophagusCapsuleObjs.size(); i++)
    {
        scene->addSceneObject(esophagusCapsuleObjs[i]);
        auto capsCollision = std::make_shared<PbdObjectCollision>(
            tissueObj, esophagusCapsuleObjs[i]);
        capsCollision->setRestitution(0.0);
        capsCollision->setFriction(0.0);
        capsCollision->setUseCorrectVelocity(false);
        capsCollision->setRigidBodyCompliance(0.00001);
        scene->addInteraction(capsCollision);
    }
    //std::vector<std::shared_ptr<CollidingObject>> esophagusCapsuleObjs;
    //makeCollidingCapsules("CapsuleObjs", esophagusCapsuleObjs, pbdModel, 0.003);
    //for (int i = 0; i < esophagusCapsuleObjs.size(); i++)
    //{
    //    scene->addSceneObject(esophagusCapsuleObjs[i]);
    //    auto capsCollision = std::make_shared<PbdObjectCollision>(
    //        tissueObj, esophagusCapsuleObjs[i]);
    //    capsCollision->setRestitution(0.0);
    //    capsCollision->setFriction(0.0);
    //    //capsCollision->setRigidBodyCompliance(tissueCapsuleCompliance);
    //    capsCollision->setDeformableStiffnessA(0.1);
    //    scene->addInteraction(capsCollision);
    //}

    /*std::vector<std::shared_ptr<CollidingObject>> fatEdgeCapsuleObjs;
    makeFatCapsules("fatCapsules", fatEdgeCapsuleObjs, pbdModel, 0.0025);
    for (int i = 0; i < fatEdgeCapsuleObjs.size(); i++)
    {
        scene->addSceneObject(fatEdgeCapsuleObjs[i]);
        auto capsCollision = std::make_shared<PbdObjectCollision>(
            tissueObj, fatEdgeCapsuleObjs[i]);
        capsCollision->setRestitution(0.0);
        capsCollision->setRigidBodyCompliance(tissueCapsuleCompliance);
        scene->addInteraction(capsCollision);
    }*/

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
    auto leftGrasping = std::make_shared<PbdObjectGrasping>(tissueObj);
    leftGrasping->setGeometryToPick(tissueObj->getVisualGeometry(),
        std::dynamic_pointer_cast<PointwiseMap>(tissueObj->getPhysicsToVisualMap()));
    leftGrasping->setStiffness(graspStiffness);
    scene->addInteraction(leftGrasping);
    auto rightGrasping = std::make_shared<PbdObjectGrasping>(tissueObj);
    rightGrasping->setGeometryToPick(tissueObj->getVisualGeometry(),
        std::dynamic_pointer_cast<PointwiseMap>(tissueObj->getPhysicsToVisualMap()));
    rightGrasping->setStiffness(graspStiffness);
    scene->addInteraction(rightGrasping);

    auto leftToolCollision = std::make_shared<PbdObjectCollision>(tissueObj, leftToolObj);
    leftToolCollision->setFriction(0.0);
    leftToolCollision->setRestitution(0.0);
    leftToolCollision->setRigidBodyCompliance(0.0000001);
    leftToolCollision->setUseCorrectVelocity(false);
    scene->addInteraction(leftToolCollision);
    /* auto leftToolFatCollision = std::make_shared<PbdObjectCollision>(fatObj, leftToolObj);
     leftToolFatCollision->setFriction(0.0);
     leftToolFatCollision->setRestitution(0.0);
     leftToolFatCollision->setRigidBodyCompliance(0.0000001);
     scene->addInteraction(leftToolFatCollision);*/
    auto rightToolCollision = std::make_shared<PbdObjectCollision>(tissueObj, rightToolObj);
    rightToolCollision->setFriction(0.0);
    rightToolCollision->setRestitution(0.0);
    rightToolCollision->setRigidBodyCompliance(0.0000001);
    rightToolCollision->setUseCorrectVelocity(false);
    scene->addInteraction(rightToolCollision);

    auto diaphgramCollision = std::make_shared<PbdObjectCollision>(tissueObj, diaphgramObj);
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

        auto ohHapticManager = DeviceManagerFactory::makeDeviceManager("OpenHapticDeviceManager");

        auto driver = std::make_shared<SimulationManager>();
        driver->addModule(viewer);
        driver->addModule(sceneManager);
        driver->addModule(ohHapticManager);
        driver->setDesiredDt(executionTimestep);

        // Setup left haptic device
        std::shared_ptr<DeviceClient> leftHapticDevice = ohHapticManager->makeDeviceClient("LeftDevice");
        auto                          leftController   = leftToolObj->getComponent<PbdObjectController>();
        leftController->setDevice(leftHapticDevice);
        leftController->setTranslationOffset(Vec3d(0.0, 0.0, -1.2));

        // Setup right haptic device
        std::shared_ptr<DeviceClient> rightHapticDevice = ohHapticManager->makeDeviceClient("RightDevice");
        auto                          rightController   = rightToolObj->getComponent<PbdObjectController>();
        rightController->setDevice(rightHapticDevice);
        rightController->setTranslationOffset(Vec3d(0.1, 0.0, -1.2));

        auto camController = cameraToolObj->getComponent<PbdObjectController>();
        camController->setTranslationOffset(Vec3d(0.0, 0.0, -1.2));

        (*rightToolObj->getPbdBody()->vertices)[0] = Vec3d(rightPortHolePos);
        auto rightPortHole = rightToolObj->addComponent<PortHoleInteraction>();
        rightPortHole->setTool(rightToolObj);
        rightPortHole->setPortHoleLocation(rightPortHolePos);
        rightPortHole->setToolGeometry(rightToolObj->getCollidingGeometry());
        rightPortHole->setCompliance(0.00000001);
        auto rightPortVisuals = rightToolObj->addComponent<VisualModel>();
        auto rightPortSphere  = std::make_shared<Sphere>(rightPortHolePos, 0.005);
        rightPortVisuals->setGeometry(rightPortSphere);

        (*leftToolObj->getPbdBody()->vertices)[0] = Vec3d(leftPortHolePos);
        auto leftPortHole = leftToolObj->addComponent<PortHoleInteraction>();
        leftPortHole->setTool(leftToolObj);
        leftPortHole->setPortHoleLocation(leftPortHolePos);
        leftPortHole->setToolGeometry(leftToolObj->getCollidingGeometry());
        leftPortHole->setCompliance(0.00000001);
        auto leftPortVisuals = leftToolObj->addComponent<VisualModel>();
        auto leftPortSphere  = std::make_shared<Sphere>(leftPortHolePos, 0.005);
        leftPortVisuals->setGeometry(leftPortSphere);

        //(*cameraToolObj->getPbdBody()->vertices)[0] = Vec3d(cameraPortHolePos);
        auto cameraPortHole = cameraToolObj->addComponent<PortHoleInteraction>();
        cameraPortHole->setTool(cameraToolObj);
        cameraPortHole->setPortHoleLocation(cameraPortHolePos);
        cameraPortHole->setToolGeometry(cameraToolObj->getCollidingGeometry());
        cameraPortHole->setCompliance(0.00000001);
        auto camPortVisuals = cameraToolObj->addComponent<VisualModel>();
        auto camPortSphere  = std::make_shared<Sphere>(cameraPortHolePos, 0.005);
        camPortVisuals->setGeometry(camPortSphere);

        // Add mouse and keyboard controls to the viewer
        {
            // Couple extra key press controls
            connect<KeyEvent>(viewer->getKeyboardDevice(), &KeyboardDeviceClient::keyPress, [&](KeyEvent* e)
                {
                    if (e->m_key == 'g')
                    {
                        if (pbdModel->getConfig()->m_gravity[0] == 0.0)
                        {
                            pbdModel->getConfig()->m_gravity = Vec3d(0.0, -10.0, 0.0);
                        }
                        else
                        {
                            pbdModel->getConfig()->m_gravity = Vec3d::Zero();
                        }
                    }
                    else if (e->m_key == 'u')
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

                    scene->getCamera("default")->setPosition(capsule->getPosition() + dir * (length * 0.48));
                    scene->getCamera("default")->setFocalPoint(capsule->getPosition() + dir * (length * 0.5));
                    scene->getCamera("default")->setViewUp(up);
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
                            leftGrasping->beginCellGrasp(
                                std::dynamic_pointer_cast<AnalyticalGeometry>(leftToolObj->getVisualGeometry()));
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
                            //rightGrasping->beginVertexGrasp(std::dynamic_pointer_cast<AnalyticalGeometry>(leftToolObj->getCollidingGeometry()));
                            //rightToolCollision->setEnabled(false);
                            rightGrasping->beginCellGrasp(
                                std::dynamic_pointer_cast<AnalyticalGeometry>(rightToolObj->getVisualGeometry()));
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

            // Update Fat capsules
            /*connect<Event>(sceneManager, &SceneManager::postUpdate, [&](Event*) {
                auto lineMesh = std::dynamic_pointer_cast<LineMesh>(fatObj->getCollidingGeometry());

                std::shared_ptr<VecDataArray<double, 3>> verticesPtr = lineMesh->getVertexPositions();
                const VecDataArray<double, 3>& vertices = *verticesPtr;
                std::shared_ptr<VecDataArray<int, 2>> indicesPtr = lineMesh->getCells();
                const VecDataArray<int, 2>& indices = *indicesPtr;

                for (int i = 0; i < indices.size(); i++)
                {
                    const Vec2i& cell = indices[i];
                    const Vec3d p = vertices[cell[0]];
                    const Vec3d q = vertices[cell[1]];
                    const Vec3d center = (q + p) * 0.5;
                    const Vec3d diff = q - p;
                    const Vec3d dir = diff.normalized();

                    auto capsuleObj = std::dynamic_pointer_cast<PbdObject>(fatEdgeCapsuleObjs[i]);
                    (*capsuleObj->getPbdBody()->vertices)[0] = center;
                    (*capsuleObj->getPbdBody()->orientations)[0] = Quatd::FromTwoVectors(Vec3d(0.0, 1.0, 0.0), dir);

                    auto capsule =
                        std::dynamic_pointer_cast<Capsule>(fatEdgeCapsuleObjs[i]->getCollidingGeometry());
                    capsule->setLength(std::max(diff.norm(), 0.0));
                }
                });*/

            // Add default mouse and keyboard controls to the viewer
            std::shared_ptr<Entity> mouseAndKeyControls =
                SimulationUtils::createDefaultSceneControl(driver);
            scene->addSceneObject(mouseAndKeyControls);
        }

        driver->start();
    }
}