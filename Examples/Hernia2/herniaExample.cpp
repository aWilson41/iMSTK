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
#include "imstkPlane.h"
#include "imstkPointwiseMap.h"
#include "imstkPortHoleInteraction.h"
#include "imstkRenderMaterial.h"
#include "imstkScene.h"
#include "imstkSceneManager.h"
#include "imstkSimulationManager.h"
#include "imstkSimulationUtils.h"
#include "imstkSphere.h"
#include "imstkVisualModel.h"
#include "imstkVTKViewer.h"

using namespace imstk;

const double tissueCapsuleCompliance = 0.0;

const double capsuleMass = 1.0;
const double capsuleInertiaScale = 0.01;
const double capsuleDistConstraintStiffness      = 100.0;
const double capsuleHingeConstraintCompliance    = 0.000001;
const double capsuleHingeConstraintEndCompliance = 0.1;

const double tissueParticleMass      = 1.0;
const double tissueDistStiffness     = 1000.0;
const double tissueDihedralStiffness = 1.0;

const double globalLinearDamping  = 0.05;
const double globalAngularDamping = 0.05;

const int    iterations        = 1;
const double simTimestep       = 0.0025;
const double executionTimestep = 0.004;

const double graspStiffness = 1.0;

const Vec3d rightPortHolePos = Vec3d(0.015, 0.092, -1.117);
const Vec3d leftPortHolePos  = Vec3d(-0.065, 0.078, -1.127);

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

void
makeCapsules(const std::string&                       name,
             std::vector<std::shared_ptr<PbdObject>>& capsuleObjs,
             std::shared_ptr<PbdModel>                model,
             const double                             capsuleRadius)
{
    capsuleObjs.clear();

    // Read in a medial line
    auto maLineMesh =
        MeshIO::read<LineMesh>(VLAHHS_DATA_PATH "/volume2/capsules.obj");
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
                auto capsule0 = std::dynamic_pointer_cast<Capsule>(capsuleObjs[i]->getVisualGeometry());
                Vec3d jp0     = Vec3d::Zero();
                {
                    const Vec3d& capsulePos         = capsule0->getPosition();
                    const double capsuleLength      = capsule0->getLength();
                    const Quatd& capsuleOrientation = capsule0->getOrientation();
                    const Vec3d& capsulePosA        = capsulePos - 0.5 * capsuleLength * capsuleOrientation.toRotationMatrix().col(1);
                    const Vec3d& capsulePosB        = capsulePos + (capsulePos - capsulePosA);
                    jp0 = capsulePosB;
                }
                auto capsule1 = std::dynamic_pointer_cast<Capsule>(capsuleObjs[i + 1]->getVisualGeometry());
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
            //{
            //    // Assumes consistent directions of the capsules
            //    auto capsule0 = std::dynamic_pointer_cast<Capsule>(capsuleObjs[0]->getVisualGeometry());
            //    Vec3d jp0 = Vec3d::Zero();
            //    {
            //        const Vec3d& capsulePos = capsule0->getPosition();
            //        const double capsuleLength = capsule0->getLength();
            //        const Quatd& capsuleOrientation = capsule0->getOrientation();
            //        const Vec3d& capsulePosA = capsulePos - 0.5 * capsuleLength * capsuleOrientation.toRotationMatrix().col(1);
            //        const Vec3d& capsulePosB = capsulePos + (capsulePos - capsulePosA);
            //        jp0 = capsulePosA;
            //    }
            //    // Assumes consistent directions of the capsules
            //    auto capsule1 = std::dynamic_pointer_cast<Capsule>(capsuleObjs.back()->getVisualGeometry());
            //    Vec3d jp1 = Vec3d::Zero();
            //    {
            //        const Vec3d& capsulePos = capsule1->getPosition();
            //        const double capsuleLength = capsule1->getLength();
            //        const Quatd& capsuleOrientation = capsule1->getOrientation();
            //        const Vec3d& capsulePosA = capsulePos - 0.5 * capsuleLength * capsuleOrientation.toRotationMatrix().col(1);
            //        const Vec3d& capsulePosB = capsulePos + (capsulePos - capsulePosA);
            //        jp1 = capsulePosB;
            //    }

            //    const PbdParticleId& fixedVid0 = model->addVirtualParticle(jp0, 0.0, Vec3d::Zero(), true);
            //    auto bodyToBodyConstraint0 = std::make_shared<PbdBodyToBodyDistanceConstraint>();
            //    bodyToBodyConstraint0->initConstraint(model->getBodies(),
            //        { capsuleObjs[0]->getPbdBody()->bodyHandle, 0 },
            //        jp0,
            //        fixedVid0,
            //        0.0, capsuleHingeConstraintEndCompliance);
            //    container.addConstraint(bodyToBodyConstraint0);

            //    const PbdParticleId& fixedVid1 = model->addVirtualParticle(jp1, 0.0, Vec3d::Zero(), true);
            //    auto bodyToBodyConstraint1 = std::make_shared<PbdBodyToBodyDistanceConstraint>();
            //    bodyToBodyConstraint1->initConstraint(model->getBodies(),
            //        { capsuleObjs.back()->getPbdBody()->bodyHandle, 0 },
            //        jp1,
            //        fixedVid1,
            //        0.0, capsuleHingeConstraintEndCompliance);
            //    container.addConstraint(bodyToBodyConstraint1);
            //}
            {
            }
        });
}

///
/// \brief Creates tissue object
///
static std::shared_ptr<PbdObject>
makeTissueObj(const std::string&        name,
              std::shared_ptr<PbdModel> model)
{
    // Setup the Geometry
    std::shared_ptr<TetrahedralMesh> tissueMesh =
        MeshIO::read<TetrahedralMesh>(VLAHHS_DATA_PATH "/volume2/stomach_surface_.msh");
    std::shared_ptr<TetrahedralMesh> esophagusMesh =
        MeshIO::read<TetrahedralMesh>(VLAHHS_DATA_PATH "/volume2/esophagus.vtk");
    std::shared_ptr<SurfaceMesh> tissueSurfMesh = tissueMesh->extractSurfaceMesh();

    // Compute coincident tets to generate a mask 0 if part of stomach, 1 if part of esophagus
    {
        int  foundCount = 0;
        auto maskPtr    = std::make_shared<DataArray<unsigned char>>(tissueMesh->getNumCells());
        maskPtr->fill(0);
        tissueMesh->setVertexAttribute("mask", maskPtr);
        for (int i = 0; i < esophagusMesh->getNumCells(); i++)
        {
            const Vec4i& esoTet  = (*esophagusMesh->getCells())[i];
            const Vec3d& a0      = (*esophagusMesh->getVertexPositions())[esoTet[0]];
            const Vec3d& a1      = (*esophagusMesh->getVertexPositions())[esoTet[1]];
            const Vec3d& a2      = (*esophagusMesh->getVertexPositions())[esoTet[2]];
            const Vec3d& a3      = (*esophagusMesh->getVertexPositions())[esoTet[3]];
            const Vec3d  centerA = (a0 + a1 + a2 + a3) * 0.25;

            // For every tet of the esophagus find corresponding tet in tissue
            for (int j = 0; j < tissueMesh->getNumCells(); j++)
            {
                const Vec4i& tissueTet = (*tissueMesh->getCells())[j];
                const Vec3d& b0      = (*tissueMesh->getVertexPositions())[tissueTet[0]];
                const Vec3d& b1      = (*tissueMesh->getVertexPositions())[tissueTet[1]];
                const Vec3d& b2      = (*tissueMesh->getVertexPositions())[tissueTet[2]];
                const Vec3d& b3      = (*tissueMesh->getVertexPositions())[tissueTet[3]];
                const Vec3d  centerB = (b0 + b1 + b2 + b3) * 0.25;

                // Compare all arrangements of vertices. Can we assume widing is the same?
                if (centerA.isApprox(centerB, 0.00001))
                {
                    (*maskPtr)[j] = 1;
                    printf("found tet %d, found count %d!\n", i, foundCount);
                    foundCount++;
                }
            }
        }
    }

    // Setup the material
    auto material = std::make_shared<RenderMaterial>();
    material->setColor(Color(93.0 / 255.0, 38.0 / 255.0, 37.0 / 255.0));
    material->setDisplayMode(RenderMaterial::DisplayMode::WireframeSurface);
    material->setBackFaceCulling(false);
    material->setOpacity(0.5);

    // Setup the Object
    auto tissueObj = std::make_shared<PbdObject>(name);
    tissueObj->setVisualGeometry(tissueMesh);
    tissueObj->setCollidingGeometry(tissueMesh);
    auto colMap = std::make_shared<PointwiseMap>(tissueMesh, tissueSurfMesh);
    colMap->setTolerance(0.0001);
    tissueObj->setPhysicsToCollidingMap(colMap);
    tissueObj->setCollidingGeometry(tissueSurfMesh);
    tissueObj->getVisualModel(0)->setRenderMaterial(material);
    tissueObj->setPhysicsGeometry(tissueMesh);

    tissueObj->setDynamicalModel(model);
    tissueObj->getPbdBody()->uniformMassValue = tissueParticleMass;

    // Setup the Parameters
    /*model->getConfig()->enableConstraint(PbdModelConfig::ConstraintGenType::Distance, tissueDistStiffness,
        tissueObj->getPbdBody()->bodyHandle);
    model->getConfig()->enableConstraint(PbdModelConfig::ConstraintGenType::Dihedral, tissueDihedralStiffness,
        tissueObj->getPbdBody()->bodyHandle);*/

    /*model->getConfig()->enableConstraint(PbdModelConfig::ConstraintGenType::Distance, 1000.0,
         tissueObj->getPbdBody()->bodyHandle);
     model->getConfig()->enableConstraint(PbdModelConfig::ConstraintGenType::Volume, 100.0,
         tissueObj->getPbdBody()->bodyHandle);*/

    /* model->getConfig()->enableFemConstraint(PbdFemConstraint::MaterialType::NeoHookean,
         tissueObj->getPbdBody()->bodyHandle);
     model->getConfig()->m_femParams->m_PoissonRatio = tetPossionsRatio;
     model->getConfig()->m_femParams->m_YoungModulus = tetYoungsModulus;*/

    auto customFunctor = std::make_shared<PbdFemTetConstraintFunctorCustom>();
    customFunctor->setBodyIndex(tissueObj->getPbdBody()->bodyHandle);
    const double possionsRatio = 0.4;
    const double stomachYoungsModulus   = 3000.0;
    const double esophagusYoungsModulus = 200000.0;
    auto         femConfigStomach       = std::make_shared<PbdFemConstraintConfig>(
        stomachYoungsModulus / 2.0 / (1.0 + possionsRatio), stomachYoungsModulus * possionsRatio / ((1.0 + possionsRatio) * (1.0 - 2.0 * possionsRatio)),
        stomachYoungsModulus, possionsRatio);
    customFunctor->setFemConfigStomach(femConfigStomach);
    auto femConfigEsophagus = std::make_shared<PbdFemConstraintConfig>(
        esophagusYoungsModulus / 2.0 / (1.0 + possionsRatio), esophagusYoungsModulus * possionsRatio / ((1.0 + possionsRatio) * (1.0 - 2.0 * possionsRatio)),
        esophagusYoungsModulus, possionsRatio);
    customFunctor->setFemConfigEsophagus(femConfigEsophagus);
    customFunctor->setMaterialTypeEsophagus(PbdFemConstraint::MaterialType::StVK);
    customFunctor->setMaterialTypeStomach(PbdFemConstraint::MaterialType::NeoHookean);
    model->getConfig()->addPbdConstraintFunctor(customFunctor);

    return tissueObj;
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
        0.003, capsuleLength,
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

    auto controllerGhost = lapTool->addComponent<ObjectControllerGhost>();
    controllerGhost->setController(controller);

    return lapTool;
}

int
main()
{
    // Write log to stdout and file
    Logger::startLogger();

    // Setup a scene
    auto scene = std::make_shared<Scene>("Hernia2");
    scene->getActiveCamera()->setFocalPoint(0.0158845, -0.00876985, -1.21677);
    scene->getActiveCamera()->setPosition(0.0784281, 0.333918, -1.33171);
    scene->getActiveCamera()->setViewUp(-0.00193505, -0.32117, -0.947019);
    scene->getConfig()->debugCamBoundingBox = false;
    *scene->getCamera("debug") = *scene->getActiveCamera();

    auto planeObject = std::make_shared<CollidingObject>();
    {
        auto plane = std::make_shared<Plane>(Vec3d(0.0, -0.05, -1.2), Vec3d(0.0, 1.0, 0.0));
        //planeObject->setVisualGeometry(plane);
        planeObject->setCollidingGeometry(plane);
    }
    scene->addSceneObject(planeObject);

    auto                            pbdModel  = std::make_shared<PbdModel>();
    std::shared_ptr<PbdModelConfig> pbdConfig = pbdModel->getConfig();
    pbdConfig->m_gravity = Vec3d(0.0, -8.0, 0.0);
    //pbdConfig->m_gravity = Vec3d(0.0, 0.0, 0.0);
    pbdConfig->m_dt = simTimestep;
    pbdConfig->m_iterations = iterations;
    pbdConfig->m_linearDampingCoeff  = globalLinearDamping;
    pbdConfig->m_angularDampingCoeff = globalAngularDamping;
    pbdConfig->m_doPartitioning      = false;

    std::shared_ptr<PbdObject> tissueObj = makeTissueObj("Tissue", pbdModel);
    scene->addSceneObject(tissueObj);

    auto tissuePlaneCollision = std::make_shared<PbdObjectCollision>(tissueObj, planeObject);
    tissuePlaneCollision->setDeformableStiffnessA(0.3);
    //tissuePlaneCollision->setDeformableStiffnessB();
    scene->addInteraction(tissuePlaneCollision);

    std::shared_ptr<PbdObject> leftToolObj = makeLapToolObj("leftToolObj", pbdModel);
    scene->addSceneObject(leftToolObj);
    std::shared_ptr<PbdObject> rightToolObj = makeLapToolObj("rightToolObj", pbdModel);
    scene->addSceneObject(rightToolObj);

    std::vector<std::shared_ptr<PbdObject>> esophagusCapsuleObjs;
    makeCapsules("CapsuleObjs", esophagusCapsuleObjs, pbdModel, 0.005);
    for (int i = 0; i < esophagusCapsuleObjs.size(); i++)
    {
        scene->addSceneObject(esophagusCapsuleObjs[i]);
        auto capsCollision = std::make_shared<PbdObjectCollision>(
            tissueObj, esophagusCapsuleObjs[i]);
        capsCollision->setRestitution(0.0);
        capsCollision->setFriction(0.0);
        capsCollision->setUseCorrectVelocity(false);
        capsCollision->setRigidBodyCompliance(0.00000001);
        scene->addInteraction(capsCollision);
    }

    // Add picking interaction for both jaws of the tool
    auto leftGrasping = std::make_shared<PbdObjectGrasping>(tissueObj);
    leftGrasping->setGeometryToPick(tissueObj->getCollidingGeometry(),
        std::dynamic_pointer_cast<PointwiseMap>(tissueObj->getPhysicsToCollidingMap()));
    leftGrasping->setStiffness(graspStiffness);
    scene->addInteraction(leftGrasping);
    auto rightGrasping = std::make_shared<PbdObjectGrasping>(tissueObj);
    rightGrasping->setGeometryToPick(tissueObj->getCollidingGeometry(),
        std::dynamic_pointer_cast<PointwiseMap>(tissueObj->getPhysicsToCollidingMap()));
    rightGrasping->setStiffness(graspStiffness);
    scene->addInteraction(rightGrasping);

    auto leftToolCollision = std::make_shared<PbdObjectCollision>(tissueObj, leftToolObj);
    leftToolCollision->setFriction(0.5);
    leftToolCollision->setRestitution(0.0);
    leftToolCollision->setRigidBodyCompliance(0.0000001);
    scene->addInteraction(leftToolCollision);
    /* auto leftToolFatCollision = std::make_shared<PbdObjectCollision>(fatObj, leftToolObj);
     leftToolFatCollision->setFriction(0.0);
     leftToolFatCollision->setRestitution(0.0);
     leftToolFatCollision->setRigidBodyCompliance(0.0000001);
     scene->addInteraction(leftToolFatCollision);*/
    auto rightToolCollision = std::make_shared<PbdObjectCollision>(tissueObj, rightToolObj);
    rightToolCollision->setFriction(0.5);
    rightToolCollision->setRestitution(0.0);
    rightToolCollision->setRigidBodyCompliance(0.0000001);
    scene->addInteraction(rightToolCollision);

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

        auto hapticManager = DeviceManagerFactory::makeDeviceManager();

        auto driver = std::make_shared<SimulationManager>();
        driver->addModule(viewer);
        driver->addModule(sceneManager);
        driver->addModule(hapticManager);
        driver->setDesiredDt(executionTimestep);

        std::shared_ptr<DeviceClient> leftHapticDevice = hapticManager->makeDeviceClient("LeftDevice");
        auto                          leftController   = leftToolObj->getComponent<PbdObjectController>();
        leftController->setDevice(leftHapticDevice);
        leftController->setTranslationOffset(Vec3d(0.0, 0.0, -1.2));

        std::shared_ptr<DeviceClient> rightHapticDevice = hapticManager->makeDeviceClient("RightDevice");

        //dummyClient->setOrientation(Quatd::FromTwoVectors(Vec3d(0.0, 0.0, 1.0), Vec3d(0.25, 2.0, 0.5).normalized()));
        auto rightController = rightToolObj->getComponent<PbdObjectController>();
        rightController->setDevice(rightHapticDevice);
        rightController->setTranslationScaling(0.5);
        rightController->setTranslationOffset(Vec3d(0.0, 0.0, -1.2));

        auto rightPortHole = rightToolObj->addComponent<PortHoleInteraction>();
        rightPortHole->setTool(rightToolObj);

        (*rightToolObj->getPbdBody()->vertices)[0] = Vec3d(rightPortHolePos);
        rightPortHole->setPortHoleLocation(rightPortHolePos);
        rightPortHole->setToolGeometry(rightToolObj->getCollidingGeometry());
        rightPortHole->setCompliance(0.00000001);
        auto sphere = std::make_shared<Sphere>(rightPortHolePos, 0.005);
        auto rightPortVisuals = rightToolObj->addComponent<VisualModel>();
        rightPortVisuals->setGeometry(sphere);

        auto leftPortHole = leftToolObj->addComponent<PortHoleInteraction>();
        leftPortHole->setTool(leftToolObj);
        (*leftToolObj->getPbdBody()->vertices)[0] = Vec3d(leftPortHolePos);
        leftPortHole->setPortHoleLocation(leftPortHolePos);
        leftPortHole->setToolGeometry(leftToolObj->getCollidingGeometry());
        leftPortHole->setCompliance(0.00000001);
        auto sphere2 = std::make_shared<Sphere>(leftPortHolePos, 0.005);
        auto leftPortVisuals = leftToolObj->addComponent<VisualModel>();
        leftPortVisuals->setGeometry(sphere2);

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
                            leftToolCollision->setEnabled(false);
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
                            leftToolCollision->setEnabled(true);
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
                            rightToolCollision->setEnabled(false);
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
                            rightToolCollision->setEnabled(true);
                        }
                    }
                });

            // Add default mouse and keyboard controls to the viewer
            std::shared_ptr<Entity> mouseAndKeyControls =
                SimulationUtils::createDefaultSceneControl(driver);
            scene->addSceneObject(mouseAndKeyControls);
        }

        driver->start();
    }
}