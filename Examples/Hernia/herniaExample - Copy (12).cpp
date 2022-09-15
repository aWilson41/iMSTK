/*
** This file is part of the Interactive Medical Simulation Toolkit (iMSTK)
** iMSTK is distributed under the Apache License, Version 2.0.
** See accompanying NOTICE for details.
*/

#include "imstkCamera.h"
#include "imstkCapsule.h"
#include "imstkCompositeImplicitGeometry.h"
#include "imstkDeviceManager.h"
#include "imstkDeviceManagerFactory.h"
#include "imstkDirectionalLight.h"
#include "imstkDummyClient.h"
#include "imstkGeometryUtilities.h"
#include "imstkImplicitGeometryToImageData.h"
#include "imstkKeyboardDeviceClient.h"
#include "imstkKeyboardSceneControl.h"
#include "imstkLineMesh.h"
#include "imstkLogger.h"
#include "imstkMeshIO.h"
#include "imstkMouseDeviceClient.h"
#include "imstkMouseSceneControl.h"
#include "imstkOrientedBox.h"
#include "imstkPbdContactConstraint.h"
#include "imstkPbdAngularConstraint.h"
#include "imstkPbdModel.h"
#include "imstkPbdModelConfig.h"
#include "imstkPbdObject.h"
#include "imstkPbdObjectCollision.h"
#include "imstkPbdObjectController.h"
#include "imstkPbdObjectGrasping.h"
#include "imstkPlane.h"
#include "imstkPointwiseMap.h"
#include "imstkRenderMaterial.h"
#include "imstkScene.h"
#include "imstkSceneManager.h"
#include "imstkSimulationManager.h"
#include "imstkSimulationUtils.h"
#include "imstkSphere.h"
#include "imstkSurfaceMeshFlyingEdges.h"
#include "imstkSurfaceMeshSubdivide.h"
#include "imstkVisualModel.h"
#include "imstkVTKViewer.h"
#include "imstkImageData.h"
#include "imstkDummyClient.h"
#include "imstkSignedDistanceField.h"
#include "imstkSurfaceMeshSubdivide.h"
#include "imstkPortHoleInteraction.h"
#include "imstkPbdSolver.h"

using namespace imstk;

#define HERNIA_RESOURCE_DIR "E:/Hernia/vlahhs/resources/"

double tissueCapsuleCompliance = 0.0;

double capsuleMass = 1.0;
double capsuleInertiaScale = 0.001;
double capsuleDistConstraintStiffness = 1000.0;
//double capsuleHingeConstraintCompliance = 0.1;

double tissueParticleMass = 0.1;
double tissueDistStiffness = 10000.0;
double tissueDihedralStiffness = 1.0;

double globalLinearDamping = 0.05;
double globalAngularDamping = 0.01;

int iterations = 2;
int collisionIterations = 1;
double simTimestep = 0.001;
double executionTimestep = 0.0025;

double graspStiffness = 0.01;

double tetYoungsModulus = 5000.0;
double tetPossionsRatio = 0.46;

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

// Stability:
// XPBD appears to handle multi contact poorly (not as well as pbd). It's hard for me to
// fully understand. In pbd if you have multi contact you can drop stiffness and you will
// get smaller steps towards the solution. It's great not to step immediately to the solution
// particularly when multiple constraints are involved. But xpbd is parameterized by dt such
// that it is consistent in time. Say the solution of x(t), x(1.4)=5 and x(1.5)=10. With dt=1
// Every step we do some 10 iterations. There is one step here at dt=1 so we do 10. Then what
// prevents us from jumping to 10 immediately from 5 on the first iteration. This is bad. I'd
// almost want something that parameterizes over the iteration count sucn that 10 is reached at
// the final iteration. In xpdd though this can be thought of exactly like taking smaller
// timesteps though. This is the motivation of the substeps PBD paper from mueller. Why not
// always just use 1 iteration and lower dt.
// 
// The problem, hard to quantify, is likely the other computation in the scene. PBD does
// reprojection, which is a nifty trick to reach better solutions with less computation.
// Other parts of the scene (namely collision detection and rendering) are not being run
// with extra iterations in the solver. We already do this with rendering and substeps.
// We can take say 1000 physics steps for every 1 render. A similar thing might be required
// with collision detection. Running it at a different rate to allow/improve the usage of single
// iterations in xpbd.

///
/// \brief Compute the a vector of points indices that are coincident with the on the parent
///
static std::vector<int>
computeFixedPtsViaMap(std::shared_ptr<PointSet> parent,
    std::shared_ptr<PointSet> child,
    const double tolerance = 0.00001)
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
/// \brief Creates tissue object
/// \param name
/// \param physical dimension of tissue
/// \param dimensions of tetrahedral grid used for tissue
/// \param center of tissue block
///
static std::shared_ptr<PbdObject>
makeTissueObj(const std::string& name,
    std::shared_ptr<PbdModel> model)
{
    // Setup the Geometry
    /*std::shared_ptr<SurfaceMesh> tissueMesh =
        MeshIO::read<SurfaceMesh>("C:/Users/Andx_/Desktop/NewHernia/stomach_surface.obj");*/
    std::shared_ptr<TetrahedralMesh> tissueMesh =
        MeshIO::read<TetrahedralMesh>("C:/Users/Andx_/Desktop/NewHernia/volume/stomach_surface_tmp_.msh");
    std::shared_ptr<SurfaceMesh> tissueSurfMesh = tissueMesh->extractSurfaceMesh();
    std::shared_ptr<SurfaceMesh> colMesh = MeshIO::read<SurfaceMesh>(
        "C:/Users/Andx_/Desktop/NewHernia/volume/stomach_surface_collision.obj");

    // Setup the material
    auto material = std::make_shared<RenderMaterial>();
    material->setColor(Color(93.0 / 255.0, 38.0 / 255.0, 37.0 / 255.0));
    //material->setDisplayMode(RenderMaterial::DisplayMode::WireframeSurface);
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

    auto fixedVerts = MeshIO::read<PointSet>("C:/Users/Andx_/Desktop/NewHernia/volume/fixed_verts.obj");
    PointwiseMap fixedMapper;
    fixedMapper.setParentGeometry(tissueObj->getPhysicsGeometry());
    fixedMapper.setChildGeometry(fixedVerts);
    fixedMapper.compute();
    const std::unordered_map<int, int>& fixedMap = fixedMapper.getMap();
    for (auto i : fixedMap)
    {
        tissueObj->getPbdBody()->fixedNodeIds.push_back(i.second);
    }

    // Setup the Parameters
  /*  model->getConfig()->enableConstraint(PbdModelConfig::ConstraintGenType::Distance, tissueDistStiffness,
        tissueObj->getPbdBody()->bodyHandle);
    model->getConfig()->enableConstraint(PbdModelConfig::ConstraintGenType::Dihedral, tissueDihedralStiffness,
        tissueObj->getPbdBody()->bodyHandle);*/
    model->getConfig()->enableFemConstraint(PbdFemConstraint::MaterialType::NeoHookean,
        tissueObj->getPbdBody()->bodyHandle);
    model->getConfig()->m_femParams->m_PoissonRatio = tetPossionsRatio;
    model->getConfig()->m_femParams->m_YoungModulus = tetYoungsModulus;

    //auto esophagusTetMesh = MeshIO::read<TetrahedralMesh>("C:/Users/Andx_/Desktop/stomach_esophagus.vtk");

    // A set of vertices on the esophagus are constrained

    return tissueObj;
}

///
/// \brief Creates and PbdObject that links/connects to the tissueObj provided
/// via computation of overlapping points
/// 
static std::shared_ptr<PbdObject>
makeLigamentObj(std::string name,
    std::string ligamentGeomFileName,
    std::string fixedVertexFileName,
    const double stiffness,
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
    model->getConfig()->addPbdConstraintFunctor([=](PbdConstraintContainer& container)
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
makeCapsules(const std::string& name,
    std::vector<std::shared_ptr<PbdObject>>& capsuleObjs,
    std::shared_ptr<PbdModel> model,
    const double capsuleRadius)
{
    capsuleObjs.clear();

    // Read in a medial line
    auto                                     maLineMesh =
        MeshIO::read<LineMesh>("C:/Users/Andx_/Desktop/NewHernia/capsules.obj");
    std::shared_ptr<VecDataArray<double, 3>> maVerticesPtr = maLineMesh->getVertexPositions();
    const VecDataArray<double, 3>& maVertices = *maVerticesPtr;
    std::shared_ptr<VecDataArray<int, 2>>    maIndicesPtr = maLineMesh->getCells();
    const VecDataArray<int, 2>& maIndices = *maIndicesPtr;

    // For every line segment in the line mesh create a rigid body capsule
    for (int i = 0; i < maIndices.size(); i++)
    {
        const Vec2i& cell = maIndices[i];
        const Vec3d  p = maVertices[cell[0]];
        const Vec3d  q = maVertices[cell[1]];
        const Vec3d  center = (q + p) * 0.5;
        const Vec3d  diff = q - p;
        const Vec3d  dir = diff.normalized();

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
        capsuleObj->getPbdBody()->setRigid(center, // Position
            capsuleMass,                                  // Mass
            orientation,                           // Orientation
            Mat3d::Identity() * capsuleInertiaScale);            // Inertia
        model->getConfig()->setBodyDamping(capsuleObj->getPbdBody()->bodyHandle, 0.1, 0.0);

        capsuleObjs.push_back(capsuleObj);
    }

    // This one just connects every single capsule to be fixed in place in positin
    // and orientation
    model->getConfig()->addPbdConstraintFunctor(
        [=](PbdConstraintContainer& container)
        {
            // For every capsule
            for (int i = 0; i < maIndices.size(); i++)
            {
                auto capsule = std::dynamic_pointer_cast<Capsule>(capsuleObjs[i]->getVisualGeometry());
                const Vec3d center = (*capsuleObjs[i]->getPbdBody()->vertices)[0];
                const Quatd orientation = (*capsuleObjs[i]->getPbdBody()->orientations)[0];
                const Vec3d dir = orientation._transformVector(Vec3d(0.0, 1.0, 0.0));
                const PbdParticleId fixedParticle =
                    model->addVirtualParticle(
                        center, // Position
                        orientation, // Orientation
                        0.0, // Mass
                        Mat3d::Identity(), // Inertia
                        Vec3d::Zero(), // Velocity
                        Vec3d::Zero(), // Angular Velocity
                        true); // Persist

                auto distConstraint = std::make_shared<PbdDistanceConstraint>();
                distConstraint->initConstraint(
                    0.0,
                    { capsuleObjs[i]->getPbdBody()->bodyHandle, 0 },
                    fixedParticle, capsuleDistConstraintStiffness);
                container.addConstraint(distConstraint);

                /*auto angularDistConstraint = std::make_shared<PbdAngularDistanceConstraint>();
                angularDistConstraint->initConstraint(
                    { capsuleObjs[i]->getPbdBody()->bodyHandle, 0 },
                    fixedParticle,
                    capsuleHingeConstraintCompliance);
                container.addConstraint(angularDistConstraint);*/
            }


            // Also add capsule->capsule joint constraints to keep them connected
            for (int i = 0; i < maIndices.size() - 1; i++)
            {
                // Assumes consistent directions of the capsules
                auto capsule0 = std::dynamic_pointer_cast<Capsule>(capsuleObjs[i]->getVisualGeometry());
                Vec3d jp0 = Vec3d::Zero();
                {
                    const Vec3d& capsulePos = capsule0->getPosition();
                    const double capsuleLength = capsule0->getLength();
                    const Quatd& capsuleOrientation = capsule0->getOrientation();
                    const Vec3d& capsulePosA = capsulePos - 0.5 * capsuleLength * capsuleOrientation.toRotationMatrix().col(1);
                    const Vec3d& capsulePosB = capsulePos + (capsulePos - capsulePosA);
                    jp0 = capsulePosB;
                }
                auto capsule1 = std::dynamic_pointer_cast<Capsule>(capsuleObjs[i + 1]->getVisualGeometry());
                Vec3d jp1 = Vec3d::Zero();
                {
                    const Vec3d& capsulePos = capsule1->getPosition();
                    const double capsuleLength = capsule1->getLength();
                    const Quatd& capsuleOrientation = capsule1->getOrientation();
                    const Vec3d& capsulePosA = capsulePos - 0.5 * capsuleLength * capsuleOrientation.toRotationMatrix().col(1);
                    const Vec3d& capsulePosB = capsulePos + (capsulePos - capsulePosA);
                    jp1 = capsulePosA;
                }
                const Vec3d jointPos = (jp0 + jp1) * 0.5;

                auto bodyToBodyConstraint = std::make_shared<PbdBodyToBodyDistanceConstraint>();
                bodyToBodyConstraint->initConstraint(model->getBodies(),
                    { capsuleObjs[i]->getPbdBody()->bodyHandle, 0 },
                    jointPos,
                    { capsuleObjs[i + 1]->getPbdBody()->bodyHandle, 0 },
                    jointPos, 0.000001);
                container.addConstraint(bodyToBodyConstraint);
            }
            // Fix the two ends
            {
                // Assumes consistent directions of the capsules
                auto capsule0 = std::dynamic_pointer_cast<Capsule>(capsuleObjs[0]->getVisualGeometry());
                Vec3d jp0 = Vec3d::Zero();
                {
                    const Vec3d& capsulePos = capsule0->getPosition();
                    const double capsuleLength = capsule0->getLength();
                    const Quatd& capsuleOrientation = capsule0->getOrientation();
                    const Vec3d& capsulePosA = capsulePos - 0.5 * capsuleLength * capsuleOrientation.toRotationMatrix().col(1);
                    const Vec3d& capsulePosB = capsulePos + (capsulePos - capsulePosA);
                    jp0 = capsulePosA;
                }
                // Assumes consistent directions of the capsules
                auto capsule1 = std::dynamic_pointer_cast<Capsule>(capsuleObjs.back()->getVisualGeometry());
                Vec3d jp1 = Vec3d::Zero();
                {
                    const Vec3d& capsulePos = capsule1->getPosition();
                    const double capsuleLength = capsule1->getLength();
                    const Quatd& capsuleOrientation = capsule1->getOrientation();
                    const Vec3d& capsulePosA = capsulePos - 0.5 * capsuleLength * capsuleOrientation.toRotationMatrix().col(1);
                    const Vec3d& capsulePosB = capsulePos + (capsulePos - capsulePosA);
                    jp1 = capsulePosB;
                }

                const PbdParticleId& fixedVid0 = model->addVirtualParticle(jp0, 0.0, Vec3d::Zero(), true);
                auto bodyToBodyConstraint0 = std::make_shared<PbdBodyToBodyDistanceConstraint>();
                bodyToBodyConstraint0->initConstraint(model->getBodies(),
                    { capsuleObjs[0]->getPbdBody()->bodyHandle, 0 },
                    jp0,
                    fixedVid0,
                    0.0, 0.00000001);
                container.addConstraint(bodyToBodyConstraint0);

                const PbdParticleId& fixedVid1 = model->addVirtualParticle(jp1, 0.0, Vec3d::Zero(), true);
                auto bodyToBodyConstraint1 = std::make_shared<PbdBodyToBodyDistanceConstraint>();
                bodyToBodyConstraint1->initConstraint(model->getBodies(),
                    { capsuleObjs.back()->getPbdBody()->bodyHandle, 0 },
                    jp1,
                    fixedVid1,
                    0.0, 0.00000001);
                container.addConstraint(bodyToBodyConstraint1);
            }
            {
            }
        });
}
void
makeCollidingCapsules(const std::string& name,
    std::vector<std::shared_ptr<CollidingObject>>& capsuleObjs,
    std::shared_ptr<PbdModel> model,
    const double capsuleRadius)
{
    capsuleObjs.clear();

    // Read in a medial line
    auto                                     maLineMesh =
        MeshIO::read<LineMesh>("C:/Users/Andx_/Desktop/NewHernia/capsules.obj");
    std::shared_ptr<VecDataArray<double, 3>> maVerticesPtr = maLineMesh->getVertexPositions();
    const VecDataArray<double, 3>& maVertices = *maVerticesPtr;
    std::shared_ptr<VecDataArray<int, 2>>    maIndicesPtr = maLineMesh->getCells();
    const VecDataArray<int, 2>& maIndices = *maIndicesPtr;

    // For every line segment in the line mesh create a rigid body capsule
    for (int i = 0; i < maIndices.size(); i++)
    {
        const Vec2i& cell = maIndices[i];
        const Vec3d  p = maVertices[cell[0]];
        const Vec3d  q = maVertices[cell[1]];
        const Vec3d  center = (q + p) * 0.5;
        const Vec3d  diff = q - p;
        const Vec3d  dir = diff.normalized();

        const Quatd orientation = Quatd::FromTwoVectors(Vec3d(0.0, 1.0, 0.0), dir);
        auto capsule = std::make_shared<Capsule>(center, capsuleRadius, diff.norm(), orientation);
        //auto capsule = std::make_shared<Sphere>(Vec3d::Zero(), capsuleRadius);

        auto capsuleObj = std::make_shared<CollidingObject>(name);
        capsuleObj->setVisualGeometry(capsule);
        capsuleObj->setCollidingGeometry(capsule);
        //capsuleObj->getVisualModel(0)->getRenderMaterial()->setColor(Color::Red);
        //capsuleObj->getVisualModel(0)->getRenderMaterial()->setOpacity(0.5);

        capsuleObjs.push_back(capsuleObj);
    }
}

void
makeFatCapsules(const std::string& name,
    std::vector<std::shared_ptr<CollidingObject>>& capsuleObjs,
    std::shared_ptr<PbdModel> model,
    const double radius)
{
    capsuleObjs.clear();

    auto lineMesh = MeshIO::read<LineMesh>(HERNIA_RESOURCE_DIR "/FundoTest4/fat_edge.obj");

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

        auto capsule = std::make_shared<Capsule>(
            Vec3d::Zero(),
            radius,
            diff.norm(),
            Quatd::Identity());

        auto capsuleObj = std::make_shared<PbdObject>(name);
        capsuleObj->setDynamicalModel(model);
        capsuleObj->setVisualGeometry(capsule);
        capsuleObj->setPhysicsGeometry(capsule);
        capsuleObj->setCollidingGeometry(capsule);
        capsuleObj->getVisualModel(0)->getRenderMaterial()->setColor(Color::Red);
        //capsuleObj->getVisualModel(0)->getRenderMaterial()->setOpacity(0.5);

        // Setup body
        const Quatd orientation = Quatd::FromTwoVectors(Vec3d(0.0, 1.0, 0.0), dir);
        capsuleObj->getPbdBody()->setRigid(
            center, // Position
            capsuleMass, // Mass
            Quatd::FromTwoVectors(Vec3d(0.0, 1.0, 0.0), dir), // Orientation
            Mat3d::Identity() * capsuleInertiaScale); // Inertia
        model->getConfig()->setBodyDamping(capsuleObj->getPbdBody()->bodyHandle, 0.1, 0.0);

        capsuleObjs.push_back(capsuleObj);
    }
}

std::shared_ptr<PbdObject>
makeLapToolObj(const std::string& name,
    std::shared_ptr<PbdModel> model)
{
    auto lapTool = std::make_shared<PbdObject>(name);

    const double capsuleLength = 0.3;
    auto         toolGeom = std::make_shared<Capsule>(Vec3d(0.0, 0.0, capsuleLength * 0.5),
        0.002, capsuleLength,
        Quatd::FromTwoVectors(Vec3d(0.0, 1.0, 0.0), Vec3d(0.0, 0.0, 1.0)));

    lapTool->setDynamicalModel(model);
    lapTool->setPhysicsGeometry(toolGeom);
    lapTool->setCollidingGeometry(toolGeom);
    lapTool->setVisualGeometry(toolGeom);

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
        MeshIO::read<ImageData>(HERNIA_RESOURCE_DIR "/FundoTest1/stomach_container_SDF.mhd")
        ->cast(IMSTK_DOUBLE);
    diaphgramObj->setCollidingGeometry(std::make_shared<SignedDistanceField>(diaphgramSdfImage));

    /*imstkNew<SurfaceMeshFlyingEdges> isoExtract;
    isoExtract->setInputImage(diaphgramSdfImage);
    isoExtract->update();
    diaphgramObj->setVisualGeometry(isoExtract->getOutputMesh());
    diaphgramObj->getVisualModel(0)->getRenderMaterial()->setOpacity(0.5);*/

    auto surfMesh = MeshIO::read<SurfaceMesh>(HERNIA_RESOURCE_DIR "/diaphragm.obj");
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
    material->addTexture(std::make_shared<Texture>(HERNIA_RESOURCE_DIR "/textures/Diaphragm_SG1.png",
        Texture::Type::Diffuse));
    material->addTexture(std::make_shared<Texture>(
        HERNIA_RESOURCE_DIR "/textures/Diaphragm_SG1_normal.png", Texture::Type::Normal));
    material->addTexture(std::make_shared<Texture>(
        HERNIA_RESOURCE_DIR "/textures/Diaphragm_SG1_normal.png", Texture::Type::CoatNormal));
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

    auto pbdModel = std::make_shared<PbdModel>();
    std::shared_ptr<PbdModelConfig> pbdConfig = pbdModel->getConfig();
    pbdConfig->m_gravity = Vec3d(0.0, 0.0, 0.0);
    //pbdConfig->m_gravity = Vec3d(0.0, -20.0, 0.0);
    //pbdConfig->m_gravity = Vec3d(9.8, 0.0, 0.0);
    pbdConfig->m_dt = simTimestep;
    pbdConfig->m_iterations = iterations;
    pbdConfig->m_collisionIterations = collisionIterations;
    pbdConfig->m_linearDampingCoeff = globalLinearDamping;
    pbdConfig->m_angularDampingCoeff = globalAngularDamping;
    pbdConfig->m_doPartitioning = false;
    //pbdConfig->m_collisionSolverType = PbdConstraint::SolverType::PBD;

    std::shared_ptr<PbdObject> tissueObj = makeTissueObj("Tissue", pbdModel);
    scene->addSceneObject(tissueObj);

    std::shared_ptr<PbdObject> leftToolObj = makeLapToolObj("leftToolObj", pbdModel);
    scene->addSceneObject(leftToolObj);
    std::shared_ptr<PbdObject> rightToolObj = makeLapToolObj("rightToolObj", pbdModel);
    scene->addSceneObject(rightToolObj);
    std::shared_ptr<PbdObject> cameraToolObj = makeLapToolObj("cameraToolObj", pbdModel);
    (*cameraToolObj->getPbdBody()->vertices)[0] = Vec3d(0.001, 0.078, -1.206);
    (*cameraToolObj->getPbdBody()->orientations)[0] = Quatd(0.826, -0.54, -0.149, 0.061);
    scene->addSceneObject(cameraToolObj);

    std::vector<std::shared_ptr<PbdObject>> esophagusCapsuleObjs;
    makeCapsules("CapsuleObjs", esophagusCapsuleObjs, pbdModel, 0.002);
    for (int i = 0; i < esophagusCapsuleObjs.size(); i++)
    {
        scene->addSceneObject(esophagusCapsuleObjs[i]);
        auto capsCollision = std::make_shared<PbdObjectCollision>(
            tissueObj, esophagusCapsuleObjs[i]);
        capsCollision->setRestitution(0.0);
        capsCollision->setFriction(0.0);
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
         HERNIA_RESOURCE_DIR "/FundoTest1/lesser_omentum_sim.obj",
         HERNIA_RESOURCE_DIR "/FundoTest1/lesser_omentum_sim_fixed.obj",
         1000.0,
         tissueObj);
     scene->addSceneObject(lesserOmentumObj);

     std::shared_ptr<PbdObject> greaterOmentumObj =
         makeLigamentObj("greater_omentum",
             HERNIA_RESOURCE_DIR "/FundoTest4/greater_omentum_sim.obj",
             HERNIA_RESOURCE_DIR "/FundoTest4/greater_omentum_sim_fixed.obj",
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

     //auto fatObj = makeLigamentObj("fat",
     //    HERNIA_RESOURCE_DIR "/FundoTest4/fat.obj",
     //    HERNIA_RESOURCE_DIR "/FundoTest4/fat_fixed.obj",
     //    1000.0,
     //    tissueObj);
     //{
     //    // Swap the collision geometry
     //    auto cdLineMesh = MeshIO::read<LineMesh>(HERNIA_RESOURCE_DIR "/FundoTest4/fat_edge.obj");
     //    fatObj->setCollidingGeometry(cdLineMesh);

     //    // Setup a mapping between physics geometry and collision
     //    auto colMap = std::make_shared<PointwiseMap>(fatObj->getPhysicsGeometry(), cdLineMesh);
     //    colMap->setTolerance(0.0001);
     //    fatObj->setPhysicsToCollidingMap(colMap);
     //}
     //scene->addSceneObject(fatObj);
#pragma endregion

    // Add picking interaction for both jaws of the tool
    auto leftGrasping = std::make_shared<PbdObjectGrasping>(tissueObj);
    /*leftGrasping->setGeometryToPick(tissueObj->getCollidingGeometry(),
        std::dynamic_pointer_cast<PointwiseMap>(tissueObj->getPhysicsToCollidingMap()));*/
    leftGrasping->setStiffness(graspStiffness);
    scene->addInteraction(leftGrasping);
    auto rightGrasping = std::make_shared<PbdObjectGrasping>(tissueObj);
   /* rightGrasping->setGeometryToPick(tissueObj->getCollidingGeometry(),
        std::dynamic_pointer_cast<PointwiseMap>(tissueObj->getPhysicsToCollidingMap()));*/
    rightGrasping->setStiffness(graspStiffness);
    scene->addInteraction(rightGrasping);

    /*auto leftToolCollision = std::make_shared<PbdObjectCollision>(tissueObj, leftToolObj);
    leftToolCollision->setFriction(0.0);
    leftToolCollision->setRestitution(0.0);
    leftToolCollision->setRigidBodyCompliance(0.0000001);
    scene->addInteraction(leftToolCollision);*/
    /* auto leftToolFatCollision = std::make_shared<PbdObjectCollision>(fatObj, leftToolObj);
     leftToolFatCollision->setFriction(0.0);
     leftToolFatCollision->setRestitution(0.0);
     leftToolFatCollision->setRigidBodyCompliance(0.0000001);
     scene->addInteraction(leftToolFatCollision);*/
     /*auto rightToolCollision = std::make_shared<PbdObjectCollision>(tissueObj, rightToolObj);
     rightToolCollision->setFriction(0.0);
     rightToolCollision->setRestitution(0.0);
     rightToolCollision->setRigidBodyCompliance(0.0000001);
     scene->addInteraction(rightToolCollision);*/

    auto diaphgramCollision = std::make_shared<PbdObjectCollision>(tissueObj, diaphgramObj);
    diaphgramCollision->setFriction(0.0);
    diaphgramCollision->setRestitution(0.0);
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

        auto driver = std::make_shared<SimulationManager>();
        driver->addModule(viewer);
        driver->addModule(sceneManager);
        driver->setDesiredDt(executionTimestep);

        // Use haptic device for the left tool (doing more complex movements)
        auto leftHapticManager = DeviceManagerFactory::makeDeviceManager("OpenHapticDeviceManager");
        driver->addModule(leftHapticManager);
        std::shared_ptr<DeviceClient> hapticDevice = leftHapticManager->makeDeviceClient();
        leftHapticManager->init();

        auto leftController = leftToolObj->getComponent<PbdObjectController>();
        leftController->setDevice(hapticDevice);
        leftController->setTranslationOffset(Vec3d(0.0, 0.0, -1.3));

        // Due to an issue with Haply and OpenHaptics, OpenHaptics is initialized first manually, then a wait, then Haply
        std::this_thread::sleep_for(std::chrono::seconds(1));

        // Use Haply device for right tool
        auto rightHapticManager = DeviceManagerFactory::makeDeviceManager("HaplyDeviceManager");
        driver->addModule(rightHapticManager);
        std::shared_ptr<DeviceClient> rightHapticDevice = rightHapticManager->makeDeviceClient();
        rightHapticManager->init();

        //dummyClient->setOrientation(Quatd::FromTwoVectors(Vec3d(0.0, 0.0, 1.0), Vec3d(0.25, 2.0, 0.5).normalized()));
        auto rightController = rightToolObj->getComponent<PbdObjectController>();
        rightController->setDevice(rightHapticDevice);
        rightController->setTranslationScaling(0.5);
        rightController->setTranslationOffset(Vec3d(0.1, 0.0, -1.23));

        auto camController = cameraToolObj->getComponent<PbdObjectController>();
        camController->setTranslationOffset(Vec3d(0., 0.0, -1.23));

        auto rightPortHole = std::make_shared<PortHoleInteraction>(rightToolObj);
        const Vec3d rightPortHolePos = Vec3d(0.015, 0.092, -1.117);
        (*rightToolObj->getPbdBody()->vertices)[0] = Vec3d(rightPortHolePos);
        rightPortHole->setPortHoleLocation(rightPortHolePos);
        auto sphere = std::make_shared<Sphere>(rightPortHolePos, 0.005);
        rightPortHole->setVisualGeometry(sphere);
        rightPortHole->setToolGeometry(rightToolObj->getCollidingGeometry());
        rightPortHole->setCompliance(0.00000001);
        scene->addInteraction(rightPortHole);

        auto leftPortHole = std::make_shared<PortHoleInteraction>(leftToolObj);
        const Vec3d leftPortHolePos = Vec3d(-0.065, 0.078, -1.127);
        (*leftToolObj->getPbdBody()->vertices)[0] = Vec3d(leftPortHolePos);
        leftPortHole->setPortHoleLocation(leftPortHolePos);
        auto sphere2 = std::make_shared<Sphere>(leftPortHolePos, 0.005);
        leftPortHole->setVisualGeometry(sphere2);
        leftPortHole->setToolGeometry(leftToolObj->getCollidingGeometry());
        leftPortHole->setCompliance(0.00000001);
        scene->addInteraction(leftPortHole);

        auto cameraPortHole = std::make_shared<PortHoleInteraction>(cameraToolObj);
        const Vec3d cameraPortHolePos = Vec3d(-0.004, 0.093, -1.2);
        //(*cameraToolObj->getPbdBody()->vertices)[0] = Vec3d(cameraPortHolePos);
        cameraPortHole->setPortHoleLocation(cameraPortHolePos);
        auto sphere3 = std::make_shared<Sphere>(cameraPortHolePos, 0.005);
        cameraPortHole->setVisualGeometry(sphere3);
        cameraPortHole->setToolGeometry(cameraToolObj->getCollidingGeometry());
        cameraPortHole->setCompliance(0.00000001);
        scene->addInteraction(cameraPortHole);

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
                    auto capsule = std::dynamic_pointer_cast<Capsule>(cameraToolObj->getCollidingGeometry());
                    const Vec3d dir = -capsule->getRotation().col(2).normalized();
                    const Vec3d up = capsule->getRotation().col(1).normalized();
                    const double r = capsule->getRadius();
                    const double length = capsule->getLength();

                    scene->getCamera("default")->setPosition(capsule->getPosition() + dir * (length * 0.48));
                    scene->getCamera("default")->setFocalPoint(capsule->getPosition() + dir * (length * 0.5));
                    scene->getCamera("default")->setViewUp(up);
                });
            connect<ButtonEvent>(hapticDevice, &DeviceClient::buttonStateChanged,
                [&](ButtonEvent* e)
                {
                    if (e->m_buttonState == BUTTON_PRESSED)
                    {
                        if (e->m_button == 0)
                        {
                            LOG(INFO) << "Left Grasp!";
                            leftGrasping->beginVertexGrasp(std::dynamic_pointer_cast<AnalyticalGeometry>(leftToolObj->getCollidingGeometry()));
                            //leftToolCollision->setEnabled(false);
                            /*leftGrasping->beginCellGrasp(
                                std::dynamic_pointer_cast<AnalyticalGeometry>(leftToolObj->getCollidingGeometry()));*/
                        }
                        else if (e->m_button == 1)
                        {
                            LOG(INFO) << "Right Grasp!";
                            rightGrasping->beginVertexGrasp(std::dynamic_pointer_cast<AnalyticalGeometry>(leftToolObj->getCollidingGeometry()));
                            //rightToolCollision->setEnabled(false);
                            /*rightGrasping->beginCellGrasp(
                                std::dynamic_pointer_cast<AnalyticalGeometry>(rightToolObj->getCollidingGeometry()));*/
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
                        else if (e->m_button == 1)
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
                SimulationUtils::createDefaultSceneControlEntity(driver);
            scene->addSceneObject(mouseAndKeyControls);
        }

        driver->start();
    }
}