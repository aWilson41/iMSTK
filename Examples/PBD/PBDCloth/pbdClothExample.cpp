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

#include "imstkCamera.h"
#include "imstkGeometryUtilities.h"
#include "imstkImageData.h"
#include "imstkKeyboardDeviceClient.h"
#include "imstkKeyboardSceneControl.h"
#include "imstkMeshIO.h"
#include "imstkMouseDeviceClient.h"
#include "imstkMouseSceneControl.h"
#include "imstkNew.h"
#include "imstkPbdModel.h"
#include "imstkPbdObject.h"
#include "imstkRenderMaterial.h"
#include "imstkScene.h"
#include "imstkSceneManager.h"
#include "imstkSimulationManager.h"
#include "imstkSurfaceMesh.h"
#include "imstkVisualModel.h"
#include "imstkVTKViewer.h"

#include "imstkLineMesh.h"
#include "imstkPbdDistanceConstraint.h"
#include "imstkPbdConstraintContainer.h"

using namespace imstk;

static void
setFabricTextures(std::shared_ptr<RenderMaterial> material)
{
    auto diffuseTex = MeshIO::read<ImageData>(iMSTK_DATA_ROOT "/textures/fabricDiffuse.jpg");
    material->addTexture(std::make_shared<Texture>(diffuseTex, Texture::Type::Diffuse));
    auto normalTex = MeshIO::read<ImageData>(iMSTK_DATA_ROOT "/textures/fabricNormal.jpg");
    material->addTexture(std::make_shared<Texture>(normalTex, Texture::Type::Normal));
    auto ormTex = MeshIO::read<ImageData>(iMSTK_DATA_ROOT "/textures/fabricORM.jpg");
    material->addTexture(std::make_shared<Texture>(ormTex, Texture::Type::ORM));
}

static void
setFleshTextures(std::shared_ptr<RenderMaterial> material)
{
    auto diffuseTex = MeshIO::read<ImageData>(iMSTK_DATA_ROOT "/textures/fleshDiffuse.jpg");
    material->addTexture(std::make_shared<Texture>(diffuseTex, Texture::Type::Diffuse));
    auto normalTex = MeshIO::read<ImageData>(iMSTK_DATA_ROOT "/textures/fleshNormal.jpg");
    material->addTexture(std::make_shared<Texture>(normalTex, Texture::Type::Normal));
    auto ormTex = MeshIO::read<ImageData>(iMSTK_DATA_ROOT "/textures/fleshORM.jpg");
    material->addTexture(std::make_shared<Texture>(ormTex, Texture::Type::ORM));
}

///
/// \brief Creates cloth object
/// \param name
/// \param cloth width
/// \param cloth height
/// \param cloth row count
/// \param cloth column count
///
static std::shared_ptr<PbdObject>
makeClothObj(const std::string& name,
             const double       width,
             const double       height,
             const int          rowCount,
             const int          colCount)
{
    imstkNew<PbdObject> clothObj(name);

    // Setup the Geometry
    std::shared_ptr<SurfaceMesh> clothMesh =
        GeometryUtils::toTriangleGrid(Vec3d::Zero(),
            Vec2d(10.0, 10.0), Vec2i(16, 16), Quatd::Identity(), 2.0);

    // Setup the Parameters
    imstkNew<PbdModelConfig> pbdParams;
    pbdParams->m_gravity        = Vec3d(0.0, -9.8, 0.0);
    pbdParams->m_dt             = 0.005;
    pbdParams->m_iterations     = 5;
    pbdParams->m_doPartitioning = false;

    // Setup the Model
    imstkNew<PbdModel> pbdModel;
    pbdModel->setModelGeometry(clothMesh);
    pbdModel->configure(pbdParams);

    // Setup the VisualModel
    imstkNew<RenderMaterial> material;
    material->setBackFaceCulling(false);
    material->setDisplayMode(RenderMaterial::DisplayMode::Surface);
    material->setShadingModel(RenderMaterial::ShadingModel::PBR);
    setFleshTextures(material);
    imstkNew<VisualModel> visualModel;
    visualModel->setGeometry(clothMesh);
    visualModel->setRenderMaterial(material);

    // Setup the Object
    clothObj->addVisualModel(visualModel);
    clothObj->setPhysicsGeometry(clothMesh);
    clothObj->setDynamicalModel(pbdModel);
    //clothObj->getPbdBody()->fixedNodeIds = { 0, colCount - 1 };
    clothObj->getPbdBody()->uniformMassValue = width * height / (rowCount * colCount);
    pbdParams->enableConstraint(PbdModelConfig::ConstraintGenType::Distance, 1.0e2,
        clothObj->getPbdBody()->bodyHandle);
    /*pbdParams->enableConstraint(PbdModelConfig::ConstraintGenType::Dihedral, 1.0e1,
        clothObj->getPbdBody()->bodyHandle);*/

    return clothObj;
}

///
/// \brief Create pbd string object
///
static std::shared_ptr<PbdObject>
makePbdString(
    const std::string& name,
    const Vec3d& pos, const Vec3d& dir, const int numVerts,
    const double stringLength, std::shared_ptr<PbdModel> model)
{
    imstkNew<PbdObject> stringObj(name);

    // Setup the Geometry
    std::shared_ptr<LineMesh> stringMesh =
        GeometryUtils::toLineGrid(pos, dir, stringLength, numVerts);

    // Setup the VisualModel
    imstkNew<RenderMaterial> material;
    material->setBackFaceCulling(false);
    material->setColor(Color::Red);
    material->setLineWidth(2.0);
    material->setPointSize(6.0);
    material->setDisplayMode(RenderMaterial::DisplayMode::Wireframe);

    imstkNew<VisualModel> visualModel;
    visualModel->setGeometry(stringMesh);
    visualModel->setRenderMaterial(material);

    // Setup the Object
    stringObj->addVisualModel(visualModel);
    stringObj->setPhysicsGeometry(stringMesh);
    stringObj->setCollidingGeometry(stringMesh);
    stringObj->setDynamicalModel(model);
    stringObj->getPbdBody()->fixedNodeIds     = { 9 };
    stringObj->getPbdBody()->uniformMassValue = 100.0;//0.002 / numVerts; // grams
    model->getConfig()->enableConstraint(PbdModelConfig::ConstraintGenType::Distance, 10000.0,
        stringObj->getPbdBody()->bodyHandle);

    return stringObj;
}

///
/// \brief This example demonstrates the cloth simulation
/// using Position based dynamics
///
int
main()
{
    // Write log to stdout and file
    Logger::startLogger();

    // Setup a scene
    imstkNew<Scene> scene("PBDCloth");
    scene->getActiveCamera()->setFocalPoint(0.0, -5.0, 0.0);
    scene->getActiveCamera()->setPosition(0.0, 1.5, 25.0);
    scene->getActiveCamera()->setViewUp(0.0, 1.0, 0.0);

    std::shared_ptr<PbdObject> clothObj = makeClothObj("Cloth", 10.0, 10.0, 16, 16);
    scene->addSceneObject(clothObj);

    std::shared_ptr<PbdObject> lineObj = makePbdString("thread",
        Vec3d(0.0, 0.0, 0.0), Vec3d(1.0, 1.0, 0.0), 10, 5.0, clothObj->getPbdModel());
    scene->addSceneObject(lineObj);

    // Add a special constraint generator to add a constraint between the two bodies
    clothObj->getPbdModel()->getConfig()->addPbdConstraintFunctor(
        [&](PbdConstraintContainer& container)
        {
            VecDataArray<double, 3>& vertices1 =
                *std::dynamic_pointer_cast<PointSet>(clothObj->getPhysicsGeometry())->getVertexPositions();
            VecDataArray<double, 3>& vertices2 =
                *std::dynamic_pointer_cast<PointSet>(lineObj->getPhysicsGeometry())->getVertexPositions();

            auto constraint = std::make_shared<PbdDistanceConstraint>();
            constraint->initConstraint(0.0,
                { 0, 8 }, { 1, 2 }, 1000.0);
            container.addConstraint(constraint);
        });

    // Run the simulation
    {
        // Setup a viewer to render
        imstkNew<VTKViewer> viewer;
        viewer->setActiveScene(scene);

        // Setup a scene manager to advance the scene
        imstkNew<SceneManager> sceneManager;
        sceneManager->setActiveScene(scene);
        sceneManager->pause(); // Start simulation paused

        imstkNew<SimulationManager> driver;
        driver->addModule(viewer);
        driver->addModule(sceneManager);
        driver->setDesiredDt(0.001);

        // Add mouse and keyboard controls to the viewer
        {
            auto mouseControl = std::make_shared<MouseSceneControl>();
            mouseControl->setDevice(viewer->getMouseDevice());
            mouseControl->setSceneManager(sceneManager);
            viewer->addControl(mouseControl);

            auto keyControl = std::make_shared<KeyboardSceneControl>();
            keyControl->setDevice(viewer->getKeyboardDevice());
            keyControl->setSceneManager(sceneManager);
            keyControl->setModuleDriver(driver);
            viewer->addControl(keyControl);
        }

        using Vec3uc = Eigen::Matrix<unsigned char, 3, 1>;
        queueConnect<KeyEvent>(viewer->getKeyboardDevice(), &KeyboardDeviceClient::keyPress, sceneManager, [&](KeyEvent* e)
            {
                // Set new textures
                if (e->m_key == '1')
                {
                    setFleshTextures(clothObj->getVisualModel(0)->getRenderMaterial());
                }
                else if (e->m_key == '2')
                {
                    setFabricTextures(clothObj->getVisualModel(0)->getRenderMaterial());
                }
                // Darken the texture pixel values
                else if (e->m_key == 'h')
                {
                    auto imageData = clothObj->getVisualModel(0)->getRenderMaterial()->getTexture(Texture::Type::Diffuse)->getImageData();
                    std::shared_ptr<VecDataArray<unsigned char, 3>> scalars = std::dynamic_pointer_cast<VecDataArray<unsigned char, 3>>(imageData->getScalars());
                    Vec3uc* scalarPtr = scalars->getPointer();
                    for (int i = 0; i < scalars->size(); i++)
                    {
                        scalarPtr[i] = (scalarPtr[i].cast<double>() * 0.8).cast<unsigned char>();
                    }
                    clothObj->getVisualModel(0)->getRenderMaterial()->getTexture(Texture::Type::Diffuse)->postModified();
                }
        });

        driver->start();
    }

    return 0;
}