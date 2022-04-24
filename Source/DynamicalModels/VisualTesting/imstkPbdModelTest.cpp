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
#include "imstkDirectionalLight.h"
#include "imstkPbdModel.h"
#include "imstkPbdObject.h"
#include "imstkScene.h"
#include "imstkSceneManager.h"
#include "imstkSimulationManager.h"
#include "imstkSurfaceMesh.h"
#include "imstkTetrahedralMesh.h"
#include "imstkVisualTestingUtils.h"

using namespace imstk;

///
/// \class PbdModelVisualTest
///
class PbdModelVisualTest : public VisualTest
{
public:
    void SetUp() override
    {
        VisualTest::SetUp();

        m_pbdObj    = std::make_shared<PbdObject>("pbdobj");
        m_pbdModel  = std::make_shared<PbdModel>();
        m_pbdConfig = std::make_shared<PbdModelConfig>();
        m_pbdModel->configure(m_pbdConfig);
        m_pbdObj->setDynamicalModel(m_pbdModel);
    }

    ///
    /// \brief Create a scene composed of two collision objects with
    /// the respective collision geometries and method
    ///
    void createScene()
    {
        // Setup the scene
        m_scene = std::make_shared<Scene>(::testing::UnitTest::GetInstance()->current_test_info()->name());
        if (m_camera != nullptr)
        {
            *m_scene->getActiveCamera() = *m_camera;
        }
        m_scene->addSceneObject(m_pbdObj);

        // Light
        auto light = std::make_shared<DirectionalLight>();
        light->setFocalPoint(Vec3d(5.0, -8.0, -5.0));
        light->setIntensity(1.0);
        m_scene->addLight("Light", light);
    }

public:
    std::shared_ptr<PbdObject>      m_pbdObj    = nullptr;
    std::shared_ptr<PbdModel>       m_pbdModel  = nullptr;
    std::shared_ptr<PbdModelConfig> m_pbdConfig = nullptr;

    std::shared_ptr<imstk::Camera> m_camera = nullptr;
};

///
/// \brief
///
TEST_F(PbdModelVisualTest, NeoHookeanTetInversion)
{
    m_camera = std::make_shared<Camera>();
    m_camera->setFocalPoint(-0.0366287, 0.420204, 0.474284);
    m_camera->setPosition(-2.60143, 1.23713, 2.42823);
    m_camera->setViewUp(0.216266, 0.968787, -0.121162);

    VecDataArray<double, 3> vertices(4);
    vertices[0] = Vec3d(0.5, 0.0, -1.0 / 3.0);
    vertices[1] = Vec3d(-0.5, 0.0, -1.0 / 3.0);
    vertices[2] = Vec3d(0.0, 0.0, 2.0 / 3.0);
    vertices[3] = Vec3d(0.0, 1.0, 0.0);

    VecDataArray<int, 4> indices(1);
    indices[0] = Vec4i(0, 1, 2, 3);

    auto tetMesh = std::make_shared<TetrahedralMesh>();
    tetMesh->initialize(
        std::make_shared<VecDataArray<double, 3>>(vertices),
        std::make_shared<VecDataArray<int, 4>>(indices));

    m_pbdConfig->enableFemConstraint(PbdFemConstraint::MaterialType::NeoHookean);
    m_pbdConfig->m_gravity = Vec3d::Zero();
    m_pbdConfig->m_dt      = 0.01;
    m_pbdConfig->m_doPartitioning            = false;
    m_pbdConfig->m_viscousDampingCoeff       = 0.0;
    m_pbdConfig->m_femParams->m_YoungModulus = 1.0;
    m_pbdConfig->m_femParams->m_PoissonRatio = 0.45;
    m_pbdConfig->m_iterations = 1;
    m_pbdModel->setModelGeometry(tetMesh);
    m_pbdObj->setDynamicalModel(m_pbdModel);
    m_pbdObj->setPhysicsGeometry(tetMesh);
    m_pbdObj->setVisualGeometry(tetMesh);
    m_pbdObj->getPbdBody()->fixedNodeIds     = { 3 };
    m_pbdObj->getPbdBody()->uniformMassValue = 1.0;

    int i = 0;
    connect<Event>(m_sceneManager, &SceneManager::postUpdate,
        [&](Event*)
        {
            if (i == 0)
            {
                (*tetMesh->getVertexPositions())[3] = Vec3d(0.0, -0.01, 0.0);
                i++;
            }

            /*const double dt = m_sceneManager->getDt();
            (*tetMesh->getVertexPositions())[3] += Vec3d(0.0, -dt * 0.01, 0.0);*/
        });

    m_sceneManager->setPaused(true);
    createScene();
    runFor(2.0);
}