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

#include "ui_MainWindow.h"
#include "imstkSimulationManager.h"
#include "imstkSceneManager.h"

#include <qmainwindow.h>
#include <QVTKOpenGLNativeWidget.h>
#include <vtkGenericOpenGLRenderWindow.h>

class MainWindow : public QMainWindow
{
Q_OBJECT
public:
    MainWindow(QWidget* parent = Q_NULLPTR) : QMainWindow(parent)
    {
        ui.setupUi(this);

        // Create a VTK widget
        ui.qvtkWidget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        renWin = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
        ui.qvtkWidget->setRenderWindow(renWin);

        connect(ui.startButton, &QPushButton::pressed, this, &MainWindow::toggleSimulation);
    }

    virtual ~MainWindow() = default;

public:
    void closeEvent(QCloseEvent* e) override
    {
        // This will only stop it on the next iteration
        m_driver->requestStatus(ModuleDriverStopped);
        // But it should be stopped before window closed here
        QMainWindow::closeEvent(e);
    }

    void toggleSimulation()
    {
        if (m_sceneManager->getPaused())
        {
            m_sceneManager->setPaused(false);
            ui.startButton->setText("Pause Simulation");
        }
        else
        {
            m_sceneManager->setPaused(true);
            ui.startButton->setText("Start Simulation");
        }
    }

    void setSimManager(std::shared_ptr<imstk::SimulationManager> driver)
    {
        m_driver = driver;

        for (auto module : m_driver->getModules())
        {
            if (m_sceneManager = std::dynamic_pointer_cast<imstk::SceneManager>(module))
            {
                break;
            }
        }
    }

    QVTKOpenGLNativeWidget* getQvtkWidget() { return ui.qvtkWidget; }

private:
    Ui::MainWindowWidget ui;
    std::shared_ptr<imstk::SimulationManager>     m_driver       = nullptr;
    std::shared_ptr<imstk::SceneManager>          m_sceneManager = nullptr;
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renWin;
};