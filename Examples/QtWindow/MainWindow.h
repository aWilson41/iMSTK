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

#include <qmainwindow.h>
#include <qpointer.h>
#include <QVTKOpenGLNativeWidget.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>

#include <vtkSphereSource.h>

class MainWindow : public QMainWindow
{
Q_OBJECT
public:
    MainWindow(QWidget* parent = Q_NULLPTR) : QMainWindow(parent)
    {
        ui.setupUi(this);

        // Create a VTK widget
        ui.qvtkWidget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        vtkSmartPointer<vtkGenericOpenGLRenderWindow> renWin = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
        ui.qvtkWidget->setRenderWindow(renWin);

        /* vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
         ren->SetBackground(0.1, 0.1, 1.0);
         ui.qvtkWidget->renderWindow()->AddRenderer(ren);*/
    }

    virtual ~MainWindow() = default;

public:
    QVTKOpenGLNativeWidget* getQvtkWidget() { return ui.qvtkWidget; }

private:
    Ui::MainWindowWidget ui;
};