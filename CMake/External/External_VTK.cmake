#-----------------------------------------------------------------------------
# Dependencies
#-----------------------------------------------------------------------------
set(VTK_DEPENDENCIES "")
if(${${PROJECT_NAME}_ENABLE_VR})
  list(APPEND VTK_DEPENDENCIES "openvr")
endif(${${PROJECT_NAME}_ENABLE_VR})

#-----------------------------------------------------------------------------
# Add External Project
#-----------------------------------------------------------------------------
include(imstkAddExternalProject)
imstk_add_external_project( VTK
  GIT_REPOSITORY https://gitlab.kitware.com/vtk/vtk.git
  GIT_TAG 9b6a039f43404053a0653f742148d123f6ada7d6
  CMAKE_ARGS
      -DBUILD_EXAMPLES:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
      -DVTK_Group_StandAlone:BOOL=OFF
      -DVTK_Group_Rendering:BOOL=OFF
      -DModule_vtkRenderingOpenGL2:BOOL=ON
      -DModule_vtkIOXML:BOOL=ON
      -DModule_vtkIOLegacy:BOOL=ON
      -DModule_vtkIOPLY:BOOL=ON
      -DModule_vtkIOGeometry:BOOL=ON
      -DModule_vtkInteractionStyle:BOOL=ON
      -DModule_vtkRenderingAnnotation:BOOL=ON
      -DModule_vtkRenderingOpenVR:BOOL=${${PROJECT_NAME}_ENABLE_VR}
      -DModule_vtkInteractionWidgets:BOOL=ON
      -DModule_vtkglew:BOOL=ON
      -DModule_vtkRenderingContext2D:BOOL=ON
      -DVTK_RENDERING_BACKEND:STRING=OpenGL2
      -DVTK_WRAP_PYTHON:BOOL=OFF
      -DVTK_OPENVR_OBJECT_FACTORY:BOOL=OFF
      -DVTK_LEGACY_REMOVE:BOOL=ON
      -DVTK_BUILD_TESTING:STRING=OFF
  DEPENDENCIES ${VTK_DEPENDENCIES}
  RELATIVE_INCLUDE_PATH ""
  #VERBOSE
  )
