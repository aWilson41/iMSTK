#-----------------------------------------------------------------------------
# Dependencies
#-----------------------------------------------------------------------------
if (iMSTK_USE_RENDERING)
  set(VTK_DEPENDENCIES "OpenVR")
endif()

set(${PROJECT_NAME}_VTK_REPO_SOURCE "9.1" CACHE STRING "Select VTK Source Branch/Tag")
set(VTK_SOURCES "9.1;master;release;nightly-master" CACHE INTERNAL "List of available VTK branch,tags to get")
set_property(CACHE ${PROJECT_NAME}_VTK_REPO_SOURCE PROPERTY STRINGS ${VTK_SOURCES})

set(VTK_MODULE_SETTINGS
  -DVTK_MODULE_ENABLE_VTK_ChartsCore:STRING=YES
  -DVTK_MODULE_ENABLE_VTK_FiltersCore:STRING=YES
  -DVTK_MODULE_ENABLE_VTK_ImagingGeneral:STRING=YES
  -DVTK_MODULE_ENABLE_VTK_ImagingStencil:STRING=YES
  -DVTK_MODULE_ENABLE_VTK_IOExport:STRING=YES
  -DVTK_MODULE_ENABLE_VTK_IOGeometry:STRING=YES
  -DVTK_MODULE_ENABLE_VTK_IOImport:STRING=YES
  -DVTK_MODULE_ENABLE_VTK_IOLegacy:STRING=YES
  -DVTK_MODULE_ENABLE_VTK_IOPLY:STRING=YES
  -DVTK_MODULE_ENABLE_VTK_IOParallel:STRING=YES
  -DVTK_MODULE_ENABLE_VTK_IOParallelXML:STRING=YES
  -DVTK_MODULE_ENABLE_VTK_IOXML:STRING=YES
  -DVTK_MODULE_ENABLE_VTK_InteractionStyle:STRING=YES
  -DVTK_MODULE_ENABLE_VTK_InteractionWidgets:STRING=YES
  -DVTK_ENABLE_WRAPPING:STRING=OFF
  -DVTK_BUILD_EXAMPLES:STRING=DONT_WANT
  -DVTK_BUILD_TESTING:STRING=OFF
  -DVTK_GROUP_ENABLE_StandAlone:STRING=DONT_WANT
  -DVTK_GROUP_ENABLE_Rendering:STRING=DONT_WANT)

# If using VTK rendering include rendering submodules
if (iMSTK_MODULE_ENABLE_RenderingVTK)
  list(APPEND VTK_MODULE_SETTINGS
    -DVTK_MODULE_ENABLE_VTK_glew:STRING=YES
    -DVTK_MODULE_ENABLE_VTK_RenderingAnnotation:STRING=YES
    -DVTK_MODULE_ENABLE_VTK_RenderingContext2D:STRING=YES
    -DVTK_MODULE_ENABLE_VTK_RenderingContextOpenGL2:STRING=YES
    -DVTK_MODULE_ENABLE_VTK_RenderingExternal:STRING=YES
    -DVTK_MODULE_ENABLE_VTK_RenderingOpenGL2:STRING=YES
    -DVTK_MODULE_ENABLE_VTK_RenderingOpenVR:STRING=YES
    -DVTK_MODULE_ENABLE_VTK_RenderingVolumeOpenGL2:STRING=YES
    -DVTK_MODULE_ENABLE_VTK_ViewsContext2D:STRING=YES)
endif()
if (iMSTK_MODULE_ENABLE_Filtering)
  list(APPEND VTK_MODULE_SETTINGS
    -DVTK_MODULE_ENABLE_VTK_ImagingMath:STRING=YES
    -DVTK_MODULE_ENABLE_VTK_FiltersExtraction:STRING=YES
    -DVTK_MODULE_ENABLE_VTK_FiltersGeneral:STRING=YES
    -DVTK_MODULE_ENABLE_VTK_FiltersModeling:STRING=YES)
endif()

# If using VTK offscreen rendering include these options
if (${PROJECT_NAME}_USE_VTK_OSMESA)
  list(APPEND VTK_MODULE_SETTINGS
    -DVTK_DEFAULT_RENDER_WINDOW_OFFSCREEN=1
    -DVTK_OPENGL_HAS_OSMESA=1
    -DVTK_USE_COCOA=0
    -DVTK_USE_X=0)
endif()

set(${PROJECT_NAME}_VTK_SOURCE GIT_REPOSITORY https://gitlab.kitware.com/vtk/vtk.git)

if(${PROJECT_NAME}_VTK_REPO_SOURCE EQUAL "9.1")
  set(${PROJECT_NAME}_VTK_HASH GIT_TAG 285daeedd58eb890cb90d6e907d822eea3d2d092)
else()
  set(${PROJECT_NAME}_VTK_HASH GIT_TAG origin/${${PROJECT_NAME}_VTK_REPO_SOURCE})
endif()

#-----------------------------------------------------------------------------
# Add External Project
#-----------------------------------------------------------------------------
include(imstkAddExternalProject)
imstk_add_external_project( VTK
  ${${PROJECT_NAME}_VTK_SOURCE}
  ${${PROJECT_NAME}_VTK_HASH}
  CMAKE_ARGS
       ${VTK_MODULE_SETTINGS}
      -DVTK_WRAP_PYTHON:BOOL=OFF
      -DVTK_LEGACY_REMOVE:BOOL=ON
      -DCMAKE_INSTALL_RPATH:PATH=$ORIGIN/../lib
  DEPENDENCIES ${VTK_DEPENDENCIES}
  RELATIVE_INCLUDE_PATH ""
  #VERBOSE
)
if(NOT USE_SYSTEM_VTK)
  set(VTK_DIR ${CMAKE_INSTALL_PREFIX}/lib/cmake/vtk-9.1)
  #message(STATUS "VTK_DIR : ${VTK_DIR}")
endif()
