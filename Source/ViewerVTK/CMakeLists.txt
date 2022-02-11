#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------

set(H_FILES
  imstkAbstractVTKViewer.h
  imstkScreenCaptureUtility.h
  imstkVTKInteractorStyle.h
  imstkVTKInteractorStyleVR.h
  imstkVTKOpenVRViewer.h
  imstkVTKScreenCaptureUtility.h
  imstkVTKTextStatusManager.h
  imstkVTKViewer.h)

set(SRC_FILES
  imstkAbstractVTKViewer.cpp
  imstkScreenCaptureUtility.cpp
  imstkVTKInteractorStyle.cpp
  imstkVTKInteractorStyleVR.cpp
  imstkVTKOpenVRViewer.cpp
  imstkVTKScreenCaptureUtility.cpp
  imstkVTKTextStatusManager.cpp
  imstkVTKViewer.cpp)

if (NOT iMSTK_VTK_REPO_SOURCE STREQUAL "9.1.202202")
  list(APPEND H_FILES imstkVtkOpenVRRenderWindowInteractorImstk.h)
  list(APPEND SRC_FILES imstkVtkOpenVRRenderWindowInteractorImstk.cpp)
endif()

include(imstkAddLibrary)
imstk_add_library(ViewerVTK
  H_FILES
    ${H_FILES}
  CPP_FILES
    ${SRC_FILES}
  DEPENDS
    ViewerCore
    RenderingVTK
  )

# Install VR actions json
if (iMSTK_VTK_REPO_SOURCE STREQUAL "9.1.202202")
  file(GLOB VTK_VR_JSON_FILES vtk_openxr_*.json)
else()
  file(GLOB VTK_VR_JSON_FILES vtk_openvr_*.json)
endif()

if (MSVC)
  # Show them in visual studios
  target_sources(ViewerVTK PRIVATE ${VTK_VR_JSON_FILES})
  set_source_files_properties(${VTK_VR_JSON_FILES} PROPERTIES HEADER_FILE_ONLY TRUE)
endif()
install(FILES ${VTK_VR_JSON_FILES} DESTINATION bin)

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
# if( ${PROJECT_NAME}_BUILD_TESTING )
#   add_subdirectory(Testing)
# endif()