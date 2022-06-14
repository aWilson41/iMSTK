#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------

include(imstkAddLibrary)
imstk_add_library(ViewerVTK
  H_FILES
    imstkAbstractVTKViewer.h
    imstkScreenCaptureUtility.h
    imstkVTKInteractorStyle.h
    imstkVTKScreenCaptureUtility.h
    imstkVTKTextStatusManager.h
    imstkVTKViewer.h
  CPP_FILES
    imstkAbstractVTKViewer.cpp
    imstkScreenCaptureUtility.cpp
    imstkVTKInteractorStyle.cpp
    imstkVTKScreenCaptureUtility.cpp
    imstkVTKTextStatusManager.cpp
    imstkVTKViewer.cpp
  DEPENDS
    ViewerCore
    RenderingVTK
    VTK::InteractionStyle
  )

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
# if( ${PROJECT_NAME}_BUILD_TESTING )
#   add_subdirectory(Testing)
# endif()