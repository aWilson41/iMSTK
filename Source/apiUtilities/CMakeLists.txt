#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------
include(imstkAddLibrary)      
imstk_add_library( imstkApiUtilities
  TYPE STATIC
  DEPENDS
    imstkCommon
    imstkSimulationManager
    imstkScene
  )
  
#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
#if( ${PROJECT_NAME}_BUILD_TESTING )
#  add_subdirectory( Testing )
#endif()