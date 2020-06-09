#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------
include(imstkAddLibrary)
imstk_add_library( imstkScene
  TYPE STATIC
  DEPENDS
    imstkCommon
    imstkCollisionDetection
    imstkCollisionHandling
    imstkSceneEntities
    imstkControllers
  )

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
#if( ${PROJECT_NAME}_BUILD_TESTING )
#  add_subdirectory( Testing )
#endif()