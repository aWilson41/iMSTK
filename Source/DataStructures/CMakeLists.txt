#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------
include(imstkAddLibrary)
imstk_add_library( imstkDataStructures
  TYPE STATIC
  DEPENDS
    imstkCommon
    imstkGeometry #TODO remove this dependency
  )

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
if( ${PROJECT_NAME}_BUILD_TESTING )
  include(imstkAddTest)
  imstk_add_test( imstkDataStructures )
endif()
