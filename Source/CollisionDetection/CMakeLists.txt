#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------
include(imstkAddLibrary)
imstk_add_library( imstkCollisionDetection
  TYPE STATIC
  DEPENDS
    imstkCommon
    imstkDataStructures
    imstkGeometry
    imstkMeshIO #TODO: remove this dependency
    SCCD
  )

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
if( ${PROJECT_NAME}_BUILD_TESTING )

  include(imstkAddTest)

  imstk_add_test(imstkCollisionDetection)
  imstk_add_data(imstkCollisionDetection ${FILE_LIST_COL_TEST})
endif()
