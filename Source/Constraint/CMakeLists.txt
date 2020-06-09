#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------
include(imstkAddLibrary)
imstk_add_library( imstkConstraints
  TYPE STATIC
  SUBDIR_LIST
    PbdConstraints
  DEPENDS
    imstkCommon
    imstkDataStructures
    imstkGeometry
  )

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
#if( ${PROJECT_NAME}_BUILD_TESTING )
#  add_subdirectory( Testing )
#endif()
