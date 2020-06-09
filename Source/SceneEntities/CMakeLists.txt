#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------
include(imstkAddLibrary)
imstk_add_library( imstkSceneEntities
  TYPE STATIC
  DEPENDS
    imstkCommon
    imstkGeometry
    imstkGeometryMappers
    imstkMaterials
    imstkMeshIO
    imstkConstraints
    imstkDevices
    imstkDynamicalModels
  )

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
#if( ${PROJECT_NAME}_BUILD_TESTING )
#  add_subdirectory( Testing )
#endif()