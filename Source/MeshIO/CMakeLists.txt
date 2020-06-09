#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------
include(imstkAddLibrary)
imstk_add_library( imstkMeshIO
  TYPE STATIC
  DEPENDS
    imstkCommon
	  imstkGeometry
    glm
    Assimp    
    VegaFEM::volumetricMesh
    ${VTK_LIBRARIES}
  )

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
#if( ${PROJECT_NAME}_BUILD_TESTING )
#  include(imstkAddTest)
#  imstk_add_test( Geometry )
#endif()
