#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------
include(imstkAddLibrary)
imstk_add_library( MeshIO
  DEPENDS
    Common
	Geometry    
    Assimp::Assimp
    VegaFEM::volumetricMesh
    VTK::IOGeometry
    VTK::IOImage
    VTK::IOLegacy
    VTK::IOPLY
    VTK::IOXML)

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
if( ${PROJECT_NAME}_BUILD_TESTING )
    add_subdirectory("Testing")
endif()
