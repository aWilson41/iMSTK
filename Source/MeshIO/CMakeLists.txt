#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------

set(Dependencies
  Common
  Geometry    
  Assimp::Assimp
  VTK::IOGeometry
  VTK::IOImage
  VTK::IOLegacy
  VTK::IOPLY
  VTK::IOXML
  )

set(H_FILES
  imstkAssimpMeshIO.h
  imstkMeshIO.h
  imstkMshMeshIO.h
  imstkVTKMeshIO.h
  )
set(SRC_FILES
  imstkAssimpMeshIO.cpp
  imstkMeshIO.cpp
  imstkMshMeshIO.cpp
  imstkVTKMeshIO.cpp
  )

if (iMSTK_USE_VegaFEM)
  list(APPEND Dependencies VegaFEM::volumetricMesh)
  list(APPEND H_FILES imstkVegaMeshIO.h)
  list(APPEND SRC_FILES imstkVegaMeshIO.cpp)
endif()

include(imstkAddLibrary)
imstk_add_library(MeshIO
  H_FILES
    ${H_FILES}
  CPP_FILES
    ${SRC_FILES}
  DEPENDS
    ${Dependencies}
  )

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
if( ${PROJECT_NAME}_BUILD_TESTING )
  add_subdirectory("Testing")
endif()
