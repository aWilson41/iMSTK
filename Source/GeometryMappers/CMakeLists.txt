#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------

include(imstkAddLibrary)
imstk_add_library( GeometryMappers
  H_FILES
    imstkGeometryMap.h
	imstkIsometricMap.h
	imstkPointToCellMap.h
	imstkPointToTetMap.h
	imstkPointwiseMap.h
	imstkTriangleToTetMap.h
  CPP_FILES
    imstkGeometryMap.cpp
	imstkIsometricMap.cpp
	imstkPointToCellMap.cpp
	imstkPointToTetMap.cpp
	imstkPointwiseMap.cpp
	imstkTriangleToTetMap.cpp
  DEPENDS
	FilteringCore
	CollisionDetection
  )

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
if( ${PROJECT_NAME}_BUILD_TESTING )
  add_subdirectory( Testing )
endif()
