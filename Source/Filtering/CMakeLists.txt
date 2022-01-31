#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------
file(GLOB hFiles *.h)
file(GLOB cppFiles *cpp)

include(imstkAddLibrary)
imstk_add_library( Filtering
  H_FILES
    ${hFiles}
  CPP_FILES
    ${cppFiles}
  DEPENDS
    FilteringCore
    VTK::ImagingMath
    VTK::FiltersExtraction
    VTK::FiltersGeneral
    VTK::FiltersModeling)

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
if( ${PROJECT_NAME}_BUILD_TESTING )
  add_subdirectory(Testing)
endif()