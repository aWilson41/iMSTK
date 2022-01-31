#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------
include(imstkAddLibrary)
imstk_add_library( Geometry
  DEPENDS
    Common
    VTK::CommonDataModel
    VTK::FiltersSources
    VTK::FiltersGeneral)

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
if( ${PROJECT_NAME}_BUILD_TESTING )
  add_subdirectory(Testing)
endif()
