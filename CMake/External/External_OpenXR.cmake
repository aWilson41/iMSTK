#-----------------------------------------------------------------------------
# Add External Project
#-----------------------------------------------------------------------------
include(imstkAddExternalProject)

# Download options
if(NOT DEFINED iMSTK_OpenXR_GIT_SHA)
  set(iMSTK_OpenXR_GIT_SHA "7450c5420789ff0ab63aaad2c4ed018807629a77")
endif()
if(NOT DEFINED iMSTK_OpenXR_GIT_REPOSITORY)
  set(EXTERNAL_PROJECT_DOWNLOAD_OPTIONS
    URL https://github.com/KhronosGroup/OpenXR-SDK-Source/archive/${iMSTK_OpenXR_GIT_SHA}.zip
    URL_HASH MD5=b77fb4ad8b8d8fe9d4174426361a9faa
    )
else()
  set(EXTERNAL_PROJECT_DOWNLOAD_OPTIONS
    GIT_REPOSITORY ${iMSTK_OpenXR_GIT_REPOSITORY}
    GIT_TAG ${iMSTK_OpenXR_GIT_SHA}
    )
endif()

imstk_add_external_project( OpenXR
  ${EXTERNAL_PROJECT_DOWNLOAD_OPTIONS}
  CMAKE_CACHE_ARGS
    -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON
    -DDYNAMIC_LOADER:BOOL=ON # Build dynamic dll instead of static
    -DBUILD_TESTS:BOOL=OFF)
