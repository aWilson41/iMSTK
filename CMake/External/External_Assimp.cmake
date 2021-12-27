#-----------------------------------------------------------------------------
# Add External Project
#-----------------------------------------------------------------------------
include(imstkAddExternalProject)

# Download options
if(NOT DEFINED iMSTK_Assimp_GIT_SHA)
  set(iMSTK_Assimp_GIT_SHA "fbcfce0cd1800eec6944325b21bbd9006ed2e9d1")
endif()
if(NOT DEFINED iMSTK_Assimp_GIT_REPOSITORY)
  set(EXTERNAL_PROJECT_DOWNLOAD_OPTIONS
    URL https://gitlab.kitware.com/iMSTK/assimp/-/archive/${iMSTK_Assimp_GIT_SHA}/assimp-${iMSTK_Assimp_GIT_SHA}.zip
    URL_HASH MD5=ed3be5abbba1b9e5224d6247fa84d6db
    )
else()
  set(EXTERNAL_PROJECT_DOWNLOAD_OPTIONS
    GIT_REPOSITORY ${iMSTK_Assimp_GIT_REPOSITORY}
    GIT_TAG ${iMSTK_Assimp_GIT_SHA}
    )
endif()

imstk_add_external_project( Assimp
  ${EXTERNAL_PROJECT_DOWNLOAD_OPTIONS}
  CMAKE_CACHE_ARGS
    #-DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON
    -DASSIMP_BUILD_ASSIMP_TOOLS:BOOL=OFF
    -DASSIMP_BUILD_TESTS:BOOL=OFF
    #-DASSIMP_NO_EXPORT:BOOL=ON
    -DLIBRARY_SUFFIX:STRING=
    -DASSIMP_BUILD_ZLIB:BOOL=ON
    #-DBUILD_SHARED_LIBS:BOOL=OFF
  RELATIVE_INCLUDE_PATH "include"
  #DEPENDENCIES ""
  #VERBOSE
)
