#-----------------------------------------------------------------------------
# Add External Project
#-----------------------------------------------------------------------------
include(imstkAddExternalProject)

# Download options
if(NOT DEFINED iMSTK_Assimp_GIT_SHA)
  set(iMSTK_Assimp_GIT_SHA "bda31ee2315f372c4f0de7c29229f5c9d075fc77") # imstk-v5.2.5-2022-10-16-9519a62dd
endif()
if(NOT DEFINED iMSTK_Assimp_GIT_REPOSITORY)
  set(iMSTK_Assimp_GIT_REPOSITORY "https://gitlab.kitware.com/iMSTK/assimp.git")
endif()

set(ASSIMP_MODULE_SETTINGS
  -DASSIMP_BUILD_ZLIB:BOOL=ON
  -DASSIMP_BUILD_ASSIMP_TOOLS:BOOL=OFF
  -DASSIMP_BUILD_TESTS:BOOL=OFF
  -DASSIMP_NO_EXPORT:BOOL=OFF
  )
if (${PROJECT_NAME}_BUILD_FOR_ANDROID)
  list(APPEND ASSIMP_MODULE_SETTINGS
    -DANDROID_NDK:STRING=${ANDROID_NDK}
    -DANDROID_ABI:STRING=${ANDROID_ARCH_ABI}
    -DCMAKE_TOOLCHAIN_FILE:PATH=${CMAKE_BINARY_DIR}/${_ANDROID_DIR}-toolchain.cmake
    -DCMAKE_MAKE_PROGRAM:FILEPATH=${CMAKE_MAKE_PROGRAM}
    )
endif()

imstk_add_external_project( Assimp
  GIT_REPOSITORY ${iMSTK_Assimp_GIT_REPOSITORY}
  GIT_TAG ${iMSTK_Assimp_GIT_SHA}
  CMAKE_CACHE_ARGS
    -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
    -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON
    ${ASSIMP_MODULE_SETTINGS}
    -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
    -DLIBRARY_SUFFIX:STRING=
  RELATIVE_INCLUDE_PATH "include"
  #DEPENDENCIES ""
  #VERBOSE
)
