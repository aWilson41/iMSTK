#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------

include(imstkAddLibrary)
imstk_add_library(Devices
  H_FILES
    imstkDeviceClient.h
    imstkDummyClient.h
    imstkKeyboardDeviceClient.h
    imstkMouseDeviceClient.h
    imstkOpenVRDeviceClient.h
  CPP_FILES
    imstkDeviceClient.cpp
    imstkDummyClient.cpp
    imstkKeyboardDeviceClient.cpp
    imstkMouseDeviceClient.cpp
    imstkOpenVRDeviceClient.cpp
  DEPENDS
    Common)

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
#if( ${PROJECT_NAME}_BUILD_TESTING )
#  add_subdirectory( Testing )
#endif()