#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------

include(imstkAddLibrary)
imstk_add_library(DevicesOpenHaptics
  H_FILES
    imstkHapticDeviceClient.h
    imstkHapticDeviceManager.h
  CPP_FILES
    imstkHapticDeviceClient.cpp
    imstkHapticDeviceManager.cpp
  DEPENDS
    Devices
    OpenHapticsSDK::OpenHapticsSDK)

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
# if( ${PROJECT_NAME}_BUILD_TESTING )
#   add_subdirectory(Testing)
# endif()