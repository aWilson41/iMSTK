#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------

include(imstkAddLibrary)
imstk_add_library(DevicesVRPN
  H_FILES
    imstkVRPNDeviceClient.h
    imstkVRPNDeviceManager.h
  CPP_FILES
    imstkVRPNDeviceClient.cpp
    imstkVRPNDeviceManager.cpp
  DEPENDS
    Devices
    VRPN::VRPN)

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
# if( ${PROJECT_NAME}_BUILD_TESTING )
#   add_subdirectory(Testing)
# endif()