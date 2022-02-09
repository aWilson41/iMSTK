include(imstkFind)
#-----------------------------------------------------------------------------
# Find All Headers and Libraries
#-----------------------------------------------------------------------------

#message("Using imstk FindOpenXR cmake")
imstk_find_header(OpenXR openxr.h)
imstk_find_libary(OpenXR openxr_loader " ")#Use same library for debug
imstk_find_package(OpenXR OpenXR::OpenXR)

#message(STATUS "OpenVR include : ${OpenXR_INCLUDE_DIRS}")
#message(STATUS "OpenVR libraries : ${OpenXR_LIBRARIES}")
