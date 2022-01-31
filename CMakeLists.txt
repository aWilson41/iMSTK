cmake_minimum_required(VERSION 3.9)
project(iMSTK VERSION 5.0.0 LANGUAGES C CXX)

#-----------------------------------------------------------------------------
# Set a default build type if none was specified
#-----------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

set(CMAKE_DEBUG_POSTFIX "d")

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

#-----------------------------------------------------------------------------
# Project install directories
#-----------------------------------------------------------------------------
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install CACHE PATH "Install location" FORCE)
endif()
set(CMAKE_PREFIX_PATH ${CMAKE_INSTALL_PREFIX})
# Let's go ahead and make these directories
file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/bin)
file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/include)
file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib) 

#-----------------------------------------------------------------------------
# Update CMake module path & CMake dir
#-----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/CMake
    ${CMAKE_CURRENT_SOURCE_DIR}/CMake/Utilities
    ${CMAKE_INSTALL_PREFIX}
    ${CMAKE_INSTALL_PREFIX}/lib/cmake # Vega and VTK
    )
set(${PROJECT_NAME}_CMAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CMake)

# Add "mark_as_superbuild" CMake function to associate option with inner build
include(imstkSolveDependencies)

#-----------------------------------------------------------------------------
# C++11 Support
#-----------------------------------------------------------------------------
set(_msg "Setting C++ standard")
message(STATUS "${_msg}")
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 14)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message(STATUS "${_msg} - C++${CMAKE_CXX_STANDARD}")
if(NOT CMAKE_CXX_STANDARD MATCHES "^(11|14|17)$")
  message(FATAL_ERROR "CMAKE_CXX_STANDARD must be set to 11, 14 or 17")
endif()

set(THREADS_PREFER_PTHREAD_FLAG ON)

# Prevents a compiler error for Visual Studio 15.8
if(MSVC)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  add_definitions(-D_SCL_SECURE_NO_WARNINGS)
  if(${MSVC_VERSION} GREATER_EQUAL 1915)
    add_definitions(-D_DISABLE_EXTENDED_ALIGNED_STORAGE)
  endif()
endif()

#-----------------------------------------------------------------------------
# Options
#-----------------------------------------------------------------------------

# CUDA
# option(${PROJECT_NAME}_ENABLE_CUDA_BACKEND "Enable iMSTK CUDA backend" OFF)
# mark_as_superbuild(${PROJECT_NAME}_ENABLE_CUDA_BACKEND)
# if (${PROJECT_NAME}_ENABLE_CUDA_BACKEND)
#   include(CheckLanguage)
#   check_language(CUDA)
#   if (CMAKE_CUDA_COMPILER)
#     enable_language(CUDA)
# 	set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
# 	set(CUDA_PROPAGATE_HOST_FLAGS ON)
# 	if(NOT DEFINED CMAKE_CUDA_STANDARD)
# 	  set(CMAKE_CUDA_STANDARD 11)
# 	  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
# 	endif()
#   else()
# 	message(STATUS "WARNING: CUDA compiler NOT FOUND; CUDA backend not enabled!")
#     set(${PROJECT_NAME}_ENABLE_CUDA_BACKEND OFF)
#   endif()
# endif ()

#-----------------------------------------------------------------------------
# Define all the modules
# Note: Can be acquired with iMSTK_MODULE_LIST
#-----------------------------------------------------------------------------
include(imstkModule)
set(iMSTK_MODULE_LIST "")
define_module(Common ON)
define_module(Materials ON)
define_module(Geometry ON)
define_module(MeshIO ON)
define_module(GeometryMappers ON)
define_module(DataStructures ON)
define_module(Constraint ON)
define_module(DevicesCore ON)
define_module(DevicesVRPN OFF)
define_module(DevicesOpenHaptics OFF)
define_module(Filtering ON)
define_module(FilteringCore ON)
define_module(RenderingCore ON)
define_module(RenderingVTK ON)
define_module(Solvers ON)
define_module(DynamicalModels ON)
define_module(SceneEntities ON)
define_module(Animation ON)
define_module(Controllers ON)
define_module(CollisionDetection ON)
define_module(CollisionHandling ON)
define_module(Scene ON)
define_module(SimulationManager ON)
define_module(ViewerCore ON)
define_module(ViewerVTK ON)
define_module(Wrappers OFF)

#-----------------------------------------------------------------------------
# Update module groups
# Note: These should not be used elsewhere. These indicate what modules the user
# wants, not what ones are present.
#-----------------------------------------------------------------------------
set(${PROJECT_NAME}_MODULE_GROUP_RENDERING "DEFAULT" CACHE STRING "Provides default rendering support (VTK)")
set_property(CACHE ${PROJECT_NAME}_MODULE_GROUP_RENDERING PROPERTY STRINGS YES DEFAULT NO)
mark_as_superbuild(${PROJECT_NAME}_MODULE_GROUP_RENDERING:STRING)
if (${PROJECT_NAME}_MODULE_GROUP_RENDERING STREQUAL "YES")
  set(${PROJECT_NAME}_MODULE_ENABLE_ViewerCore ON CACHE BOOL "" FORCE)
  set(${PROJECT_NAME}_MODULE_ENABLE_ViewerVTK ON CACHE BOOL "" FORCE)
  set(${PROJECT_NAME}_MODULE_ENABLE_RenderingCore ON CACHE BOOL "" FORCE)
  set(${PROJECT_NAME}_MODULE_ENABLE_RenderingVTK ON CACHE BOOL "" FORCE)
  set(${PROJECT_NAME}_MODULE_ENABLE_SimulationManager ON CACHE BOOL "" FORCE)
elseif(${PROJECT_NAME}_MODULE_GROUP_RENDERING STREQUAL "NO")
  set(${PROJECT_NAME}_MODULE_ENABLE_ViewerCore OFF CACHE BOOL "" FORCE)
  set(${PROJECT_NAME}_MODULE_ENABLE_ViewerVTK OFF CACHE BOOL "" FORCE)
  set(${PROJECT_NAME}_MODULE_ENABLE_RenderingCore OFF CACHE BOOL "" FORCE)
  set(${PROJECT_NAME}_MODULE_ENABLE_RenderingVTK OFF CACHE BOOL "" FORCE)
  set(${PROJECT_NAME}_MODULE_ENABLE_SimulationManager OFF CACHE BOOL "" FORCE)
endif()

# option(${PROJECT_NAME}_MODULE_GROUP_AUDIO "Provides default audio support (SFML)" OFF)
# mark_as_superbuild(${PROJECT_NAME}_MODULE_GROUP_AUDIO:BOOL)
# if (${PROJECT_NAME}_MODULE_GROUP_AUDIO)
#   set(${PROJECT_NAME}_MODULE_ENABLE_SFML ON CACHE BOOL "")
# endif()

set(${PROJECT_NAME}_MODULE_GROUP_HAPTICS "DEFAULT" CACHE STRING "Provides default haptic support (OpenHaptics)")
set_property(CACHE ${PROJECT_NAME}_MODULE_GROUP_HAPTICS PROPERTY STRINGS YES DEFAULT NO)
mark_as_superbuild(${PROJECT_NAME}_MODULE_GROUP_HAPTICS:STRING)
if (${PROJECT_NAME}_MODULE_GROUP_HAPTICS STREQUAL "YES")
  set(${PROJECT_NAME}_MODULE_ENABLE_DevicesCore ON CACHE BOOL "" FORCE)
  set(${PROJECT_NAME}_MODULE_ENABLE_DevicesOpenHaptics ON CACHE BOOL "" FORCE)
elseif(${PROJECT_NAME}_MODULE_GROUP_HAPTICS STREQUAL "NO")
  set(${PROJECT_NAME}_MODULE_ENABLE_DevicesCore OFF CACHE BOOL "" FORCE)
  set(${PROJECT_NAME}_MODULE_ENABLE_DevicesOpenHaptics OFF CACHE BOOL "" FORCE)
endif()

#-----------------------------------------------------------------------------
# Check module configuration, throw messages and correct invalid configurations
#-----------------------------------------------------------------------------
if (${PROJECT_NAME}_MODULE_ENABLE_DevicesOpenHaptics AND WIN32)
  message(FATAL_ERROR "Setting ${PROJECT_NAME}_MODULE_ENABLE_DevicesOpenHaptics to ON is only supported on Windows")
  set(${PROJECT_NAME}_MODULE_ENABLE_DevicesOpenHaptics OFF CACHE BOOL "" FORCE)
endif()

#-----------------------------------------------------------------------------
# Update module configurations
# Note: We changed some in states in the previous blocks so we should update our list of provided modules
#-----------------------------------------------------------------------------
update_module_list()
mark_as_superbuild(${iMSTK_MODULE_LIST}:STRING)

#-----------------------------------------------------------------------------
# Configure other iMSTK options
#-----------------------------------------------------------------------------
# Build with examples
option(${PROJECT_NAME}_BUILD_EXAMPLES "Build iMSTK examples" ON)
mark_as_superbuild(${PROJECT_NAME}_BUILD_EXAMPLES:BOOL)

# Build with testing
option(${PROJECT_NAME}_BUILD_TESTING "Build iMSTK tests" ON)
mark_as_superbuild(${PROJECT_NAME}_BUILD_TESTING:BOOL)

set(BUILD_TESTING OFF)
if (${PROJECT_NAME}_BUILD_TESTING)
  set(BUILD_TESTING ON)
endif ()
if (UNIX)
  option(${PROJECT_NAME}_USE_MODEL_REDUCTION "Build with model reduction, requires a VegaFEM built with Intel MKL and arpack" OFF)
  mark_as_superbuild(${PROJECT_NAME}_USE_MODEL_REDUCTION:BOOL)

  option(${PROJECT_NAME}_COLOR_OUTPUT "For coloring output in gcc" OFF)
  mark_as_superbuild(${PROJECT_NAME}_COLOR_OUTPUT:BOOL)

  option(${PROJECT_NAME}_COVERAGE "Build for calculating coverage" OFF)
  mark_as_superbuild(${PROJECT_NAME}_COVERAGE:BOOL)

  option(${PROJECT_NAME}_USE_VTK_OSMESA "Build with VTK OSMesa rendering" OFF)
  mark_as_superbuild(${PROJECT_NAME}_USE_VTK_OSMESA:BOOL)

  if (${PROJECT_NAME}_COVERAGE AND NOT ${CMAKE_BUILD_TYPE} MATCHES [Dd][Ee][Bb][Uu][Gg])
    message(WARNING "Coverage requires a debug build forcing debug")
    set(CMAKE_BUILD_TYPE "DEBUG")
  endif()  
endif()

# Build with CSharp wrapper
if (${PROJECT_NAME}_MODULE_ENABLE_Wrappers)
  find_package(SWIG REQUIRED)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)

  include_directories(${CMAKE_CURRENT_SOURCE_DIR})

  set(CMAKE_SWIG_FLAGS "")
endif()


#-----------------------------------------------------------------------------
# CTest/Dashboards
#-----------------------------------------------------------------------------
if (${PROJECT_NAME}_BUILD_TESTING)
  include(CTest)
  set_property(CACHE BUILD_TESTING PROPERTY TYPE INTERNAL)
endif ()


#-----------------------------------------------------------------------------
# SUPERBUILD
#-----------------------------------------------------------------------------
option(${PROJECT_NAME}_SUPERBUILD "Build ${PROJECT_NAME} and the projects it depends on." ON)

if(${PROJECT_NAME}_SUPERBUILD)

  #-----------------------------------------------------------------------------
  # Define External dependencies
  #-----------------------------------------------------------------------------
  macro(imstk_define_dependency extProj)
    list(APPEND ${PROJECT_NAME}_DEPENDENCIES ${extProj})
    option(USE_SYSTEM_${extProj} "Exclude ${extProj} from superbuild and use an existing build." OFF)
    mark_as_advanced(USE_SYSTEM_${extProj})
  endmacro()

  option(${PROJECT_NAME}_USE_Uncrustify "Use Uncrustify as a code style beautifier." ON)
  if(${PROJECT_NAME}_USE_Uncrustify)
    find_package (Python COMPONENTS Interpreter)
    if(Python_Interpreter_FOUND)
      imstk_define_dependency(Uncrustify)
    else()
      message(WARN "Uncrustify requires a python interpreter, please install python and rerun configure")
    endif()
  endif()
  
  imstk_define_dependency(Assimp)
  imstk_define_dependency(Eigen3)
  imstk_define_dependency(g3log)

  if (iMSTK_MODULE_ENABLE_DevicesVRPN)
    imstk_define_dependency(LibNiFalcon)
  endif()

  if (iMSTK_MODULE_ENABLE_RenderingVTK)
    imstk_define_dependency(OpenVR)
  endif()

  imstk_define_dependency(TBB)
  imstk_define_dependency(VegaFEM)
  imstk_define_dependency(VTK)

  if(iMSTK_MODULE_ENABLE_DevicesOpenHaptics)
    imstk_define_dependency(OpenHaptics)
  endif()

  if (iMSTK_MODULE_ENABLE_DevicesVRPN)
    imstk_define_dependency(VRPN)
    imstk_define_dependency(Libusb) #for VRPN
    imstk_define_dependency(FTD2XX) #for LibNiFalcon (which is for VRPN)
  endif()

  # if (iMSTK_MODULE_ENABLE_SFML)
  #   imstk_define_dependency(SFML)
  # endif()

  if(${PROJECT_NAME}_BUILD_TESTING)
    imstk_define_dependency(GTest)

    #-----------------------------------------------------------------------------
    # Allow CTest to cover Innerbuild
    #-----------------------------------------------------------------------------
    configure_file(
      "${CMAKE_CURRENT_LIST_DIR}/CMake/Utilities/imstkCTestAddInnerbuild.cmake.in"
      "${CMAKE_CURRENT_BINARY_DIR}/imstkCTestAddInnerbuild.cmake"
      @ONLY
    )
    set_directory_properties(PROPERTIES TEST_INCLUDE_FILE
      "${CMAKE_CURRENT_BINARY_DIR}/imstkCTestAddInnerbuild.cmake"
    )
  endif()

  if (${PROJECT_NAME}_BUILD_TESTING OR ${PROJECT_NAME}_BUILD_EXAMPLES)
    imstk_define_dependency(iMSTKData)
  endif()

  #-----------------------------------------------------------------------------
  # Solve project dependencies
  #-----------------------------------------------------------------------------
  # Call CMakeLists.txt in CMake/External which will solve the dependencies
  # and add the External projects, including this one: this top-level
  # CMakeLists.txt will be called back with SUPERBUILD=OFF, to execute
  # the rest of the code below (INNERBUILD), which explains the `return`
  add_subdirectory(CMake/External)

  return()

endif()

#-----------------------------------------------------------------------------
#                               INNERBUILD
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Find external dependencies
#-----------------------------------------------------------------------------
include(imstkFind)
# Set where to look for packages (If not using system)
set(PACKAGE_PREFIX_DIR ${CMAKE_INSTALL_PREFIX})

# Assimp
find_package( Assimp REQUIRED )

# Eigen
find_package( Eigen3 3.1.2 REQUIRED )

if(iMSTK_MODULE_ENABLE_DevicesVRPN)
  # FTD2XX
  find_package( FTD2XX REQUIRED )
endif()

# g3log
find_package( g3log REQUIRED )

if (iMSTK_MODULE_ENABLE_DevicesVRPN)
  # LibNiFalcon
  find_package( LibNiFalcon REQUIRED)

  # Libusb
  find_package( Libusb REQUIRED)
endif()

# OpenVR
if (iMSTK_MODULE_ENABLE_RenderingVTK)
  find_package( OpenVR REQUIRED )
endif()

# SFML
# if (iMSTK_MODULE_ENABLE_SFML)
#   find_package( SFML REQUIRED )
# endif()

# TBB
find_package(TBB REQUIRED)

# VegaFEM
find_package( VegaFEM REQUIRED CONFIG )

# VRPN
if (iMSTK_MODULE_ENABLE_DevicesVRPN)
  find_package( VRPN REQUIRED )
endif()

# OpenHaptics
if(iMSTK_MODULE_ENABLE_DevicesOpenHaptics)
  find_package( OpenHapticsSDK REQUIRED )
endif()

# VTK
find_package(VTK REQUIRED CONFIG)

if (iMSTK_MODULE_ENABLE_RenderingVTK)
  find_package(vtkRenderingOpenVR QUIET)
  find_package(vtkRenderingExternal QUIET)
endif()

list(APPEND iMSTK_VTK_DEPENDENCIES
  ChartsCore
  CommonCore
  CommonDataModel
  FiltersCore
  FiltersSources
  ImagingCore
  ImagingGeneral
  ImagingStencil
  InteractionStyle)
if(iMSTK_MODULE_ENABLE_MeshIO)
  list(APPEND VTK_MODULE_SETTINGS
    IOGeometry
    IOImage
    IOImport
    IOLegacy
    IOPLY
    IOXML)
endif()
if (iMSTK_MODULE_ENABLE_RenderingVTK)
  list(APPEND iMSTK_VTK_DEPENDENCIES
    RenderingAnnotation
    RenderingContext2D
    RenderingContextOpenGL2
    RenderingCore
    RenderingExternal
    RenderingOpenGL2
    RenderingOpenVR
    RenderingVolume
    RenderingVolumeOpenGL2
    ViewsContext2D)
endif()
if (iMSTK_MODULE_ENABLE_Filtering)
  list(APPEND iMSTK_VTK_DEPENDENCIES
    ImagingMath
    FiltersExtraction
    FiltersGeneral
    FiltersModeling)
endif()
find_package(VTK COMPONENTS ${iMSTK_VTK_DEPENDENCIES} REQUIRED)
add_definitions(-D${iMSTK_MODULE_ENABLE_RenderingVTK})


#--------------------------------------------------------------------------
# External Utility Packages
#--------------------------------------------------------------------------
include(imstkAddExecutable)
# Uncrustify
find_program(Uncrustify_EXECUTABLE uncrustify)
# include(SetupUncrustifyConfig)
if(Uncrustify_EXECUTABLE)
  include(imstkAddUncrustifyCustomTarget)
else(Uncrustify_EXECUTABLE)
  message(WARNING "uncrustify not found! Cannot run code-style test.")
endif(Uncrustify_EXECUTABLE)


#-----------------------------------------------------------------------------
# Data
#----------------------------------------------------------------------
if (${PROJECT_NAME}_BUILD_TESTING OR ${PROJECT_NAME}_BUILD_EXAMPLES)
  add_definitions( -DiMSTK_DATA_ROOT=\"${CMAKE_INSTALL_PREFIX}/data/\")
  
  # Note if the target name or data name changes this has to be changed as well
  add_custom_target(CopyDataFiles ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_BINARY_DIR}/../External/iMSTKData/src/Data
    ${CMAKE_INSTALL_PREFIX}/data
    )
endif()


#-----------------------------------------------------------------------------
# Testing
#----------------------------------------------------------------------
if(${PROJECT_NAME}_BUILD_TESTING)
  find_package( GTest CONFIG REQUIRED )
endif()


#-----------------------------------------------------------------------------
# Code Coverage
#----------------------------------------------------------------------
if(${PROJECT_NAME}_COVERAGE)
  # Note HS Had to add ALL in setup_target_for_coverage to force the code coverage pass 
  set(CODE_COVERAGE_VERBOSE ON)
  include(CodeCoverage)
  append_coverage_compiler_flags()
  setup_target_for_coverage_lcov(
    NAME Coverage
    EXECUTABLE ctest -E "(CodeFormatTest|RenderTest.*)"
    EXCLUDE "${CMAKE_BINARY_DIR}/../install/include/*" "/usr/*" "*/Testing/*"
  )
endif()

# Folder name to put our headers/cmake config files under
set(${PROJECT_NAME}_INSTALL_FOLDER ${PROJECT_NAME}-${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR})


#--------------------------------------------------------------------------
# Add Source code subdirectories
#--------------------------------------------------------------------------
foreach(module ${iMSTK_MODULE_LIST})
  if (iMSTK_MODULE_ENABLE_${module})
    add_subdirectory(Source/${module})
  endif()
endforeach()
add_subdirectory(Source/Testing)

#--------------------------------------------------------------------------
# Add Examples subdirectories
#--------------------------------------------------------------------------
if(${PROJECT_NAME}_BUILD_EXAMPLES)
    add_subdirectory(Examples)
endif()

#--------------------------------------------------------------------------
# Testing
#--------------------------------------------------------------------------
if (${PROJECT_NAME}_BUILD_TESTING)
  add_subdirectory(Testing)
endif()

#--------------------------------------------------------------------------
# Add setup script for *nix systems
#--------------------------------------------------------------------------
if(NOT WIN32)
  # Create setup shell script to create an environment for running examples
  set(LIBRARY_PATH_VAR "LD_LIBRARY_PATH")
  if( APPLE )
    set(LIBRARY_PATH_VAR "DYLD_FALLBACK_LIBRARY_PATH")
  endif()
  configure_file(
    ${CMAKE_SOURCE_DIR}/CMake/setup_iMSTK.sh.in
    ${CMAKE_INSTALL_PREFIX}/setup_iMSTK.sh
    @ONLY)
endif()

#--------------------------------------------------------------------------
# Install Find Modules
#--------------------------------------------------------------------------
set(${PROJECT_NAME}_INSTALL_CONFIG_DIR "lib/cmake/${${PROJECT_NAME}_INSTALL_FOLDER}")
file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_CONFIG_DIR}/modules)
file(GLOB modules ${CMAKE_SOURCE_DIR}/CMake/Find*.cmake)
foreach(module ${modules})
  install(FILES ${module} DESTINATION "${${PROJECT_NAME}_INSTALL_CONFIG_DIR}/modules")
endforeach()
install(FILES ${CMAKE_SOURCE_DIR}/CMake/Utilities/imstkFind.cmake DESTINATION "${${PROJECT_NAME}_INSTALL_CONFIG_DIR}/modules")

#--------------------------------------------------------------------------
# Innerbuild dummy test
#--------------------------------------------------------------------------
add_test(
  NAME imstkDummyTest
  COMMAND ${CMAKE_COMMAND} -E echo "Success"
)

#--------------------------------------------------------------------------
# Export Targets
#--------------------------------------------------------------------------
# This variable controls the prefix used to generate the following files:
#  ${PROJECT_NAME}ConfigVersion.cmake
#  ${PROJECT_NAME}Config.cmake
#  ${PROJECT_NAME}Targets.cmake
# and it also used to initialize ${PROJECT_NAME}_INSTALL_CONFIG_DIR value.
set(export_config_name ${PROJECT_NAME})
string(TOLOWER "${PROJECT_NAME}" PROJECT_NAMESPACE)
set(PROJECT_NAMESPACE "${PROJECT_NAMESPACE}::")
#------------------------------------------------------------------------------
# Configure ${PROJECT_NAME}ConfigVersion.cmake common to build and install tree
include(CMakePackageConfigHelpers)
set(config_version_file ${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake)
write_basic_package_version_file(
  ${config_version_file}
  VERSION "${${PROJECT_NAME}_VERSION}"
  COMPATIBILITY AnyNewerVersion
  )
#------------------------------------------------------------------------------
# Export '${PROJECT_NAME}Targets.cmake' for a build tree
export(
  EXPORT ${PROJECT_NAME}Targets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake
  NAMESPACE ${PROJECT_NAMESPACE}
  )
# Configure '${PROJECT_NAME}Config.cmake' for a build tree
set(CONFIG_CODE "####### Expanded from \@CONFIG_CODE\@ #######\n")

set(CONFIG_CODE "${CONFIG_CODE}set(Assimp_ROOT_DIR \"${Assimp_ROOT_DIR}\")\n")
set(CONFIG_CODE "${CONFIG_CODE}set(Assimp_LIB_DIR \"${Assimp_LIB_DIR}\")\n")

set(CONFIG_CODE "${CONFIG_CODE}set(Eigen3_DIR \"${Eigen3_DIR}\")\n")

set(CONFIG_CODE "${CONFIG_CODE}set(g3log_ROOT_DIR \"${g3log_ROOT_DIR}\")\n")
set(CONFIG_CODE "${CONFIG_CODE}set(g3log_LIB_DIR \"${g3log_LIB_DIR}\")\n")

set(CONFIG_CODE "${CONFIG_CODE}set(GTest_DIR \"${GTest_DIR}\")\n")

if(iMSTK_MODULE_ENABLE_DevicesVRPN)
  set(CONFIG_CODE "${CONFIG_CODE}# Required by FindFTD2XX CMake module used in and provided by the iMSTK project\n")
  set(CONFIG_CODE "${CONFIG_CODE}set(FTD2XX_ROOT_DIR \"${FTD2XX_ROOT_DIR}\")\n")
  set(CONFIG_CODE "${CONFIG_CODE}set(FTD2XX_LIB_DIR \"${FTD2XX_LIB_DIR}\")\n")
  set(CONFIG_CODE "${CONFIG_CODE}# Required by Findlibftd2xx CMake module provided and used by the LibNiFalcon project\n")
  set(CONFIG_CODE "${CONFIG_CODE}set(LIBFTD2XX_INCLUDE_DIR \"${LIBFTD2XX_INCLUDE_DIR}\")\n")
  set(CONFIG_CODE "${CONFIG_CODE}set(LIBFTD2XX_LIBRARY \"${LIBFTD2XX_LIBRARY}\")\n")
endif()

set(CONFIG_CODE "${CONFIG_CODE}set(LibNiFalcon_ROOT_DIR \"${LibNiFalcon_ROOT_DIR}\")\n")
set(CONFIG_CODE "${CONFIG_CODE}set(LibNiFalcon_LIB_DIR \"${LibNiFalcon_LIB_DIR}\")\n")

if(iMSTK_MODULE_ENABLE_DevicesVRPN)
  set(CONFIG_CODE "${CONFIG_CODE}set(Libusb_ROOT_DIR \"${Libusb_ROOT_DIR}\")\n")
  set(CONFIG_CODE "${CONFIG_CODE}set(Libusb_LIB_DIR \"${Libusb_LIB_DIR}\")\n")
else()
  set(CONFIG_CODE "${CONFIG_CODE}set(Libusb_INCLUDE_DIR \"${Libusb_INCLUDE_DIR}\")\n")
  set(CONFIG_CODE "${CONFIG_CODE}set(Libusb_LIBRARY_libusb-1.0-RELEASE \"${Libusb_LIBRARY_libusb-1.0-RELEASE}\")\n")
  set(CONFIG_CODE "${CONFIG_CODE}set(Libusb_LIBRARY_libusb-1.0-DEBUG \"${Libusb_LIBRARY_libusb-1.0-DEBUG}\")\n")
endif()

if (iMSTK_MODULE_ENABLE_RenderingVTK)
  set(CONFIG_CODE "${CONFIG_CODE}# Add variables to CACHE to workaround issue related to use of\n")
  set(CONFIG_CODE "${CONFIG_CODE}# mark_as_advanced command in FindOpenVR and described in CMP0102\n")
  set(CONFIG_CODE "${CONFIG_CODE}# See https://cmake.org/cmake/help/latest/policy/CMP0102.html\n")
  set(CONFIG_CODE "${CONFIG_CODE}set(OpenVR_INCLUDE_DIR \"${OpenVR_INCLUDE_DIR}\" CACHE PATH \"\")\n")
  set(CONFIG_CODE "${CONFIG_CODE}set(OpenVR_LIBRARY \"${OpenVR_LIBRARY}\" CACHE FILEPATH \"\")\n")
  set(CONFIG_CODE "${CONFIG_CODE}set(OpenVR_ROOT_DIR \"${OpenVR_ROOT_DIR}\")\n")
  set(CONFIG_CODE "${CONFIG_CODE}set(OpenVR_LIB_DIR \"${OpenVR_LIB_DIR}\")\n")
endif()

set(CONFIG_CODE "${CONFIG_CODE}set(TBB_DIR \"${TBB_DIR}\")\n")

set(CONFIG_CODE "${CONFIG_CODE}set(VegaFEM_DIR \"${VegaFEM_DIR}\")\n")
set(CONFIG_CODE "${CONFIG_CODE}set(VTK_DIR \"${VTK_DIR}\")\n")
if (iMSTK_MODULE_ENABLE_RenderingVTK)
  set(CONFIG_CODE "${CONFIG_CODE}set(vtkRenderingOpenVR_DIR \"${vtkRenderingOpenVR_DIR}\")\n")
  set(CONFIG_CODE "${CONFIG_CODE}set(vtkRenderingExternal_DIR \"${vtkRenderingExternal_DIR}\")\n")
endif()

set(CONFIG_CODE "${CONFIG_CODE}##################################################")
set(build_config ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Config.cmake)
configure_package_config_file(
  CMake/${PROJECT_NAME}Config.cmake.in
  ${build_config}
  INSTALL_DESTINATION "${PROJECT_BINARY_DIR}"
  )

#------------------------------------------------------------------------------
# Export '${PROJECT_NAME}Targets.cmake' for an install tree
install(
  EXPORT ${PROJECT_NAME}Targets
  FILE ${PROJECT_NAME}Targets.cmake
  NAMESPACE ${PROJECT_NAMESPACE}
  DESTINATION ${${PROJECT_NAME}_INSTALL_CONFIG_DIR}
  )
# Configure '${PROJECT_NAME}Config.cmake' for an install tree
set(CONFIG_CODE "####### Expanded from \@CONFIG_CODE\@ #######\n")
# Note: iMSTK specific variables likes <proj>_ROOT_DIR and <proj>_LIB_DIR
# are not set because the location of the corresponding project is automatically
# inferred by Find<proj> CMake module provided by iSMTK.
set(CONFIG_CODE "${CONFIG_CODE}set(Eigen3_DIR \"${Eigen3_DIR}\")\n")
set(CONFIG_CODE "${CONFIG_CODE}set(GTest_DIR \"${GTest_DIR}\")\n")
set(CONFIG_CODE "${CONFIG_CODE}set(TBB_DIR \"${TBB_DIR}\")\n")
set(CONFIG_CODE "${CONFIG_CODE}set(VegaFEM_DIR \"${VegaFEM_DIR}\")\n")
set(CONFIG_CODE "${CONFIG_CODE}set(VTK_DIR \"${VTK_DIR}\")\n")
set(CONFIG_CODE "${CONFIG_CODE}# Add variables to CACHE to workaround issue related to use of\n")
set(CONFIG_CODE "${CONFIG_CODE}# mark_as_advanced command in FindOpenVR and described in CMP0102\n")
set(CONFIG_CODE "${CONFIG_CODE}# See https://cmake.org/cmake/help/latest/policy/CMP0102.html\n")
set(CONFIG_CODE "${CONFIG_CODE}set(OpenVR_INCLUDE_DIR \"${OpenVR_INCLUDE_DIR}\" CACHE PATH \"\")\n")
set(CONFIG_CODE "${CONFIG_CODE}set(OpenVR_LIBRARY \"${OpenVR_LIBRARY}\" CACHE FILEPATH \"\")\n")
set(CONFIG_CODE "${CONFIG_CODE}##################################################")

set(install_config ${PROJECT_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}Config.cmake)
configure_package_config_file(
  CMake/${PROJECT_NAME}Config.cmake.in 
  ${install_config}
  INSTALL_DESTINATION ${${PROJECT_NAME}_INSTALL_CONFIG_DIR}
  )

# Install config files
install(
  FILES ${config_version_file} ${install_config}
  DESTINATION "${${PROJECT_NAME}_INSTALL_CONFIG_DIR}"
  )
