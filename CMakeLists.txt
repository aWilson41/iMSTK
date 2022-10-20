cmake_minimum_required(VERSION 3.9)

option(iMSTK_BUILD_FOR_ANDROID "Build with Android NDK toolchain" OFF)

#-----------------------------------------------------------------------------
# Android configuration
#-----------------------------------------------------------------------------
# Android config must set the toolchain file before project is defined
if (iMSTK_BUILD_FOR_ANDROID)
  # Android options
  set (_ANDROID_NDK_DEFAULT "/opt/android-ndk")
  if (DEFINED ENV{ANDROID_NDK})
    set (_ANDROID_NDK_DEFAULT "$ENV{ANDROID_NDK}")
  endif()
  set(ANDROID_NDK ${_ANDROID_NDK_DEFAULT} CACHE PATH
    "Set to the absolute path of the Android NDK root directory.\
     A \$\{ANDROID_NDK\}/sources/android directory must exist."
  )
  # Verify ANDROID_NDK path
  if (NOT EXISTS "${ANDROID_NDK}/sources/android")
    message(FATAL_ERROR "Please set a valid ANDROID_NDK path")
  endif()
  set(ANDROID_NATIVE_API_LEVEL "27" CACHE STRING "Android Native API Level")
  set(ANDROID_ARCH_ABI "arm64-v8a" CACHE STRING "Target Android architecture/abi")

  # Setup a toolchain file
  set(_ANDROID_DIR "imstk-android-${ANDROID_NATIVE_API_LEVEL}-${ANDROID_ARCH_ABI}")
  set(_ANDROID_TOOLCHAIN ${CMAKE_CURRENT_BINARY_DIR}/${_ANDROID_DIR}-toolchain.cmake)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMake/iMSTKAndroid-toolchain.cmake.in ${_ANDROID_TOOLCHAIN} @ONLY)
  set(CMAKE_TOOLCHAIN_FILE ${CMAKE_CURRENT_BINARY_DIR}/${_ANDROID_DIR}-toolchain.cmake CACHE STRING "Override option" FORCE)
endif()

project(iMSTK VERSION 6.0.0)

include(CMakeDependentOption)

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

# Unity Switch
if (WIN32 OR ANDROID)
  option(${PROJECT_NAME}_BUILD_FOR_UNITY "Build iMSTK libraries for the unity asset" OFF)
  mark_as_superbuild(${PROJECT_NAME}_BUILD_FOR_UNITY:BOOL)
endif()

option(${PROJECT_NAME}_BUILD_TESTING "Build iMSTK tests" ON)
mark_as_superbuild(${PROJECT_NAME}_BUILD_TESTING:BOOL)

# Couldn't do earlier since CMAKE_MODULE_PATH hadn't given solve dependencies include
mark_as_superbuild(${PROJECT_NAME}_BUILD_FOR_ANDROID:BOOl)
mark_as_superbuild(ANDROID_NDK:STRING)
mark_as_superbuild(ANDROID_NATIVE_API_LEVEL:STRING)
mark_as_superbuild(ANDROID_ARCH_ABI:STRING)

# SWIG 
cmake_dependent_option(${PROJECT_NAME}_WRAP_CSHARP "Build iMSTK-C# wrapper code and lib" OFF "NOT ${PROJECT_NAME}_BUILD_FOR_UNITY" ON)
mark_as_superbuild(${PROJECT_NAME}_WRAP_CSHARP:BOOL)

# SWIG will generate managed array using pinning if this option is ON. Otherwise using P/Invoke default array marshalling.
option(${PROJECT_NAME}_SWIG_PINNED_ARRAY "Managed C# arrays using pinning" ON)
mark_as_advanced(${PROJECT_NAME}_SWIG_PINNED_ARRAY)
mark_as_superbuild(${PROJECT_NAME}_SWIG_PINNED_ARRAY:BOOL)

if (${PROJECT_NAME}_WRAP_CSHARP)
  find_package(SWIG REQUIRED)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)

  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
endif()

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

# If on, everything is built as shared libs that can be
# If off, everything is built as static libs
# Force off if build for android
cmake_dependent_option(BUILD_SHARED_LIBS "Build as many libraries as possible as shared or static" ON "NOT ${PROJECT_NAME}_BUILD_FOR_ANDROID" OFF)
mark_as_superbuild(BUILD_SHARED_LIBS:BOOL)

# General Options
option(${PROJECT_NAME}_BUILD_EXAMPLES "Build iMSTK examples" ON)
mark_as_superbuild(${PROJECT_NAME}_BUILD_EXAMPLES:BOOL)

# Build Benchmark
option(${PROJECT_NAME}_BUILD_BENCHMARK "Build iMSTK benchmarks" OFF)
mark_as_superbuild(${PROJECT_NAME}_BUILD_BENCHMARK:BOOL)

# Default on, force off if not win32 or on unity
cmake_dependent_option(${PROJECT_NAME}_USE_SFML "If on, imstk is built with SFML audio support" ON "WIN32; NOT ${PROJECT_NAME}_BUILD_FOR_UNITY" OFF)
mark_as_superbuild(${PROJECT_NAME}_USE_SFML:BOOL)

# Force off if build for android
cmake_dependent_option(${PROJECT_NAME}_USE_VegaFEM "Build with VegaFEM models" ON "NOT ${PROJECT_NAME}_BUILD_FOR_ANDROID" OFF)
mark_as_superbuild(${PROJECT_NAME}_USE_VegaFEM:BOOL)

# Default off, if not building VegaFEM or not unix force off
cmake_dependent_option(${PROJECT_NAME}_USE_MODEL_REDUCTION "Build with model reduction, requires a VegaFEM built with Intel MKL and arpack" OFF "${PROJECT_NAME}_USE_VegaFEM;UNIX" OFF)
mark_as_superbuild(${PROJECT_NAME}_USE_MODEL_REDUCTION:BOOL)

if (UNIX)
  option(${PROJECT_NAME}_COLOR_OUTPUT "For coloring output in gcc" OFF)
  mark_as_superbuild(${PROJECT_NAME}_COLOR_OUTPUT:BOOL)

  option(${PROJECT_NAME}_COVERAGE "Build for calculating coverage" OFF)
  mark_as_superbuild(${PROJECT_NAME}_COVERAGE:BOOL)

  if (${PROJECT_NAME}_COVERAGE AND NOT ${CMAKE_BUILD_TYPE} MATCHES [Dd][Ee][Bb][Uu][Gg])
    message(WARNING "Coverage requires a debug build forcing debug")
    set(CMAKE_BUILD_TYPE "DEBUG")
  endif()  
endif()

# Default on, force off if build for android
cmake_dependent_option(${PROJECT_NAME}_USE_TBB "Build with TBB parallelization" ON "NOT ${PROJECT_NAME}_BUILD_FOR_ANDROID" OFF)
mark_as_superbuild(${PROJECT_NAME}_USE_TBB:BOOL)

option(${PROJECT_NAME}_USE_VRPN "Build with VRPN support" OFF)
mark_as_superbuild(${PROJECT_NAME}_USE_VRPN:BOOL)

# Default off, force off if build for android
cmake_dependent_option(${PROJECT_NAME}_USE_HAPLY "Build with Haply support" OFF "NOT ${PROJECT_NAME}_BUILD_FOR_ANDROID" OFF)
if (${PROJECT_NAME}_USE_HAPLY)
  message(WARNING "iMSTK_USE_HAPLY set to true, this code is based of an early version of the Haply API")
endif()

# Default on, if building unity force off
cmake_dependent_option(${PROJECT_NAME}_USE_RENDERING_VTK "Build with/out VTK rendering support" ON "NOT ${PROJECT_NAME}_BUILD_FOR_UNITY" OFF)
mark_as_superbuild(${PROJECT_NAME}_USE_RENDERING_VTK:BOOL)

# Default off, if not rendering or unix force off
cmake_dependent_option(${PROJECT_NAME}_USE_VTK_OSMESA "Build with VTK OSMesa rendering" OFF "${PROJECT_NAME}_USE_RENDERING_VTK;UNIX" OFF)
mark_as_superbuild(${PROJECT_NAME}_USE_VTK_OSMESA:BOOL)

# Default off, if not win32 or building for android force off
cmake_dependent_option(${PROJECT_NAME}_USE_OpenHaptics "Use OpenHaptic Support" OFF "WIN32;NOT ${PROJECT_NAME}_BUILD_FOR_ANDROID" OFF)
mark_as_superbuild(${PROJECT_NAME}_USE_OpenHaptics:BOOL)

# If at least one haptic API is provided
if (${PROJECT_NAME}_USE_HAPLY OR ${PROJECT_NAME}_USE_OpenHaptics OR ${PROJECT_NAME}_USE_VRPN)
  set(${PROJECT_NAME}_USE_HAPTICS ON)
endif()
mark_as_superbuild(${PROJECT_NAME}_USE_HAPTICS:BOOL)

# Default on, if building for unity or android force off
cmake_dependent_option(${PROJECT_NAME}_USE_G3LOG "Use g3 async logger. If off, builtin sync logger is used" ON "NOT ${PROJECT_NAME}_BUILD_FOR_UNITY;NOT ${PROJECT_NAME}_BUILD_FOR_ANDROID" OFF)
mark_as_superbuild(${PROJECT_NAME}_USE_G3LOG:BOOL)

# Default on, if not rendering or not testing force off
cmake_dependent_option(${PROJECT_NAME}_BUILD_VISUAL_TESTING "Build iMSTK visual tests" ON "${PROJECT_NAME}_USE_RENDERING_VTK;${PROJECT_NAME}_BUILD_TESTING" OFF)
mark_as_superbuild(${PROJECT_NAME}_BUILD_VISUAL_TESTING:BOOL)

#-----------------------------------------------------------------------------
# CTest/Dashboards
#-----------------------------------------------------------------------------
set(BUILD_TESTING OFF)
if (${PROJECT_NAME}_BUILD_TESTING)
  set(BUILD_TESTING ON)

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

  # Default on, force off if build with android
  cmake_dependent_option(${PROJECT_NAME}_USE_Uncrustify "Use Uncrustify as a code style beautifier." ON "NOT ${PROJECT_NAME}_BUILD_FOR_ANDROID" OFF)
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
  if (${PROJECT_NAME}_USE_G3LOG)
    imstk_define_dependency(g3log)
  endif()
  if (${PROJECT_NAME}_USE_RENDERING_VTK)
    imstk_define_dependency(OpenVR)
  endif()
  if (${PROJECT_NAME}_USE_TBB)
    imstk_define_dependency(TBB)
  endif()
  if (${PROJECT_NAME}_USE_VegaFEM)
    imstk_define_dependency(VegaFEM)
  endif()
  imstk_define_dependency(VTK)
  if(${PROJECT_NAME}_USE_OpenHaptics)
    imstk_define_dependency(OpenHaptics)
  endif()
  if (${PROJECT_NAME}_USE_VRPN)
    if(WIN32)
      imstk_define_dependency(Libusb) #for VRPN
    endif()
    imstk_define_dependency(VRPN)
  endif()
  if (${PROJECT_NAME}_USE_HAPLY)
    imstk_define_dependency(Haply)
  endif()

  if(${PROJECT_NAME}_USE_SFML)
    imstk_define_dependency(SFML)
  endif()

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

  if (${PROJECT_NAME}_BUILD_TESTING OR 
      ${PROJECT_NAME}_BUILD_EXAMPLES)
    imstk_define_dependency(iMSTKData)
  endif()


  # Build with Google Benchmark for benchmarking
  if (${PROJECT_NAME}_BUILD_BENCHMARK)
    imstk_define_dependency(benchmark)
  endif()

  #-----------------------------------------------------------------------------
  # Solve project dependencies
  #-----------------------------------------------------------------------------
  # Call CMakeLists.txt in CMake/External which will solve the dependencies
  # and add the External projects, including this one: this top-level
  # CMakeLists.txt will be called back with SUPERBUILD=OFF, to execute
  # the rest of the code below (INNERBUILD), which explains the `return`
  add_subdirectory(CMake/External)

  # Wrappers are built after everything
  if (${PROJECT_NAME}_WRAP_CSHARP)
    include(ExternalProject)
    ExternalProject_Add(CWrapper
      DEPENDS
        iMSTK
      DOWNLOAD_COMMAND ""
      PREFIX "${CMAKE_CURRENT_BINARY_DIR}/CWrapper"
      SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/Wrappers"
      BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/CWrapper/build"
      STAMP_DIR "${CMAKE_CURRENT_BINARY_DIR}/CWrapper/stamp"
      TMP_DIR "${CMAKE_CURRENT_BINARY_DIR}/CWrapper/tmp"
      INSTALL_DIR "install"
      CMAKE_CACHE_ARGS
        -DBUILD_SHARED_LIBS:BOOL=True # You can only dllimport shared libraries
        -DiMSTK_DIR:PATH=${CMAKE_BINARY_DIR}/Innerbuild
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_PREFIX_PATH}
        -DiMSTK_SWIG_PINNED_ARRAY:BOOL=${iMSTK_SWIG_PINNED_ARRAY}
        -DiMSTK_USE_RENDERING_VTK:BOOL=${iMSTK_USE_RENDERING_VTK}
        -DiMSTK_USE_HAPTICS:BOOL=${iMSTK_USE_HAPTICS}
        -DiMSTK_USE_OpenHaptics:BOOL=${iMSTK_USE_OpenHaptics}
        -DiMSTK_USE_VRPN:BOOL=${iMSTK_USE_VRPN}
        -DSWIG_DIR:PATH=${SWIG_DIR}
        -DSWIG_EXECUTABLE:PATH=${SWIG_EXECUTABLE}
        -DSWIG_CSHARP_OUTPUT_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/CSharpWrapper/src
        -DCMAKE_TOOLCHAIN_FILE:PATH=${CMAKE_TOOLCHAIN_FILE}
        -DANDROID_NDK:PATH=${ANDROID_NDK}
        -DANDROID_NATIVE_API_LEVEL:STRING=${ANDROID_NATIVE_API_LEVEL}
        -DANDROID_ARCH_ABI:STRING=${ANDROID_ARCH_ABI}
      )
    set_target_properties(CWrapper PROPERTIES FOLDER iMSTK)
    # CMake only supports .NET csharp, not mono which could be used on linux
    if (MSVC OR ${PROJECT_NAME}_BUILD_FOR_ANDROID)
      # If building for android we have an android toolchain but C# .NET needs visual studio
      # So we need to switch toolchains here
      #if (${PROJECT_NAME}_BUILD_FOR_ANDROID)
      #  message(WARNING "Cannot build iMSTKSharp for android yet")
      #else()
        ExternalProject_Add(CSharpWrapper
          DEPENDS
            CWrapper
          DOWNLOAD_COMMAND ""
          PREFIX "${CMAKE_CURRENT_BINARY_DIR}/CSharpWrapper"
          SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Wrappers/iMSTKCSharpWrapper" # Should prob install the CMake file into src along with ".cs" files
          BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/CSharpWrapper/build"
          STAMP_DIR "${CMAKE_CURRENT_BINARY_DIR}/CSharpWrapper/stamp"
          TMP_DIR "${CMAKE_CURRENT_BINARY_DIR}/CSharpWrapper/tmp"
          INSTALL_DIR "install"
          CMAKE_GENERATOR "Visual Studio 17 2022" # How can I ask for the latest version of visual studio available?
          CMAKE_CACHE_ARGS
            -DBUILD_SHARED_LIBS:BOOL=TRUE  # C# cannot dll import static libs, so always build shared
            -DiMSTKCWrapper_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/CWrapper/build
            -DSWIG_OUTPUT_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/CSharpWrapper/src
            -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
            -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_PREFIX_PATH}
            -DiMSTK_SWIG_PINNED_ARRAY:BOOL=${iMSTK_SWIG_PINNED_ARRAY}
            -DiMSTK_USE_RENDERING_VTK:BOOL=${iMSTK_USE_RENDERING_VTK}
            -DiMSTK_USE_HAPTICS:BOOL=${iMSTK_USE_HAPTICS}
            -DiMSTK_USE_OpenHaptics:BOOL=${iMSTK_USE_OpenHaptics}
            -DiMSTK_USE_VRPN:BOOL=${iMSTK_USE_VRPN}
          )
        set_target_properties(CSharpWrapper PROPERTIES FOLDER iMSTK)
      #endif()
    elseif(UNIX)
      message(WARNING "C# library/.so cannot be built on unix systems through cmake, please build manually using generated \".cs\" files")
    #find_package(Mono) # Provided with cmake to find mono bits
    endif()
  endif()

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

if (${PROJECT_NAME}_USE_G3LOG)
  # g3log
  find_package( g3log REQUIRED )
  remove_definitions(-DiMSTK_SYNCHRONOUS_LOGGING)
else()
  add_definitions(-DiMSTK_SYNCHRONOUS_LOGGING)
endif()

if (${PROJECT_NAME}_USE_VRPN)
  # Libusb
  find_package( Libusb REQUIRED)
endif()

if (${PROJECT_NAME}_USE_RENDERING_VTK)
  # OpenVR
  find_package( OpenVR REQUIRED )
endif()

# SFML
if (${PROJECT_NAME}_USE_SFML)
  find_package(SFML REQUIRED)
  add_definitions(-DiMSTK_USE_SFML)
else()
  remove_definitions(-DiMSTK_USE_SFML)
endif()

# TBB
if (${PROJECT_NAME}_USE_TBB)
  find_package(TBB REQUIRED)
  add_definitions(-DiMSTK_USE_TBB)
else()
  remove_definitions(-DiMSTK_USE_TBB)
endif()

# VegaFEM
if (${PROJECT_NAME}_USE_VegaFEM)
  find_package(VegaFEM REQUIRED CONFIG)
  add_definitions(-DiMSTK_USE_VegaFEM)
else()
  remove_definitions(-DiMSTK_USE_VegaFEM)
endif()

# VRPN
if (${PROJECT_NAME}_USE_VRPN)
  find_package( VRPN REQUIRED )
  add_definitions( -DiMSTK_USE_VRPN )
else()
  remove_definitions( -DiMSTK_USE_VRPN )
endif()

if (${PROJECT_NAME}_USE_HAPLY)
  find_package( Haply REQUIRED )
  add_definitions( -DiMSTK_USE_HAPLY )
else()
  remove_definitions( -DiMSTK_USE_HAPLY )
endif()

# OpenHaptics
if(${PROJECT_NAME}_USE_OpenHaptics)
  find_package( OpenHapticsSDK REQUIRED )
  add_definitions( -DiMSTK_USE_OPENHAPTICS )
else()
  remove_definitions( -DiMSTK_USE_OPENHAPTICS )
endif()

# VTK
find_package(VTK REQUIRED CONFIG)
# modules are linked via `VTK::CommonCore`
# vtk_module_autoinit is needed

if (${PROJECT_NAME}_USE_RENDERING_VTK)
  find_package(vtkRenderingOpenVR QUIET)
  find_package(vtkRenderingExternal QUIET)
endif()

if (${PROJECT_NAME}_USE_HAPTICS)
  add_definitions(-DiMSTK_USE_HAPTICS)
else()
  remove_definitions(-DiMSTK_USE_HAPTICS)
endif()

list(APPEND iMSTK_VTK_DEPENDENCIES
  ChartsCore
  CommonCore
  CommonDataModel
  FiltersCore
  FiltersExtraction
  FiltersGeneral
  FiltersModeling
  FiltersSources
  IOExport
  IOImport
  IOPLY
  IOParallel
  IOParallelXML
  ImagingCore
  ImagingGeneral
  ImagingStencil
  InteractionStyle)
if (${PROJECT_NAME}_USE_RENDERING_VTK)
  list(APPEND iMSTK_VTK_DEPENDENCIES
    RenderingAnnotation
    RenderingContext2D
    RenderingContextOpenGL2
    RenderingCore
    RenderingExternal
    RenderingLabel
    RenderingOpenGL2
    RenderingOpenVR
    RenderingVolume
    RenderingVolumeOpenGL2
    ViewsContext2D)
endif()
find_package(VTK COMPONENTS ${iMSTK_VTK_DEPENDENCIES} REQUIRED)
if(${PROJECT_NAME}_USE_VTK_OSMESA)
  add_definitions(-DiMSTK_USE_VTK_OSMESA)
else()
  remove_definitions(-DiMSTK_USE_VTK_OSMESA)
endif()

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
  if (UNIX)
    # Don't use relative directory on Linux
    add_definitions( -DiMSTK_DATA_ROOT=\"${CMAKE_INSTALL_PREFIX}/data/\")
  else()
    # Use relative directory to bin for movement of the install directory
    add_definitions( -DiMSTK_DATA_ROOT=\"../data/\")
  endif()
  
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

#------------------------------------------------------------------
# Benchmarking
#------------------------------------------------------------------
if (${PROJECT_NAME}_BUILD_BENCHMARK)
  find_package( benchmark REQUIRED )
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
    EXECUTABLE ctest -E "(CodeFormatTest)"
    EXCLUDE "${CMAKE_BINARY_DIR}/../install/include/*" "/usr/*" "*/Testing/*" "*/Examples/*" "*/Wrappers/*"
  )
endif()

# Folder name to put our headers/cmake config files under
set(${PROJECT_NAME}_INSTALL_FOLDER ${PROJECT_NAME}-${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR})

#--------------------------------------------------------------------------
# Add Source code subdirectories
#--------------------------------------------------------------------------
add_subdirectory(Source/Common)
add_subdirectory(Source/ComponentModel)
add_subdirectory(Source/Materials)
add_subdirectory(Source/Geometry)
add_subdirectory(Source/MeshIO)
add_subdirectory(Source/GeometryMappers)
add_subdirectory(Source/DataStructures)
add_subdirectory(Source/Constraint)
add_subdirectory(Source/Devices)
add_subdirectory(Source/RenderingCore)
if (${PROJECT_NAME}_USE_RENDERING_VTK)
  add_subdirectory(Source/RenderingVTK)
endif()
add_subdirectory(Source/Solvers)
add_subdirectory(Source/DynamicalModels)
if (${PROJECT_NAME}_USE_VegaFEM)
  add_subdirectory(Source/DynamicalModelsVegaFEM)
endif()
add_subdirectory(Source/SceneEntities)
add_subdirectory(Source/Animation)
add_subdirectory(Source/Controllers)
add_subdirectory(Source/CollisionDetection)
add_subdirectory(Source/CollisionHandling)
add_subdirectory(Source/Needle)
add_subdirectory(Source/Scene)
add_subdirectory(Source/SimulationManager)
add_subdirectory(Source/ViewerCore)
if (${PROJECT_NAME}_USE_RENDERING_VTK)
  add_subdirectory(Source/ViewerVTK)
endif()
if( ${PROJECT_NAME}_BUILD_TESTING )
  add_subdirectory(Source/Testing)
  if (${PROJECT_NAME}_BUILD_VISUAL_TESTING)
    add_subdirectory(Source/VisualTesting)
  endif()
endif()
add_subdirectory(Source/Filtering)
add_subdirectory(Source/FilteringCore)

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

if (${PROJECT_NAME}_USE_G3LOG)
  set(CONFIG_CODE "${CONFIG_CODE}set(g3log_ROOT_DIR \"${g3log_ROOT_DIR}\")\n")
  set(CONFIG_CODE "${CONFIG_CODE}set(g3log_LIB_DIR \"${g3log_LIB_DIR}\")\n")
endif()

set(CONFIG_CODE "${CONFIG_CODE}set(GTest_DIR \"${GTest_DIR}\")\n")

set(CONFIG_CODE "${CONFIG_CODE}set(benchmark_DIR \"${benchmark_DIR}\")\n")

if(WIN32 AND ${PROJECT_NAME}_USE_VRPN)
  set(CONFIG_CODE "${CONFIG_CODE}set(Libusb_ROOT_DIR \"${Libusb_ROOT_DIR}\")\n")
  set(CONFIG_CODE "${CONFIG_CODE}set(Libusb_LIB_DIR \"${Libusb_LIB_DIR}\")\n")
endif()

if (${PROJECT_NAME}_USE_RENDERING_VTK)
  set(CONFIG_CODE "${CONFIG_CODE}# Add variables to CACHE to workaround issue related to use of\n")
  set(CONFIG_CODE "${CONFIG_CODE}# mark_as_advanced command in FindOpenVR and described in CMP0102\n")
  set(CONFIG_CODE "${CONFIG_CODE}# See https://cmake.org/cmake/help/latest/policy/CMP0102.html\n")
  set(CONFIG_CODE "${CONFIG_CODE}set(OpenVR_INCLUDE_DIR \"${OpenVR_INCLUDE_DIR}\" CACHE PATH \"\")\n")
  set(CONFIG_CODE "${CONFIG_CODE}set(OpenVR_LIBRARY \"${OpenVR_LIBRARY}\" CACHE FILEPATH \"\")\n")
  set(CONFIG_CODE "${CONFIG_CODE}set(OpenVR_ROOT_DIR \"${OpenVR_ROOT_DIR}\")\n")
  set(CONFIG_CODE "${CONFIG_CODE}set(OpenVR_LIB_DIR \"${OpenVR_LIB_DIR}\")\n")
endif()

if (${PROJECT_NAME}_USE_TBB)
  set(CONFIG_CODE "${CONFIG_CODE}set(TBB_DIR \"${TBB_DIR}\")\n")
endif()

if (${PROJECT_NAME}_USE_VegaFEM)
  set(CONFIG_CODE "${CONFIG_CODE}set(VegaFEM_DIR \"${VegaFEM_DIR}\")\n")
endif()
set(CONFIG_CODE "${CONFIG_CODE}set(VTK_DIR \"${VTK_DIR}\")\n")
if (${PROJECT_NAME}_USE_RENDERING_VTK)
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
if (${PROJECT_NAME}_USE_TBB)
  set(CONFIG_CODE "${CONFIG_CODE}set(TBB_DIR \"${TBB_DIR}\")\n")
endif()
if (${PROJECT_NAME}_USE_VegaFEM)
  set(CONFIG_CODE "${CONFIG_CODE}set(VegaFEM_DIR \"${VegaFEM_DIR}\")\n")
endif()
set(CONFIG_CODE "${CONFIG_CODE}set(VTK_DIR \"${VTK_DIR}\")\n")
set(CONFIG_CODE "${CONFIG_CODE}set(benchmark_DIR \"${benchmark_DIR}\")\n")
if (${PROJECT_NAME}_USE_RENDERING_VTK)
  set(CONFIG_CODE "${CONFIG_CODE}# Add variables to CACHE to workaround issue related to use of\n")
  set(CONFIG_CODE "${CONFIG_CODE}# mark_as_advanced command in FindOpenVR and described in CMP0102\n")
  set(CONFIG_CODE "${CONFIG_CODE}# See https://cmake.org/cmake/help/latest/policy/CMP0102.html\n")
  set(CONFIG_CODE "${CONFIG_CODE}set(OpenVR_INCLUDE_DIR \"${OpenVR_INCLUDE_DIR}\" CACHE PATH \"\")\n")
  set(CONFIG_CODE "${CONFIG_CODE}set(OpenVR_LIBRARY \"${OpenVR_LIBRARY}\" CACHE FILEPATH \"\")\n")
  set(CONFIG_CODE "${CONFIG_CODE}##################################################")
endif()

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
