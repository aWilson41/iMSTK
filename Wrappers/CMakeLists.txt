cmake_minimum_required(VERSION 3.9)
project(iMSTKWrappers)

message("size: ${CMAKE_SIZEOF_VOID_P}")
find_package(iMSTK REQUIRED)

# \todo: A different c wrapper is required to be generated per language
set(SWIG_INTERFACE_DIR ${CMAKE_CURRENT_LIST_DIR}/SwigInterface)

# include_directories(${CMAKE_CURRENT_SOURCE_DIR})
set_source_files_properties(${SWIG_INTERFACE_DIR}/imstkCWrapper.i PROPERTIES CPLUSPLUS ON)

set(LIB_TYPE SHARED)
if (NOT ${BUILD_SHARED_LIBS})
  set(LIB_TYPE STATIC)
endif()

include(UseSWIG)
swig_add_library(iMSTKCWrapper
    TYPE ${LIB_TYPE}
    OUTPUT_DIR ${SWIG_CSHARP_OUTPUT_DIR}
    LANGUAGE csharp
    SOURCES
    ${SWIG_INTERFACE_DIR}/imstkCWrapper.i)

# if (NOT iMSTK_BUILD_FOR_ANDROID)
#     # Always clear the ".cs" files from the build directory before generating
#     add_custom_command(
#         TARGET iMSTKCWrapper
#         PRE_BUILD
#         COMMAND ${CMAKE_COMMAND} -Ddir:PATH=${SWIG_CSHARP_OUTPUT_DIR}/ -P ${CMAKE_CURRENT_LIST_DIR}/ClearDirectory.cmake
#         DEPENDS ${PROJECT_NAME})
# endif()

# These files should not be provided via swig_add_library SOURCES
# But still useful to display in visual studios
if (MSVC) 
    set(OTHER_FILES
		${SWIG_INTERFACE_DIR}/callback.i
		${SWIG_INTERFACE_DIR}/common.i
		${SWIG_INTERFACE_DIR}/debug.i
		${SWIG_INTERFACE_DIR}/except.i
		${SWIG_INTERFACE_DIR}/ignored.i
		${SWIG_INTERFACE_DIR}/modifiers.i
		${SWIG_INTERFACE_DIR}/shared_ptr_instantiation.i
		${SWIG_INTERFACE_DIR}/std_function.i
		${SWIG_INTERFACE_DIR}/type_cast.i
		${SWIG_INTERFACE_DIR}/weak_ptr.i)

    # \todo: I think there's some property we can use to indicate to visual studios
    # that it should recompile if one of them changes, this must internally
    # happen in swig_add_library somewhere
    target_sources(iMSTKCWrapper PRIVATE ${OTHER_FILES})
    set_source_files_properties(${OTHER_FILES} PROPERTIES HEADER_FILE_ONLY TRUE)
endif()

set_target_properties(iMSTKCWrapper PROPERTIES
    SWIG_COMPILE_OPTIONS "-namespace;Imstk;-Wall"
    FOLDER Wrappers)

# Forward iMSTK variables to SWIG
set(imstk_SWIG_DEFINITIONS "")
if(iMSTK_SWIG_PINNED_ARRAY)
    list(APPEND imstk_SWIG_DEFINITIONS "SWIG_PINNED_ARRAY")
endif()
if (iMSTK_USE_RENDERING_VTK)
    list(APPEND imstk_SWIG_DEFINITIONS "iMSTK_USE_RENDERING_VTK")
endif()
if(iMSTK_USE_OpenHaptics)
    list(APPEND imstk_SWIG_DEFINITIONS "iMSTK_USE_OpenHaptics")
endif()
if(iMSTK_USE_HAPLY)
    list(APPEND imstk_SWIG_DEFINITIONS "iMSTK_USE_HAPLY")
endif()
if(iMSTK_USE_HAPTICS)
    list(APPEND imstk_SWIG_DEFINITIONS "iMSTK_USE_HAPTICS")
endif()
if(iMSTK_USE_VRPN)
    list(APPEND imstk_SWIG_DEFINITIONS "iMSTK_USE_VRPN")
endif()
if(iMSTK_USE_VegaFEM)
    list(APPEND imstk_SWIG_DEFINITIONS "iMSTK_USE_VegaFEM")
endif()
if(NOT iMSTK_USE_G3LOG)
    list(APPEND imstk_SWIG_DEFINITIONS "iMSTK_SYNCHRONOUS_LOGGING")
endif()

set_property(TARGET iMSTKCWrapper PROPERTY SWIG_COMPILE_DEFINITIONS ${imstk_SWIG_DEFINITIONS})

target_compile_definitions(iMSTKCWrapper PRIVATE  ${imstk_SWIG_DEFINITIONS})

# Link C Wrapper to the imstk libraries
set(CWRAPPER_DEPENDENCIES
    imstk::Common
    imstk::Materials
	imstk::Geometry
	imstk::DataStructures
	imstk::Devices
	imstk::Filtering
	imstk::FilteringCore
	imstk::Materials
	imstk::RenderingCore
	imstk::Solvers
	imstk::DynamicalModels
	imstk::CollisionDetection
	imstk::CollisionHandling
	imstk::SceneEntities
    imstk::Needle
	imstk::Scene
	imstk::SimulationManager
    imstk::ViewerCore
	imstk::Constraints
	imstk::Animation
    )
if (iMSTK_USE_VegaFEM)
    list(APPEND CWRAPPER_DEPENDENCIES DynamicalModelsVegaFEM)
endif()
if (iMSTK_USE_RENDERING_VTK)
    list(APPEND CWRAPPER_DEPENDENCIES
        RenderingVTK
        ViewerVTK
        )
endif()
target_link_libraries(iMSTKCWrapper PRIVATE ${CWRAPPER_DEPENDENCIES})
target_compile_options(iMSTKCWrapper PRIVATE
    $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
        -Wall -Wno-unused-function>
    $<$<CXX_COMPILER_ID:MSVC>:
        -W4 -MP -wd4505 -wd4127 /bigobj>)

# Ignoring warnings 
# MSVC C4505 Unreferenced function has been removed
# MSVC C4127 Caused by eigen "conditional expression is constant"


# Install the c wrapper
install( TARGETS iMSTKCWrapper EXPORT iMSTKCWrapperTargets
    RUNTIME DESTINATION bin COMPONENT RuntimeLibraries
    LIBRARY DESTINATION lib COMPONENT RuntimeLibraries
    ARCHIVE DESTINATION lib COMPONENT Development)

# if (UNIX AND NOT APPLE AND NOT iMSTK_BUILD_ANDROID)
# 	# Why use CODE and then execute COMMAND?
#     install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink iMSTKCWrapper.so ${CMAKE_INSTALL_PREFIX}/lib/libiMSTKCWrapper.so)")
# endif()

# Export the target file for the config file to use later			
export(EXPORT iMSTKCWrapperTargets
	FILE ${CMAKE_CURRENT_BINARY_DIR}/iMSTKCWrapperTargets.cmake)

configure_file(iMSTKCWrapperConfig.cmake.in
	"${CMAKE_CURRENT_BINARY_DIR}/iMSTKCWrapperConfig.cmake"
	@ONLY)

# Clear the ".cs" wrapper code in the install directory and install new ones
install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -Dsource_dir:PATH=${SWIG_CSHARP_OUTPUT_DIR} -Ddest_dir:PATH=${CMAKE_INSTALL_PREFIX}/include/iMSTKSharp/ -P ${CMAKE_CURRENT_LIST_DIR}/ClearInstall.cmake)")