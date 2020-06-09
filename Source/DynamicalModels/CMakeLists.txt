#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------
include(imstkAddLibrary)

imstk_add_library( imstkDynamicalModels
  TYPE STATIC
  SUBDIR_LIST
    ObjectStates
    TimeIntegrators
    ObjectModels
    InternalForceModel
  DEPENDS    
    imstkDataStructures
    imstkConstraints
    imstkGeometry
    imstkSolvers
    imstkMeshIO
    VegaFEM::massSpringSystem
    VegaFEM::corotationalLinearFEM
    VegaFEM::isotropicHyperelasticFEM
    VegaFEM::forceModel
    VegaFEM::stvk
    VegaFEM::graph
    VegaFEM::volumetricMesh
    PhysX
  )
 
#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
#if( ${PROJECT_NAME}_BUILD_TESTING )
#  add_subdirectory( Testing )
#endif()