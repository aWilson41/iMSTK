#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------

list(APPEND Dependencies
  DynamicalModels
  MeshIO
  VegaFEM::massSpringSystem
  VegaFEM::corotationalLinearFEM
  VegaFEM::isotropicHyperelasticFEM
  VegaFEM::forceModel
  VegaFEM::stvk
  VegaFEM::graph
  VegaFEM::volumetricMesh
  )

set(H_FILES
  InternalForceModel/imstkCorotationalFemForceModel.h
  InternalForceModel/imstkInternalForceModel.h
  InternalForceModel/imstkInternalForceModelTypes.h
  InternalForceModel/imstkIsotropicHyperelasticFeForceModel.h
  InternalForceModel/imstkLinearFemForceModel.h
  InternalForceModel/imstkStVKForceModel.h
  imstkFemDeformableBodyModel.h
  imstkVectorizedState.h
  TimeIntegrators/imstkBackwardEuler.h
  TimeIntegrators/imstkNewmarkBeta.h
  TimeIntegrators/imstkTimeIntegrator.h
  )

set(SRC_FILES
  InternalForceModel/imstkCorotationalFemForceModel.cpp
  InternalForceModel/imstkInternalForceModel.cpp
  InternalForceModel/imstkIsotropicHyperelasticFeForceModel.cpp
  InternalForceModel/imstkLinearFemForceModel.cpp
  InternalForceModel/imstkStVKForceModel.cpp
  imstkFemDeformableBodyModel.cpp
  imstkVectorizedState.cpp
  TimeIntegrators/imstkBackwardEuler.cpp
  )

# If using model reduction link to these libs
if (iMSTK_USE_MODEL_REDUCTION)
  list(APPEND Dependencies
    VegaFEM::reducedStvk  
    VegaFEM::reducedElasticForceModel
    )

  list(APPEND H_FILES imstkReducedStVKBodyModel.h)
  list(APPEND SRC_FILES imstkReducedStVKBodyModel.cpp)
endif()

include(imstkAddLibrary)
imstk_add_library(DynamicalModelsVegaFEM
  H_FILES
    ${H_FILES}
  CPP_FILES
    ${SRC_FILES}
  DEPENDS
    ${Dependencies}
  )
 
#-----------------------------------------------------------------------------
# Testing and benchmarking
#-----------------------------------------------------------------------------
# if( ${PROJECT_NAME}_BUILD_TESTING )
#   add_subdirectory(Testing)
# endif()

# if( ${PROJECT_NAME}_BUILD_BENCHMARK )
#   add_subdirectory(Benchmarking)
# endif()
