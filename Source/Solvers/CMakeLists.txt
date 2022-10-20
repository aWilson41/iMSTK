#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------

include(imstkAddLibrary)

set(H_FILES
  imstkPbdSolver.h
  imstkProjectedGaussSeidelSolver.h
  imstkSolverBase.h
  )

set(SRC_FILES
  imstkPbdSolver.cpp
  )

if (NOT iMSTK_BUILD_FOR_ANDROID)
  list(APPEND H_FILES
    imstkConjugateGradient.h
    imstkDirectLinearSolver.h
    imstkGaussSeidel.h
    imstkIterativeLinearSolver.h
    imstkJacobi.h
    imstkLinearSolver.h
    imstkNewtonSolver.h
    imstkNonLinearSolver.h
    imstkNonLinearSystem.h
    imstkSOR.h
    )
  list(APPEND SRC_FILES
    imstkConjugateGradient.cpp
    imstkDirectLinearSolver.cpp
    imstkGaussSeidel.cpp
    imstkIterativeLinearSolver.cpp
    imstkJacobi.cpp
    imstkNewtonSolver.cpp
    imstkNonLinearSolver.cpp
    imstkNonLinearSystem.cpp
    imstkSOR.cpp
    )
endif()

imstk_add_library( Solvers
  H_FILES
    ${H_FILES}
  CPP_FILES
    ${SRC_FILES}
  DEPENDS
    Constraints
  )

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
if( ${PROJECT_NAME}_BUILD_TESTING )
  add_subdirectory(Testing)
endif()
