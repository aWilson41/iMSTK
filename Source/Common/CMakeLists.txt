#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------

list(APPEND Dependencies
  Eigen3::Eigen
  )

set(H_FILES
  imstkAbstractDataArray.h
  imstkColor.h
  imstkColorFunction.h
  imstkDataArray.h
  imstkEventObject.h
  imstkFactory.h
  imstkLogger.h
  imstkMacros.h
  imstkMath.h
  imstkModule.h
  imstkModuleDriver.h
  imstkNew.h
  imstkTypes.h
  imstkVecDataArray.h
  Parallel/imstkAtomicOperations.h
  Parallel/imstkParallelFor.h
  Parallel/imstkParallelReduce.h
  Parallel/imstkParallelUtils.h
  Parallel/imstkSpinLock.h
  Parallel/imstkThreadManager.h
  TaskGraph/imstkSequentialTaskGraphController.h
  TaskGraph/imstkTaskGraph.h
  TaskGraph/imstkTaskGraphController.h
  TaskGraph/imstkTaskGraphVizWriter.h
  TaskGraph/imstkTaskNode.h
  Utils/imstkTimer.h
  )
set(SRC_FILES
  imstkColor.cpp
  imstkModule.cpp
  imstkModuleDriver.cpp
  Parallel/imstkThreadManager.cpp
  TaskGraph/imstkSequentialTaskGraphController.cpp
  TaskGraph/imstkTaskGraph.cpp
  TaskGraph/imstkTaskGraphController.cpp
  TaskGraph/imstkTaskGraphVizWriter.cpp
  TaskGraph/imstkTaskNode.cpp
  Utils/imstkTimer.cpp
  )

if (iMSTK_USE_TBB)
  list(APPEND Dependencies
    TBB::tbb
    )
  list(APPEND H_FILES
    TaskGraph/imstkTbbTaskGraphController.h
    )
  list(APPEND SRC_FILES
    TaskGraph/imstkTbbTaskGraphController.cpp
    )
endif()

if (iMSTK_USE_G3LOG)
  list(APPEND Dependencies g3log::g3log)
  list(APPEND H_FILES imstkLoggerG3.h)
  list(APPEND SRC_FILES imstkLoggerG3.cpp)
else()
  list(APPEND H_FILES imstkLoggerSynchronous.h)
  list(APPEND SRC_FILES imstkLoggerSynchronous.cpp)
endif()

include(imstkAddLibrary)
imstk_add_library( Common
  H_FILES
    ${H_FILES}
  CPP_FILES
    ${SRC_FILES}
  DEPENDS
    ${Dependencies}
  )

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
if( ${PROJECT_NAME}_BUILD_TESTING )
  add_subdirectory(Testing)
endif()
