include(imstkAddTest)

set(SRC_FILES
  imstkDataArrayTest.cpp
  imstkEventObjectTest.cpp
  imstkFactoryTest.cpp
  imstkMathTest.cpp
  imstkModuleTest.cpp
  imstkSequentialTaskGraphControllerTest.cpp
  imstkTaskGraphTest.cpp
  imstkVecDataArrayTest.cpp
  )

if (iMSTK_USE_TBB)
  list(APPEND SRC_FILES
    imstkTbbTaskGraphControllerTest.cpp
    )
endif()

imstk_add_test(Common
  CPP_FILES
    ${SRC_FILES}
    )