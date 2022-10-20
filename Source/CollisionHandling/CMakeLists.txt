#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------

set(H_FILES
  imstkCollisionHandling.h
  imstkLevelSetCH.h
  imstkPbdCollisionHandling.h
  imstkRigidBodyCH.h
  imstkSphCollisionHandling.h
  )
set(SRC_FILES
  imstkCollisionHandling.cpp
  imstkLevelSetCH.cpp
  imstkPbdCollisionHandling.cpp
  imstkRigidBodyCH.cpp
  imstkSphCollisionHandling.cpp
  )

if (iMSTK_USE_VegaFEM)
  list(APPEND H_FILES
    imstkBoneDrillingCH.h
    imstkPenaltyCH.h
    )
  list(APPEND SRC_FILES
    imstkBoneDrillingCH.cpp
    imstkPenaltyCH.cpp
    )
endif()

include(imstkAddLibrary)
imstk_add_library(CollisionHandling
  H_FILES
    ${H_FILES}
  CPP_FILES
    ${SRC_FILES}
  DEPENDS
    CollisionDetection
    SceneEntities
    Controllers
  )