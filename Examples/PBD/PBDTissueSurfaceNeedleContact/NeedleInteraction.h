/*
** This file is part of the Interactive Medical Simulation Toolkit (iMSTK)
** iMSTK is distributed under the Apache License, Version 2.0.
** See accompanying NOTICE for details.
*/

#pragma once

#include "imstkPbdRigidObjectCollision.h"

#include "NeedlePbdCH.h"
#include "NeedleRigidBodyCH.h"

using namespace imstk;

///
/// \class NeedleInteraction
///
/// \brief Defines interaction between NeedleObject and PbdObject
///
class NeedleInteraction : public PbdRigidObjectCollision
{
public:
    NeedleInteraction(std::shared_ptr<PbdObject> tissueObj, std::shared_ptr<NeedleObject> needleObj) : PbdRigidObjectCollision(tissueObj, needleObj)
    {
        if (std::dynamic_pointer_cast<LineMesh>(needleObj->getCollidingGeometry()) == nullptr)
        {
            LOG(WARNING) << "NeedleInteraction only works with LineMesh collision geometry on NeedleObject";
        }

        auto needleRbdCH = std::make_shared<NeedleRigidBodyCH>();
        needleRbdCH->setInputRigidObjectA(needleObj);
        needleRbdCH->setInputCollidingObjectB(tissueObj);
        needleRbdCH->setInputCollisionData(getCollisionDetection()->getCollisionData());
        needleRbdCH->setBaumgarteStabilization(0.001);
        setCollisionHandlingB(needleRbdCH);

        auto needlePbdCH = std::make_shared<NeedlePbdCH>();
        needlePbdCH->setInputObjectA(tissueObj);
        needlePbdCH->setInputObjectB(needleObj);
        needlePbdCH->setInputCollisionData(getCollisionDetection()->getCollisionData());
        // These two can control compliance
        needlePbdCH->setDeformableStiffnessA(1.0);
        needlePbdCH->setDeformableStiffnessB(0.01);
        setCollisionHandlingA(needlePbdCH);
    }

    ~NeedleInteraction() override = default;
};