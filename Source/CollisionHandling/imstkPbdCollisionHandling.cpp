/*=========================================================================

   Library: iMSTK

   Copyright (c) Kitware, Inc. & Center for Modeling, Simulation,
   & Imaging in Medicine, Rensselaer Polytechnic Institute.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0.txt

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

=========================================================================*/

#include "imstkPbdCollisionHandling.h"
#include "imstkPbdContactConstraint.h"
#include "imstkPbdEdgeEdgeConstraint.h"
#include "imstkPbdEdgeEdgeCCDConstraint.h"
#include "imstkPbdModel.h"
#include "imstkPbdObject.h"
#include "imstkPbdPointEdgeConstraint.h"
#include "imstkPbdPointPointConstraint.h"
#include "imstkPbdPointTriangleConstraint.h"
#include "imstkPbdSolver.h"
#include "imstkPointwiseMap.h"
#include "imstkLineMesh.h"
#include "imstkSurfaceMesh.h"

namespace imstk
{
#define REGISTER_CASE(case0, case1, func)                                                      \
    m_funcTable[{ case0, case1 }] = [this](const ColElemSide& sideA, const ColElemSide& sideB) \
                                    { func(sideA, sideB); }

std::pair<PbdParticleId, Vec3d>
PbdCollisionHandling::getBodyAndContactPoint(const CollisionElement& elem, const CollisionSideData& data)
{
    if (elem.m_type == CollisionElementType::PointDirection)
    {
        return { { data.bodyId, 0 }, elem.m_element.m_PointDirectionElement.pt };
    }
    else
    {
        const PbdParticleId& ptBv = getVertex(elem, data)[0];
        return { { data.bodyId, 0 }, (*data.vertices)[ptBv.second] };
    }
}

AbstractDataArray*
PbdCollisionHandling::getCellArray(std::shared_ptr<Geometry> pointSet)
{
    if (auto cellMesh = std::dynamic_pointer_cast<AbstractCellMesh>(pointSet))
    {
        return cellMesh->getIndices().get();
    }
    else
    {
        return nullptr;
    }
}

///
/// \brief Gets triangle, edge, or vertex from the mesh given the CollisionElement
/// \tparam Array type for vertex index storage, ex: VecDataArray<int, 4>, VecDataArray<int, 3>
/// \tparam Cell type the index array is representing
///
template<typename ArrType, int cellType>
static std::array<PbdParticleId, ArrType::NumComponents>
getElementVertIds(const CollisionElement& elem, const PbdCollisionHandling::CollisionSideData& side)
{
    // Note: The unrolling of this functions loops could be important too performance
    typename ArrType::ValueType cell;
    if (elem.m_type == CollisionElementType::CellIndex && elem.m_element.m_CellIndexElement.cellType == cellType)
    {
        // If one index it refers to the cell
        if (elem.m_element.m_CellIndexElement.idCount == 1)
        {
            cell = (*dynamic_cast<ArrType*>(side.indicesPtr))[elem.m_element.m_CellIndexElement.ids[0]];
        }
        else if (elem.m_element.m_CellIndexElement.idCount == ArrType::NumComponents)
        {
            for (int i = 0; i < ArrType::NumComponents; i++)
            {
                cell[i] = elem.m_element.m_CellIndexElement.ids[i];
            }
        }
    }

    std::array<PbdParticleId, ArrType::NumComponents> results;
    if (cell[0] != -1)
    {
        if (side.mapPtr != nullptr)
        {
            for (int i = 0; i < ArrType::NumComponents; i++)
            {
                cell[i] = side.mapPtr->getParentVertexId(cell[i]);
            }
        }
        for (int i = 0; i < ArrType::NumComponents; i++)
        {
            int vid = cell[i];
            if (side.bodyId == 0)
            {
                vid = side.model->addVirtualParticle((*side.vertices)[vid], 0.0).second;
            }
            results[i] = { side.bodyId, vid };
        }
    }
    return results;
}

std::array<PbdParticleId, 2>
PbdCollisionHandling::getEdge(const CollisionElement& elem, const CollisionSideData& side)
{
    return getElementVertIds<VecDataArray<int, 2>, IMSTK_EDGE>(elem, side);
}

std::array<PbdParticleId, 3>
PbdCollisionHandling::getTriangle(const CollisionElement& elem, const CollisionSideData& side)
{
    return getElementVertIds<VecDataArray<int, 3>, IMSTK_TRIANGLE>(elem, side);
}

std::array<PbdParticleId, 1>
PbdCollisionHandling::getVertex(const CollisionElement& elem, const CollisionSideData& side)
{
    std::array<PbdParticleId, 1> results;
    int                          ptId = -1;
    if (elem.m_type == CollisionElementType::CellIndex && elem.m_element.m_CellIndexElement.cellType == IMSTK_VERTEX)
    {
        ptId = elem.m_element.m_CellIndexElement.ids[0];
    }
    else if (elem.m_type == CollisionElementType::PointIndexDirection)
    {
        ptId = elem.m_element.m_PointIndexDirectionElement.ptIndex;
    }
    if (ptId != -1)
    {
        if (side.mapPtr != nullptr)
        {
            ptId = side.mapPtr->getParentVertexId(ptId);
        }
        if (side.bodyId == 0)
        {
            ptId = side.model->addVirtualParticle((*side.vertices)[ptId], 0.0).second;
        }
        results[0] = { side.bodyId, ptId };
    }
    return results;
}

std::ostream&
operator<<(std::ostream& os, const PbdCHTableKey& key)
{
    os << getContactCaseStr(key.elemAType) << " vs " << getContactCaseStr(key.elemBType);
    return os;
}

PbdCollisionHandling::PbdCollisionHandling()
{
    REGISTER_CASE(PbdContactCase::Vertex, PbdContactCase::Vertex, V_V);
    REGISTER_CASE(PbdContactCase::Vertex, PbdContactCase::Edge, V_E);
    REGISTER_CASE(PbdContactCase::Edge, PbdContactCase::Edge, E_E);
    REGISTER_CASE(PbdContactCase::Vertex, PbdContactCase::Triangle, V_T);

    REGISTER_CASE(PbdContactCase::Triangle, PbdContactCase::Body, T_Body);
    REGISTER_CASE(PbdContactCase::Edge, PbdContactCase::Body, E_Body);
    REGISTER_CASE(PbdContactCase::Vertex, PbdContactCase::Body, V_Body);

    // If swap occurs the colliding object could be on the LHS causing issues
    REGISTER_CASE(PbdContactCase::Primitive, PbdContactCase::Triangle, V_T);
    REGISTER_CASE(PbdContactCase::Primitive, PbdContactCase::Edge, V_E);
    REGISTER_CASE(PbdContactCase::Primitive, PbdContactCase::Vertex, V_V);

    // One way point direction resolution
    REGISTER_CASE(PbdContactCase::Vertex, PbdContactCase::None, V_V);
}

PbdCollisionHandling::~PbdCollisionHandling()
{
    for (const auto ptr: m_constraints)
    {
        delete ptr;
    }
}

PbdCollisionHandling::CollisionSideData
PbdCollisionHandling::getDataFromObject(std::shared_ptr<CollidingObject> obj)
{
    // Pack info into struct, gives some contextual hints as well
    CollisionSideData side;
    side.pbdObj    = std::dynamic_pointer_cast<PbdObject>(obj).get(); // Garunteed
    side.colObj    = obj.get();
    side.objType   = ObjType::Colliding;
    side.model     = nullptr;
    side.mapPtr    = nullptr;
    side.geometry  = side.colObj->getCollidingGeometry().get();
    side.bodyId    = 0;
    side.stiffness = 0.0;

    if (side.pbdObj != nullptr)
    {
        if (side.pbdObj->getPbdBody()->bodyType == PbdBody::Type::RIGID)
        {
            side.objType = ObjType::PbdRigid;
        }
        else
        {
            side.objType = ObjType::PbdDeformable;
        }
        side.model     = side.pbdObj->getPbdModel().get();
        side.stiffness = side.model->getConfig()->m_contactStiffness;
        // If a physics geometry is provided always use that because
        // Either:
        //  A.) Physics geometry == Collision Geometry
        //  B.) A PointwiseMap is used and map should refer us back to physics geometry
        side.geometry = side.pbdObj->getPhysicsGeometry().get();
        side.bodyId   = side.pbdObj->getPbdBody()->bodyHandle;

        auto map = std::dynamic_pointer_cast<PointwiseMap>(side.pbdObj->getPhysicsToCollidingMap());
        side.mapPtr = map.get();
    }
    side.pointSet = dynamic_cast<PointSet*>(side.geometry);
    side.vertices = (side.pointSet == nullptr) ? nullptr : side.pointSet->getVertexPositions().get();
    if (side.objType == ObjType::PbdRigid)
    {
        if (auto pointSet = std::dynamic_pointer_cast<PointSet>(side.colObj->getCollidingGeometry()))
        {
            side.vertices = pointSet->getVertexPositions().get();
        }
        //side.vertices = side.colObj->getCollidingGeometry()
    }
    side.indicesPtr = getCellArray(side.colObj->getCollidingGeometry());

    return side;
}

PbdContactCase
PbdCollisionHandling::getCaseFromElement(const ColElemSide& elem)
{
    if (elem.data->objType == ObjType::PbdRigid)
    {
        return PbdContactCase::Body;
    }
    else
    {
        if (elem.elem->m_type == CollisionElementType::PointDirection)
        {
            return PbdContactCase::Primitive;
        }
        else if (elem.elem->m_type == CollisionElementType::CellIndex)
        {
            // 0 - vertex, 1 - edge, 2 - triangle
            return static_cast<PbdContactCase>(
                elem.elem->m_element.m_CellIndexElement.cellType + 1);
        }
        else if (elem.elem->m_type == CollisionElementType::CellVertex)
        {
            return static_cast<PbdContactCase>(
                elem.elem->m_element.m_CellVertexElement.size);
        }
        // If not rigid, this must be asking to resolve the vertex
        else if (elem.elem->m_type == CollisionElementType::PointIndexDirection)
        {
            return PbdContactCase::Vertex;
        }
    }
    return PbdContactCase::None;
}

void
PbdCollisionHandling::handle(
    const std::vector<CollisionElement>& elementsA,
    const std::vector<CollisionElement>& elementsB)
{
    // Remove constraints without actually clearing
    // There could be a large variance in constraints/contacts count,
    // 10s to 100s of constraints changing frequently over few frames.
    m_constraints.resize(0);

    // Break early if no collision elements
    if (elementsA.size() == 0 && elementsB.size() == 0)
    {
        return;
    }

    // Pack all the data needed for a particular side into a struct so we can
    // swap it with the contact & pass it around
    CollisionSideData dataSideA = getDataFromObject(getInputObjectA());
    CollisionSideData dataSideB = getDataFromObject(getInputObjectB());

    // Share the model with both sides even if B is not pbd, which makes it easier
    // to acquire the model without knowing which object is pbd
    if (dataSideB.model == nullptr)
    {
        dataSideB.model = dataSideA.model;
    }

    // If obj B is also pbd simulated, make sure they share the same model
    CHECK(dataSideB.pbdObj == nullptr || dataSideA.model == dataSideB.model) <<
        "PbdCollisionHandling input objects must share PbdModel";

    // For CCD (store if available)
    dataSideA.prevGeometry = m_colData->prevGeomA.get();
    dataSideB.prevGeometry = m_colData->prevGeomB.get();

    if (elementsA.size() == elementsB.size())
    {
        // Deal with two way contacts
        for (size_t i = 0; i < elementsA.size(); i++)
        {
            handleElementPair(
                { &elementsA[i], &dataSideA },
                { &elementsB[i], &dataSideB });
        }
    }
    else
    {
        // Deal with one way contacts (only one side is needed)
        for (size_t i = 0; i < elementsA.size(); i++)
        {
            handleElementPair(
                { &elementsA[i], &dataSideA },
                { nullptr, nullptr });
        }
        for (size_t i = 0; i < elementsB.size(); i++)
        {
            handleElementPair(
                { &elementsB[i], &dataSideB },
                { nullptr, nullptr });
        }
    }

    // Copy constraints
    for (int i = 0; i < m_solverConstraints.size(); i++)
    {
        delete m_solverConstraints[i];
    }
    m_solverConstraints.resize(0);
    m_solverConstraints.reserve(m_constraints.size());

    for (size_t i = 0; i < m_constraints.size(); i++)
    {
        m_solverConstraints.push_back(static_cast<PbdConstraint*>(m_constraints[i]));
    }

    if (m_solverConstraints.size() == 0)
    {
        return;
    }

    // ObjA garunteed to be PbdObject
    auto pbdObjectA = std::dynamic_pointer_cast<PbdObject>(getInputObjectA());
    pbdObjectA->getPbdModel()->getCollisionSolver()->addConstraints(&m_constraints);
}

void
PbdCollisionHandling::handleElementPair(ColElemSide sideA, ColElemSide sideB)
{
    PbdCHTableKey key;
    key.elemAType = getCaseFromElement(sideA);
    key.elemBType = PbdContactCase::None;

    // Data for sideB may not be present if handling one-way
    if (sideB.data != nullptr)
    {
        key.elemBType = getCaseFromElement(sideB);

        // Avoid a couple cases by swapping
        // Only V_T, no TV
        // Only V_E, no EV
        // Always Body on the right
        // Always Primitive on left
        if ((key.elemAType == PbdContactCase::Triangle && key.elemBType == PbdContactCase::Vertex)
            || (key.elemAType == PbdContactCase::Edge && key.elemBType == PbdContactCase::Vertex)
            || (key.elemAType == PbdContactCase::Body && key.elemBType != PbdContactCase::Body)
            || (key.elemAType != PbdContactCase::Primitive && key.elemBType == PbdContactCase::Primitive))
        {
            std::swap(key.elemAType, key.elemBType);
            std::swap(sideA, sideB);
        }
    }

    auto iter = m_funcTable.find(key);
    if (iter != m_funcTable.end())
    {
        iter->second(sideA, sideB);
    }
    else
    {
        // CH's may not handle all CollisionElements types
        LOG(INFO) << "Could not find handling case " << key;
    }
}

void
PbdCollisionHandling::V_Body(const ColElemSide& sideA, const ColElemSide& sideB)
{
    const PbdParticleId&                   ptsA = getVertex(*sideA.elem, *sideA.data)[0];
    const std::pair<PbdParticleId, Vec3d>& ptBAndContact = getBodyAndContactPoint(*sideB.elem, *sideB.data);

    PbdVertexToBodyConstraint* constraint = new PbdVertexToBodyConstraint();
    constraint->initConstraint(sideA.data->model->getBodies(),
        ptBAndContact.first,
        ptBAndContact.second,
        ptsA,
        sideA.data->stiffness); // stiffness choice?
    constraint->setFriction(m_friction);
    constraint->setRestitution(m_restitution);
    m_constraints.push_back(constraint);
}

void
PbdCollisionHandling::E_Body(const ColElemSide& sideA, const ColElemSide& sideB)
{
    std::array<PbdParticleId, 2>           ptsA = getEdge(*sideA.elem, *sideA.data);
    const std::pair<PbdParticleId, Vec3d>& ptBAndContact = getBodyAndContactPoint(*sideB.elem, *sideB.data);

    PbdEdgeToBodyConstraint* constraint = new PbdEdgeToBodyConstraint();
    constraint->initConstraint(sideA.data->model->getBodies(),
        ptBAndContact.first,
        ptBAndContact.second,
        ptsA[0], ptsA[1],
        sideA.data->stiffness); // stiffness choice?
    constraint->setFriction(m_friction);
    constraint->setRestitution(m_restitution);
    m_constraints.push_back(constraint);
}

void
PbdCollisionHandling::T_Body(const ColElemSide& sideA, const ColElemSide& sideB)
{
    std::array<PbdParticleId, 3>           ptsA = getTriangle(*sideA.elem, *sideA.data);
    const std::pair<PbdParticleId, Vec3d>& ptBAndContact = getBodyAndContactPoint(*sideB.elem, *sideB.data);

    PbdTriangleToBodyConstraint* constraint = new PbdTriangleToBodyConstraint();
    constraint->initConstraint(sideA.data->model->getBodies(),
        ptBAndContact.first,
        ptBAndContact.second,
        ptsA[0], ptsA[1], ptsA[2],
        sideA.data->stiffness); // stiffness choice?
    constraint->setFriction(m_friction);
    constraint->setRestitution(m_restitution);
    m_constraints.push_back(constraint);
}

void
PbdCollisionHandling::Body_Body(const ColElemSide& sideA, const ColElemSide& sideB)
{
}

void
PbdCollisionHandling::V_T(const ColElemSide& sideA, const ColElemSide& sideB)
{
    const PbdParticleId&         ptA  = getVertex(*sideA.elem, *sideA.data)[0];
    std::array<PbdParticleId, 3> ptsB = getTriangle(*sideB.elem, *sideB.data);

    PbdPointTriangleConstraint* constraint = new PbdPointTriangleConstraint();
    constraint->initConstraint(ptA, ptsB[0], ptsB[1], ptsB[2],
        sideA.data->stiffness, sideB.data->stiffness);
    constraint->setFriction(m_friction);
    constraint->setRestitution(m_restitution);
    m_constraints.push_back(constraint);
}

void
PbdCollisionHandling::E_E(const ColElemSide& sideA, const ColElemSide& sideB)
{
    std::array<PbdParticleId, 2> ptsA = getEdge(*sideA.elem, *sideA.data);
    std::array<PbdParticleId, 2> ptsB = getEdge(*sideB.elem, *sideB.data);

    //if (sideA.elem->m_ccdData && sideB.elem->m_ccdData
    //    && prevSideA.isValid() && prevSideB.isValid())
    //{
    //    PointSet* prevPointSetA = dynamic_cast<PointSet*>(sideA.data->prevGeometry);
    //    PointSet* prevPointSetB = dynamic_cast<PointSet*>(sideB.data->prevGeometry);

    //    // Create prev timestep MeshSide structures if the previous geometries are available.
    //    MeshSide prevSideA = MeshSide::create(prevPointSetA.get(), pbdObjectA->getPhysicsToCollidingMap().get());
    //    MeshSide prevSideB = MeshSide::create(prevPointSetB.get(), (pbdObjectB == nullptr) ? nullptr : pbdObjectB->getPhysicsToCollidingMap().get());

    //    std::array<VertexMassPair, 2> prevVertexMassA = getEdge(colElemA.m_element.m_CellIndexElement, prevSideA);
    //    std::array<VertexMassPair, 2> prevVertexMassB = getEdge(colElemB.m_element.m_CellIndexElement, prevSideB);

    //    PbdEdgeEdgeCCDConstraint* constraint = new PbdEdgeEdgeCCDConstraint();
    //    constraint->initConstraint(
    //        prev_ptA1, prev_ptA2, prev_ptB1, prev_ptB2,
    //        ptA1, ptA2, ptB1, ptB2,
    //        stiffnessA, stiffnessB);
    //    m_constraints.push_back(constraint);
    //}
    //else
    {
        PbdEdgeEdgeConstraint* constraint = new PbdEdgeEdgeConstraint();
        constraint->initConstraint(ptsA[0], ptsA[1], ptsB[0], ptsB[1],
            sideA.data->stiffness, sideB.data->stiffness);
        constraint->setFriction(m_friction);
        constraint->setRestitution(m_restitution);
        m_constraints.push_back(constraint);
    }
}

void
PbdCollisionHandling::V_E(const ColElemSide& sideA, const ColElemSide& sideB)
{
    const PbdParticleId&         ptA  = getVertex(*sideA.elem, *sideA.data)[0];
    std::array<PbdParticleId, 2> ptsB = getEdge(*sideB.elem, *sideB.data);

    PbdPointEdgeConstraint* constraint = new PbdPointEdgeConstraint();
    constraint->initConstraint(ptA, ptsB[0], ptsB[1],
        sideA.data->stiffness, sideB.data->stiffness);
    constraint->setFriction(m_friction);
    constraint->setRestitution(m_restitution);
    m_constraints.push_back(constraint);
}

void
PbdCollisionHandling::V_V(const ColElemSide& sideA, const ColElemSide& sideB)
{
    // One special case with one-way
    const PbdParticleId& ptA = getVertex(*sideA.elem, *sideA.data)[0];
    PbdParticleId        ptB;
    if (sideB.data == nullptr && sideA.elem->m_type == CollisionElementType::PointIndexDirection)
    {
        const Vec3d  pos         = sideA.data->model->getBodies().getPosition(ptA);
        const Vec3d  dir         = sideA.elem->m_element.m_PointIndexDirectionElement.dir;
        const double depth       = sideA.elem->m_element.m_PointIndexDirectionElement.penetrationDepth;
        const Vec3d  resolvedPos = pos + dir * depth;
        ptB = sideA.data->model->addVirtualParticle(resolvedPos, 0.0);
    }
    else
    {
        ptB = getVertex(*sideB.elem, *sideB.data)[0];
    }

    PbdPointPointConstraint* constraint = new PbdPointPointConstraint();
    constraint->initConstraint(ptA, ptB,
        sideA.data->stiffness, (sideB.data == nullptr) ? 0.0 : sideB.data->stiffness);
    constraint->setFriction(m_friction);
    constraint->setRestitution(m_restitution);
    m_constraints.push_back(constraint);
}
} // namespace imstk