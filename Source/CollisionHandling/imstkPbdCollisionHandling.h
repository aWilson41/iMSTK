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

#pragma once

#include "imstkCollisionHandling.h"
#include "imstkPbdConstraint.h"

#include <unordered_map>
#include <iostream>

namespace imstk
{
enum PbdContactCase
{
    None,
    Vertex,   // Mesh vertex either from Pbd or CollidingObject
    Edge,     // Mesh edge either from Pbd or CollidingObject
    Triangle, // Mesh triangle either from Pbd or CollidingObject
    Body,     // Body, usually a PD contact
    Primitive // Not a mesh, could be V,E,T, or PD additionally
};
static std::string
getContactCaseStr(PbdContactCase contactCase)
{
    switch (contactCase)
    {
    case PbdContactCase::Body:
        return "Body";
    case PbdContactCase::Edge:
        return "Edge";
    case PbdContactCase::Primitive:
        return "Primitive";
    case PbdContactCase::Triangle:
        return "Triangle";
    case PbdContactCase::Vertex:
        return "Vertex";
    default:
        return "None";
    }
    ;
}

///
/// \struct PbdCHFuncTableKey
///
/// \brief Used as a key in a function table to decide how to handle
/// resulting collision.
///
struct PbdCHTableKey
{
    PbdContactCase elemAType;
    PbdContactCase elemBType;
    //enum { IsPbdObjectB, IsCollidingObjectB } objectType; // \todo: edge cases
    //bool ccdA;
    //bool ccdB;

    friend std::ostream& operator<<(std::ostream& os, const PbdCHTableKey& dt);

    bool operator==(const PbdCHTableKey& other) const
    {
        return
            (elemAType == other.elemAType)
            && (elemBType == other.elemBType);
    }
};
} // namespace imstk

namespace std
{
///
/// \struct hash<imstk::PbdCHFuncTableKey>
///
/// \brief Gives hashing function for PbdCHFuncTableKey
/// complete unique/garunteed no collisions
///
template<>
struct hash<imstk::PbdCHTableKey>
{
    std::size_t operator()(const imstk::PbdCHTableKey& k) const
    {
        using std::size_t;
        using std::hash;

        // Base on the bit width of each value
        std::size_t v0 = static_cast<std::size_t>(k.elemAType); // first 2 bits
        std::size_t v1 = static_cast<std::size_t>(k.elemBType); // Next 2 bits
        return v0 ^ (v1 << 3);
    }
};
} // namespace std

namespace imstk
{
class PbdObject;
class PbdModel;
class PointSet;
class PointwiseMap;

struct MeshSide;
///
/// \class PbdCollisionHandling
///
/// \brief Implements PBD based collision handling. Given an input PbdObject
/// and CollisionData it creates & adds constraints in the PbdModel to be solved
/// in order to resolve the collision. This solve is ultimately implemented
/// in PbdModel, shared among all the collisions.
///
/// This supports PointDirection (PD) collision data as well as contacting feature
/// collision data (ie: EE, VT, VV, VE). The VV and VE are often redundant but
/// handled anyways for robustness. The PD is often reported for point contacts,
/// most common on primitive vs mesh collisions.
///
class PbdCollisionHandling : public CollisionHandling
{
public:
    enum class ObjType
    {
        PbdDeformable,
        PbdRigid,
        Colliding
    };
    struct CollisionSideData
    {
        CollisionSideData() = default;

        // Objects
        PbdObject* pbdObj       = nullptr;
        CollidingObject* colObj = nullptr;
        ObjType objType = ObjType::Colliding;

        PbdModel* model    = nullptr;
        double stiffness   = 0.0;
        Geometry* geometry = nullptr;
        PointSet* pointSet = nullptr;
        VecDataArray<double, 3>* vertices = nullptr;
        PointwiseMap* mapPtr = nullptr;
        AbstractDataArray* indicesPtr = nullptr;
        int bodyId = 0;

        Geometry* prevGeometry = nullptr;
    };
    ///
    /// \brief Packs the collision element together with the
    /// data it will need to process it (for swapping)
    ///
    struct ColElemSide
    {
        const CollisionElement* elem  = nullptr;
        const CollisionSideData* data = nullptr;
    };

    PbdCollisionHandling();
    ~PbdCollisionHandling() override;

    IMSTK_TYPE_NAME(PbdCollisionHandling)

    ///
    /// \brief Get/Set the restitution, which gives how much velocity is
    /// removed along the contact normals during contact
    /// @{
    double getRestitution() const { return m_restitution; }
    void setRestitution(const double restitution) { m_restitution = restitution; }
    /// @}

    ///
    /// \brief Get/Set the friction, which gives how much velocity is
    /// removed along the tangents during contact
    /// @{
    double getFriction() const { return m_friction; }
    void setFriction(const double friction) { m_friction = friction; }
    /// @}

    std::pair<PbdParticleId, Vec3d> getBodyAndContactPoint(const CollisionElement&  elem,
                                                           const CollisionSideData& data);

    ///
    /// \brief If geometry is a pointset returns cell array, otherwise returns nullptr
    ///
    AbstractDataArray* getCellArray(std::shared_ptr<Geometry> pointSet);

protected:
    std::array<PbdParticleId, 2> getEdge(const CollisionElement&  elem,
                                         const CollisionSideData& side);
    std::array<PbdParticleId, 3> getTriangle(const CollisionElement&  elem,
                                             const CollisionSideData& side);
    ///
    /// \brief getVertex takes slightly differing paths than the others, as the
    /// cell vertex directly refers to the vertex buffer, not an index buffer
    ///
    std::array<PbdParticleId, 1> getVertex(const CollisionElement&  elem,
                                           const CollisionSideData& side);

    ///
    /// \brief Creates a CollisionSideData struct from the provided object, this
    /// gives all the info needed to response to collision
    ///
    CollisionSideData getDataFromObject(std::shared_ptr<CollidingObject> obj);

    ///
    /// \brief Get the contact case from the collision element and data as
    /// additional context
    ///
    PbdContactCase getCaseFromElement(const ColElemSide& elem);

    ///
    /// \brief Add collision constraints based off contact data
    ///
    void handle(
        const std::vector<CollisionElement>& elementsA,
        const std::vector<CollisionElement>& elementsB) override;

    ///
    /// \brief Handle a single element
    ///
    void handleElementPair(ColElemSide sideA, ColElemSide sideB);

    // -----------------One-Way Rigid on X Cases-----------------
    virtual void V_Body(
        const ColElemSide& sideA,
        const ColElemSide& sideB);
    virtual void E_Body(
        const ColElemSide& sideA,
        const ColElemSide& sideB);
    virtual void T_Body(
        const ColElemSide& sideA,
        const ColElemSide& sideB);
    // ---------------Two-Way Rigid on Rigid Cases---------------
    virtual void Body_Body(
        const ColElemSide& sideA,
        const ColElemSide& sideB);

    // ----------DeformableMesh on DeformableMesh Cases----------
    virtual void V_T(
        const ColElemSide& sideA,
        const ColElemSide& sideB);
    virtual void E_E(
        const ColElemSide& sideA,
        const ColElemSide& sideB);
    virtual void V_E(
        const ColElemSide& sideA,
        const ColElemSide& sideB);
    virtual void V_V(
        const ColElemSide& sideA,
        const ColElemSide& sideB);

private:
    ///< Constraints used in the solver (user should not add too these as they are cleared)
    std::vector<PbdConstraint*> m_solverConstraints;

    double m_restitution = 0.0; ///< Coefficient of restitution (1.0 = perfect elastic, 0.0 = inelastic)
    double m_friction    = 0.0; ///< Coefficient of friction (1.0 = full frictional force, 0.0 = none)

protected:
    std::vector<PbdConstraint*> m_constraints; ///< Constraints users can add too

    std::unordered_map<PbdCHTableKey, std::function<void(
                                                        const ColElemSide& elemA, const ColElemSide& elemB)>> m_funcTable;
};
} // namespace imstk
