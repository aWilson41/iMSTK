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

#include "imstkLineMesh.h"
#include "imstkLogger.h"
#include "imstkParallelUtils.h"
#include "imstkPbdAreaConstraint.h"
#include "imstkPbdBendConstraint.h"
#include "imstkPbdConstantDensityConstraint.h"
#include "imstkPbdConstraint.h"
#include "imstkPbdDihedralConstraint.h"
#include "imstkPbdDistanceConstraint.h"
#include "imstkPbdFemConstraint.h"
#include "imstkPbdFemTetConstraint.h"
#include "imstkPbdVolumeConstraint.h"
#include "imstkPointSet.h"
#include "imstkSurfaceMesh.h"
#include "imstkTetrahedralMesh.h"
#include "imstkPbdConstraintContainer.h"

#include <set>

namespace imstk
{
///
/// \struct PbdConstraintFunctor
///
/// \brief PbdConstraintFunctor take input geometry and produce constraints.
/// It exists to allow extensible constraint generation. Any number of constraints
/// could be generated by a functor for a single geometry with any number of parameters
///
///
struct PbdConstraintFunctor
{
    public:
        PbdConstraintFunctor() = default;
        virtual ~PbdConstraintFunctor() = default;

        ///
        /// \brief Appends a set of constraint to the container given a geometry & body
        ///
        virtual void operator()(PbdConstraintContainer& constraints) = 0;

        ///
        /// \brief Appends a set of constraint to the container given a geometry
        /// and a set of newly inserted vertices
        /// This is for dealing with topology diffs during runtime
        ///
        virtual void addConstraints(
            PbdConstraintContainer&                     imstkNotUsed(constraints),
            std::shared_ptr<std::unordered_set<size_t>> imstkNotUsed(vertices)) { }
};

///
/// \class PbdConstraintFunctorLambda
///
/// \brief PbdConstraintFunctor that can be forwarded out to a function pointer
///
struct PbdConstraintFunctorLambda : public PbdConstraintFunctor
{
    public:
        PbdConstraintFunctorLambda() = default;
        PbdConstraintFunctorLambda(std::function<void(PbdConstraintContainer&)> f) :
            m_func(f) { }

        void operator()(PbdConstraintContainer& constraints) override
        {
            m_func(constraints);
        }

        std::function<void(PbdConstraintContainer&)> m_func = nullptr;
};

///
/// \class PbdBodyConstraintFunctor
///
/// \brief Abstract PbdConstraintFunctor that assumes usage of a singular body & geometry
///
struct PbdBodyConstraintFunctor : public PbdConstraintFunctor
{
    public:
        void setGeometry(std::shared_ptr<PointSet> geom)
        {
            m_geom = geom;
        }

        void setBodyIndex(const int bodyIndex) { m_bodyIndex = bodyIndex; }

        int m_bodyIndex = 0;
        std::shared_ptr<PointSet> m_geom = nullptr;
};

///
/// \struct PbdDistanceConstraintFunctor
///
/// \brief PbdDistanceConstraintFunctor generates constraints between the edges of the input
/// TetrahedralMesh, SurfaceMesh, or LineMesh
///
struct PbdDistanceConstraintFunctor : public PbdBodyConstraintFunctor
{
    public:
        PbdDistanceConstraintFunctor() = default;
        ~PbdDistanceConstraintFunctor() override = default;

        ///
        /// \brief Create the distance constraint
        ///
        virtual std::shared_ptr<PbdDistanceConstraint> makeDistConstraint(
            const VecDataArray<double, 3>& vertices,
            size_t i1, size_t i2)
        {
            auto constraint = std::make_shared<PbdDistanceConstraint>();
            constraint->initConstraint(vertices[i1], vertices[i2],
                { m_bodyIndex, i1 }, { m_bodyIndex, i2 }, m_stiffness);
            return constraint;
        }

        void operator()(PbdConstraintContainer& constraints) override
        {
            std::shared_ptr<VecDataArray<double, 3>> verticesPtr = m_geom->getVertexPositions();
            const VecDataArray<double, 3>&           vertices    = *verticesPtr;

            auto addDistConstraint = [&](
                std::vector<std::vector<bool>>& E, size_t i1, size_t i2)
                                     {
                                         // Make sure i1 is always smaller than i2 for duplicate testing
                                         if (i1 > i2)
                                         {
                                             std::swap(i1, i2);
                                         }
                                         if (E[i1][i2])
                                         {
                                             E[i1][i2] = 0;

                                             auto c = makeDistConstraint(vertices, i1, i2);
                                             constraints.addConstraint(c);
                                         }
                                     };

            if (auto tetMesh = std::dynamic_pointer_cast<TetrahedralMesh>(m_geom))
            {
                std::shared_ptr<VecDataArray<int, 4>> elementsPtr = tetMesh->getIndices();
                const VecDataArray<int, 4>&           elements    = *elementsPtr;
                const auto                            nV = tetMesh->getNumVertices();
                std::vector<std::vector<bool>>        E(nV, std::vector<bool>(nV, 1)); // To test for duplicates

                for (int k = 0; k < elements.size(); ++k)
                {
                    auto& tet = elements[k];
                    addDistConstraint(E, tet[0], tet[1]);
                    addDistConstraint(E, tet[0], tet[2]);
                    addDistConstraint(E, tet[0], tet[3]);
                    addDistConstraint(E, tet[1], tet[2]);
                    addDistConstraint(E, tet[1], tet[3]);
                    addDistConstraint(E, tet[2], tet[3]);
                }
            }
            else if (auto triMesh = std::dynamic_pointer_cast<SurfaceMesh>(m_geom))
            {
                std::shared_ptr<VecDataArray<int, 3>> elementsPtr = triMesh->getIndices();
                const VecDataArray<int, 3>&           elements    = *elementsPtr;
                const auto                            nV = triMesh->getNumVertices();
                std::vector<std::vector<bool>>        E(nV, std::vector<bool>(nV, 1)); // To test for duplicates

                for (int k = 0; k < elements.size(); ++k)
                {
                    auto& tri = elements[k];
                    addDistConstraint(E, tri[0], tri[1]);
                    addDistConstraint(E, tri[0], tri[2]);
                    addDistConstraint(E, tri[1], tri[2]);
                }
            }
            else if (auto lineMesh = std::dynamic_pointer_cast<LineMesh>(m_geom))
            {
                std::shared_ptr<VecDataArray<int, 2>> elementsPtr = lineMesh->getIndices();
                const VecDataArray<int, 2>&           elements    = *elementsPtr;
                const auto&                           nV = lineMesh->getNumVertices();
                std::vector<std::vector<bool>>        E(nV, std::vector<bool>(nV, 1)); // To test for duplicates

                for (int k = 0; k < elements.size(); k++)
                {
                    auto& seg = elements[k];
                    addDistConstraint(E, seg[0], seg[1]);
                }
            }
            else
            {
                LOG(WARNING) << "PbdDistanceConstraint can only be generated with a TetrahedralMesh, SurfaceMesh, or LineMesh";
            }
        }

        void addConstraints(PbdConstraintContainer&                     constraints,
                            std::shared_ptr<std::unordered_set<size_t>> vertices) override
        {
            // Check for correct mesh type
            CHECK(std::dynamic_pointer_cast<SurfaceMesh>(m_geom) != nullptr)
                << "Add element constraints does not support current mesh type.";

            auto                                     triMesh = std::dynamic_pointer_cast<SurfaceMesh>(m_geom);
            std::shared_ptr<VecDataArray<double, 3>> initVerticesPtr = m_geom->getInitialVertexPositions();
            std::shared_ptr<VecDataArray<int, 3>>    elementsPtr     = triMesh->getIndices();

            const VecDataArray<double, 3>& initVertices = *initVerticesPtr;
            const VecDataArray<int, 3>&    elements     = *elementsPtr;

            // Build vertex to tri map
            std::vector<std::vector<size_t>> vertexToTriMap(triMesh->getNumVertices());
            for (int k = 0; k < elements.size(); k++)
            {
                const Vec3i& tri = elements[k];
                vertexToTriMap[tri[0]].push_back(k);
                vertexToTriMap[tri[1]].push_back(k);
                vertexToTriMap[tri[2]].push_back(k);
            }
            for (std::vector<size_t>& faceIds : vertexToTriMap)
            {
                std::sort(faceIds.begin(), faceIds.end());
            }

            std::set<std::pair<size_t, size_t>> distanceSet;
            for (const size_t& vertIdx : *vertices)
            {
                for (const size_t& triIdx : vertexToTriMap[vertIdx])
                {
                    const Vec3i& tri = elements[triIdx];
                    size_t       i1  = 0;
                    size_t       i2  = 0;
                    for (size_t i = 0; i < 3; i++)
                    {
                        if (tri[i] == vertIdx)
                        {
                            i1 = tri[(i + 1) % 3];
                            i2 = tri[(i + 2) % 3];
                            break;
                        }
                    }
                    auto pair1 = std::make_pair(std::min(vertIdx, i1), std::max(vertIdx, i1));
                    auto pair2 = std::make_pair(std::min(vertIdx, i2), std::max(vertIdx, i2));
                    distanceSet.insert(pair1);
                    distanceSet.insert(pair2);
                }
            }

            constraints.reserve(constraints.getConstraints().size() + distanceSet.size());
            for (auto& c : distanceSet)
            {
                auto constraint = makeDistConstraint(initVertices, c.first, c.second);
                constraints.addConstraint(constraint);
            }
        }

        ///
        /// \brief Get/Set the stiffness, how hard the constraint is
        ///@{
        void setStiffness(const double stiffness) { m_stiffness = stiffness; }
        double getStiffness() const { return m_stiffness; }
    ///@}

    protected:
        double m_stiffness = 0.0;
};

///
/// \struct PbdFemTetConstraintFunctor
///
/// \brief PbdFemTetConstraintFunctor generates constraints per cell of a TetrahedralMesh
///
struct PbdFemTetConstraintFunctor : public PbdBodyConstraintFunctor
{
    public:
        PbdFemTetConstraintFunctor() = default;
        ~PbdFemTetConstraintFunctor() override = default;

        void operator()(PbdConstraintContainer& constraints) override
        {
            // Check for correct mesh type
            CHECK(std::dynamic_pointer_cast<TetrahedralMesh>(m_geom) != nullptr)
                << "PbdFemTetConstraint can only be generated with a TetrahedralMesh";

            // Create constraints
            auto                                     tetMesh     = std::dynamic_pointer_cast<TetrahedralMesh>(m_geom);
            std::shared_ptr<VecDataArray<double, 3>> verticesPtr = m_geom->getVertexPositions();
            const VecDataArray<double, 3>&           vertices    = *verticesPtr;
            std::shared_ptr<VecDataArray<int, 4>>    elementsPtr = tetMesh->getIndices();
            const VecDataArray<int, 4>&              elements    = *elementsPtr;

            ParallelUtils::parallelFor(elements.size(),
                [&](const size_t k)
                {
                    const Vec4i& tet = elements[k];
                    auto c = std::make_shared<PbdFemTetConstraint>(m_matType);
                    c->initConstraint(
                        vertices[tet[0]], vertices[tet[1]], vertices[tet[2]], vertices[tet[3]],
                        { m_bodyIndex, tet[0] }, { m_bodyIndex, tet[1] },
                        { m_bodyIndex, tet[2] }, { m_bodyIndex, tet[3] },
                        m_femConfig);
                    constraints.addConstraint(c);
            }, elements.size() > 100);
        }

        void setMaterialType(const PbdFemTetConstraint::MaterialType materialType) { m_matType = materialType; }
        PbdFemTetConstraint::MaterialType getMaterialType() const { return m_matType; }

        void setFemConfig(std::shared_ptr<PbdFemConstraintConfig> femConfig) { m_femConfig = femConfig; }
        std::shared_ptr<PbdFemConstraintConfig> getFemConfig() const { return m_femConfig; }

    protected:
        PbdFemTetConstraint::MaterialType m_matType = PbdFemTetConstraint::MaterialType::StVK;
        std::shared_ptr<PbdFemConstraintConfig> m_femConfig = nullptr;
};

///
/// \struct PbdVolumeConstraintFunctor
///
/// \brief PbdVolumeConstraintFunctor generates constraints per cell of a TetrahedralMesh
///
struct PbdVolumeConstraintFunctor : public PbdBodyConstraintFunctor
{
    public:
        PbdVolumeConstraintFunctor() = default;
        ~PbdVolumeConstraintFunctor() override = default;

        void operator()(PbdConstraintContainer& constraints) override
        {
            // Check for correct mesh type
            CHECK(std::dynamic_pointer_cast<TetrahedralMesh>(m_geom) != nullptr)
                << "PbdVolumeConstraint can only be generated with a TetrahedralMesh";

            // Create constraints
            auto                                     tetMesh     = std::dynamic_pointer_cast<TetrahedralMesh>(m_geom);
            std::shared_ptr<VecDataArray<double, 3>> verticesPtr = m_geom->getVertexPositions();
            const VecDataArray<double, 3>&           vertices    = *verticesPtr;
            std::shared_ptr<VecDataArray<int, 4>>    elementsPtr = tetMesh->getIndices();
            const VecDataArray<int, 4>&              elements    = *elementsPtr;

            ParallelUtils::parallelFor(elements.size(),
                [&](const size_t k)
                {
                    auto& tet = elements[k];
                    auto c    = std::make_shared<PbdVolumeConstraint>();
                    c->initConstraint(
                        vertices[tet[0]], vertices[tet[1]], vertices[tet[2]], vertices[tet[3]],
                        { m_bodyIndex, tet[0] }, { m_bodyIndex, tet[1] },
                        { m_bodyIndex, tet[2] }, { m_bodyIndex, tet[3] },
                        m_stiffness);
                    constraints.addConstraint(c);
            });
        }

        ///
        /// \brief Get/Set the stiffness, how hard the constraint is
        ///@{
        void setStiffness(const double stiffness) { m_stiffness = stiffness; }
        double getStiffness() const { return m_stiffness; }
    ///@}

    protected:
        double m_stiffness = 0.0;
};

///
/// \struct PbdAreaConstraintFunctor
///
/// \brief PbdAreaConstraintFunctor generates constraints per cell of a SurfaceMesh
///
struct PbdAreaConstraintFunctor : public PbdBodyConstraintFunctor
{
    public:
        PbdAreaConstraintFunctor() = default;
        ~PbdAreaConstraintFunctor() override = default;

        void operator()(PbdConstraintContainer& constraints) override
        {
            // Check for correct mesh type
            CHECK(std::dynamic_pointer_cast<SurfaceMesh>(m_geom) != nullptr)
                << "PbdAreaConstraint can only be generated with a SurfaceMesh";

            // ok, now create constraints
            auto                                     triMesh     = std::dynamic_pointer_cast<SurfaceMesh>(m_geom);
            std::shared_ptr<VecDataArray<double, 3>> verticesPtr = m_geom->getVertexPositions();
            const VecDataArray<double, 3>&           vertices    = *verticesPtr;
            std::shared_ptr<VecDataArray<int, 3>>    elemenstPtr = triMesh->getIndices();
            const VecDataArray<int, 3>&              elements    = *elemenstPtr;

            ParallelUtils::parallelFor(elements.size(),
                [&](const size_t k)
                {
                    auto& tri = elements[k];
                    auto c    = std::make_shared<PbdAreaConstraint>();
                    c->initConstraint(vertices[tri[0]], vertices[tri[1]], vertices[tri[2]],
                        { m_bodyIndex, tri[0] }, { m_bodyIndex, tri[1] }, { m_bodyIndex, tri[2] },
                        m_stiffness);
                    constraints.addConstraint(c);
            });
        }

        void addConstraints(PbdConstraintContainer&                     constraints,
                            std::shared_ptr<std::unordered_set<size_t>> vertices) override
        {
            // Check for correct mesh type
            CHECK(std::dynamic_pointer_cast<SurfaceMesh>(m_geom) != nullptr)
                << "Add element constraints does not support current mesh type.";

            auto                                     triMesh = std::dynamic_pointer_cast<SurfaceMesh>(m_geom);
            std::shared_ptr<VecDataArray<double, 3>> initVerticesPtr = m_geom->getInitialVertexPositions();
            std::shared_ptr<VecDataArray<int, 3>>    elementsPtr     = triMesh->getIndices();

            const VecDataArray<double, 3>& initVertices = *initVerticesPtr;
            const VecDataArray<int, 3>&    elements     = *elementsPtr;

            // Build vertex to tri map
            std::vector<std::vector<size_t>> vertexToTriMap(triMesh->getNumVertices());
            for (int k = 0; k < elements.size(); k++)
            {
                const Vec3i& tri = elements[k];
                vertexToTriMap[tri[0]].push_back(k);
                vertexToTriMap[tri[1]].push_back(k);
                vertexToTriMap[tri[2]].push_back(k);
            }
            for (std::vector<size_t>& faceIds : vertexToTriMap)
            {
                std::sort(faceIds.begin(), faceIds.end());
            }

            std::unordered_set<size_t> areaSet;
            for (const size_t& vertIdx : *vertices)
            {
                for (const size_t& triIdx : vertexToTriMap[vertIdx])
                {
                    areaSet.insert(triIdx);
                }
            }

            auto addAreaConstraint =
                [&](size_t k, double stiffness)
                {
                    const Vec3i& tri = elements[k];
                    auto         c   = std::make_shared<PbdAreaConstraint>();
                    c->initConstraint(initVertices[tri[0]], initVertices[tri[1]], initVertices[tri[2]],
                        { m_bodyIndex, tri[0] }, { m_bodyIndex, tri[1] }, { m_bodyIndex, tri[2] },
                        stiffness);
                    constraints.addConstraint(c);
                };

            constraints.reserve(constraints.getConstraints().size() + areaSet.size());
            for (auto& c : areaSet)
            {
                addAreaConstraint(c, m_stiffness);
            }
        }

        ///
        /// \brief Get/Set the stiffness, how hard the constraint is
        ///@{
        void setStiffness(const double stiffness) { m_stiffness = stiffness; }
        double getStiffness() const { return m_stiffness; }
    ///@}

    protected:
        double m_stiffness = 0.0;
};

///
/// \struct PbdBendConstraintFunctor
///
/// \brief PbdBendConstraintFunctor generates constraints per cell of a LineMesh
///
struct PbdBendConstraintFunctor : public PbdBodyConstraintFunctor
{
    public:
        PbdBendConstraintFunctor() = default;
        ~PbdBendConstraintFunctor() override = default;

        void operator()(PbdConstraintContainer& constraints) override
        {
            // Check for correct mesh type
            CHECK(std::dynamic_pointer_cast<LineMesh>(m_geom) != nullptr)
                << "PbdBendConstraint can only be generated with a LineMesh";

            auto                                     lineMesh    = std::dynamic_pointer_cast<LineMesh>(m_geom);
            std::shared_ptr<VecDataArray<double, 3>> verticesPtr = m_geom->getVertexPositions();
            const VecDataArray<double, 3>&           vertices    = *verticesPtr;
            /*std::shared_ptr< VecDataArray<int, 2>> indicesPtr = lineMesh->getIndices();
            const VecDataArray<int, 2>& indices = *indicesPtr*/

            auto addBendConstraint =
                [&](const double k, size_t i1, size_t i2, size_t i3)
                {
                    // i1 should always come first
                    if (i2 < i1)
                    {
                        std::swap(i1, i2);
                    }
                    // i3 should always come last
                    if (i2 > i3)
                    {
                        std::swap(i2, i3);
                    }

                    auto c = std::make_shared<PbdBendConstraint>();
                    if (m_restLengthOverride >= 0.0)
                    {
                        c->initConstraint(vertices[i1], vertices[i2], vertices[i3],
                        { m_bodyIndex, i1 }, { m_bodyIndex, i2 }, { m_bodyIndex, i3 }, m_restLengthOverride, k);
                    }
                    else
                    {
                        c->initConstraint(vertices[i1], vertices[i2], vertices[i3],
                        { m_bodyIndex, i1 }, { m_bodyIndex, i2 }, { m_bodyIndex, i3 }, k);
                    }
                    constraints.addConstraint(c);
                };

            // Iterate sets of stride # of segments
            for (int k = 0; k < vertices.size() - m_stride * 2; k++)
            {
                addBendConstraint(m_stiffness, k, k + m_stride, k + 2 * m_stride);
            }
        }

        ///
        /// \brief Get/Set the stiffness, how hard the constraint is
        ///@{
        void setStiffness(const double stiffness) { m_stiffness = stiffness; }
        double getStiffness() const { return m_stiffness; }
        ///@}

        ///
        /// \brief Get/Set the stride, that is the amount of lines between every bend
        /// constraint
        ///@{
        void setStride(const int stride)
        {
            m_stride = stride;
            CHECK(m_stride >= 1) << "Stride should be at least 1.";
        }

        int getStride() const { return m_stride; }
        ///@}

        ///
        /// \brief Get/Set the rest length, if not specified or negative
        /// the rest length will default to what it was when initialized
        /// @{
        double getRestLength() const { return m_restLengthOverride; }
        void setRestLength(const double restLength) { m_restLengthOverride = restLength; }
    /// @}

    protected:
        double m_stiffness = 0.0;
        int m_stride       = 3;
        double m_restLengthOverride = -1.0;
};

///
/// \struct PbdDihedralConstraintFunctor
///
/// \brief PbdDihedralConstraintFunctor generates constraints per pair of triangle neighbors in a
/// SurfaceMesh
///
struct PbdDihedralConstraintFunctor : public PbdBodyConstraintFunctor
{
    public:
        PbdDihedralConstraintFunctor() = default;
        ~PbdDihedralConstraintFunctor() override = default;

        void operator()(PbdConstraintContainer& constraints) override
        {
            // Check for correct mesh type
            CHECK(std::dynamic_pointer_cast<SurfaceMesh>(m_geom) != nullptr)
                << "PbdDihedralConstraint can only be generated with a SurfaceMesh";

            // Create constraints
            auto                                     triMesh     = std::dynamic_pointer_cast<SurfaceMesh>(m_geom);
            std::shared_ptr<VecDataArray<double, 3>> verticesPtr = triMesh->getVertexPositions();
            const VecDataArray<double, 3>&           vertices    = *verticesPtr;
            std::shared_ptr<VecDataArray<int, 3>>    elementsPtr = triMesh->getIndices();
            const VecDataArray<int, 3>&              elements    = *elementsPtr;
            const int                                nV = triMesh->getNumVertices();
            std::vector<std::vector<int>>            vertIdsToTriangleIds(nV);

            for (int k = 0; k < elements.size(); ++k)
            {
                const Vec3i& tri = elements[k];
                vertIdsToTriangleIds[tri[0]].push_back(k);
                vertIdsToTriangleIds[tri[1]].push_back(k);
                vertIdsToTriangleIds[tri[2]].push_back(k);
            }

            // Used to resolve duplicates
            std::vector<std::vector<bool>> E(nV, std::vector<bool>(nV, 1));

            auto addDihedralConstraint =
                [&](const std::vector<int>& r1, const std::vector<int>& r2,
                    const int k, int i1, int i2)
                {
                    if (i1 > i2) // Make sure i1 is always smaller than i2
                    {
                        std::swap(i1, i2);
                    }
                    if (E[i1][i2])
                    {
                        E[i1][i2] = 0;

                        // Find the shared edge
                        std::vector<size_t> rs(2);
                        auto                it = std::set_intersection(r1.begin(), r1.end(), r2.begin(), r2.end(), rs.begin());
                        rs.resize(static_cast<size_t>(it - rs.begin()));
                        if (rs.size() > 1)
                        {
                            size_t      idx  = (rs[0] == k) ? 1 : 0;
                            const auto& tri0 = elements[k];
                            const auto& tri1 = elements[rs[idx]];
                            size_t      idx0 = 0;
                            size_t      idx1 = 0;
                            for (size_t i = 0; i < 3; ++i)
                            {
                                if (tri0[i] != i1 && tri0[i] != i2)
                                {
                                    idx0 = tri0[i];
                                }
                                if (tri1[i] != i1 && tri1[i] != i2)
                                {
                                    idx1 = tri1[i];
                                }
                            }
                            auto c = std::make_shared<PbdDihedralConstraint>();
                            c->initConstraint(vertices[idx0], vertices[idx1], vertices[i1], vertices[i2],
                                { m_bodyIndex, idx0 }, { m_bodyIndex, idx1 },
                                { m_bodyIndex, i1 }, { m_bodyIndex, i2 }, m_stiffness);
                            constraints.addConstraint(c);
                        }
                    }
                };

            for (size_t i = 0; i < vertIdsToTriangleIds.size(); i++)
            {
                std::sort(vertIdsToTriangleIds[i].begin(), vertIdsToTriangleIds[i].end());
            }

            // For every triangle
            for (int k = 0; k < elements.size(); ++k)
            {
                const Vec3i& tri = elements[k];

                // Get all the neighbor triangles (to the vertices)
                std::vector<int>& neighborTriangles0 = vertIdsToTriangleIds[tri[0]];
                std::vector<int>& neighborTriangles1 = vertIdsToTriangleIds[tri[1]];
                std::vector<int>& neighborTriangles2 = vertIdsToTriangleIds[tri[2]];

                // Add constraints between all the triangles
                addDihedralConstraint(neighborTriangles0, neighborTriangles1, k, tri[0], tri[1]);
                addDihedralConstraint(neighborTriangles0, neighborTriangles2, k, tri[0], tri[2]);
                addDihedralConstraint(neighborTriangles1, neighborTriangles2, k, tri[1], tri[2]);
            }
        }

        void addConstraints(PbdConstraintContainer&                     constraints,
                            std::shared_ptr<std::unordered_set<size_t>> vertices) override
        {
            // Check for correct mesh type
            CHECK(std::dynamic_pointer_cast<SurfaceMesh>(m_geom) != nullptr)
                << "Add element constraints does not support current mesh type.";

            auto                                     triMesh = std::dynamic_pointer_cast<SurfaceMesh>(m_geom);
            std::shared_ptr<VecDataArray<double, 3>> initVerticesPtr = m_geom->getInitialVertexPositions();
            std::shared_ptr<VecDataArray<int, 3>>    elementsPtr     = triMesh->getIndices();

            const VecDataArray<double, 3>& initVertices = *initVerticesPtr;
            const VecDataArray<int, 3>&    elements     = *elementsPtr;

            // Build vertex to tri map
            std::vector<std::vector<size_t>> vertexToTriMap(triMesh->getNumVertices());
            for (int k = 0; k < elements.size(); k++)
            {
                const Vec3i& tri = elements[k];
                vertexToTriMap[tri[0]].push_back(k);
                vertexToTriMap[tri[1]].push_back(k);
                vertexToTriMap[tri[2]].push_back(k);
            }
            for (std::vector<size_t>& faceIds : vertexToTriMap)
            {
                std::sort(faceIds.begin(), faceIds.end());
            }

            std::map<std::pair<size_t, size_t>, std::pair<size_t, size_t>> dihedralSet;
            for (const size_t& vertIdx : *vertices)
            {
                for (const size_t& triIdx : vertexToTriMap[vertIdx])
                {
                    const Vec3i& tri = elements[triIdx];
                    for (size_t i = 0; i < 3; i++)
                    {
                        size_t j  = (i + 1) % 3;
                        size_t i0 = tri[i];
                        size_t i1 = tri[j];
                        if (i0 > i1)
                        {
                            std::swap(i0, i1);
                        }
                        const std::vector<size_t>& r0 = vertexToTriMap[i0];
                        const std::vector<size_t>& r1 = vertexToTriMap[i1];
                        std::vector<size_t>        rs(2);
                        auto                       it = std::set_intersection(r0.begin(), r0.end(), r1.begin(), r1.end(), rs.begin());
                        rs.resize(static_cast<size_t>(it - rs.begin()));
                        if (rs.size() > 1)
                        {
                            dihedralSet[std::make_pair(i0, i1)] = std::make_pair(rs[0], rs[1]);
                        }
                    }
                }
            }

            auto addDihedralConstraint =
                [&](int t0, int t1, int i1, int i2, double stiffness)
                {
                    const Vec3i& tri0 = elements[t0];
                    const Vec3i& tri1 = elements[t1];
                    size_t       idx0 = 0;
                    size_t       idx1 = 0;
                    for (int i = 0; i < 3; ++i)
                    {
                        if (tri0[i] != i1 && tri0[i] != i2)
                        {
                            idx0 = static_cast<size_t>(tri0[i]);
                        }
                        if (tri1[i] != i1 && tri1[i] != i2)
                        {
                            idx1 = static_cast<size_t>(tri1[i]);
                        }
                    }
                    auto c = std::make_shared<PbdDihedralConstraint>();
                    c->initConstraint(initVertices[idx0], initVertices[idx1], initVertices[i1], initVertices[i2],
                        { m_bodyIndex, idx0 }, { m_bodyIndex, idx1 },
                        { m_bodyIndex, static_cast<size_t>(i1) }, { m_bodyIndex, static_cast<size_t>(i2) },
                        stiffness);
                    constraints.addConstraint(c);
                };

            constraints.reserve(constraints.getConstraints().size() + dihedralSet.size());
            for (auto& c : dihedralSet)
            {
                addDihedralConstraint(
                    static_cast<int>(c.second.first), static_cast<int>(c.second.second),
                    static_cast<int>(c.first.first), static_cast<int>(c.first.second),
                    m_stiffness);
            }
        }

        ///
        /// \brief Get/Set the stiffness, how hard the constraint is
        ///@{
        void setStiffness(const double stiffness) { m_stiffness = stiffness; }
        double getStiffness() const { return m_stiffness; }
    ///@}

    protected:
        double m_stiffness = 0.0;
};

///
/// \struct PbdConstantDensityConstraintFunctor
///
/// \brief PbdConstantDensityConstraintFunctor generates a global constant density constraint
/// for an entire PointSet
///
struct PbdConstantDensityConstraintFunctor : public PbdBodyConstraintFunctor
{
    public:
        PbdConstantDensityConstraintFunctor() = default;
        ~PbdConstantDensityConstraintFunctor() override = default;

        void operator()(PbdConstraintContainer& constraints) override
        {
            // Check for correct mesh type
            CHECK(std::dynamic_pointer_cast<PointSet>(m_geom) != nullptr)
                << "PbdConstantDensityConstraint can only be generated with a PointSet";

            auto c = std::make_shared<PbdConstantDensityConstraint>();
            c->initConstraint(*m_geom->getVertexPositions(), m_stiffness);
            constraints.addConstraint(c);
        }

        ///
        /// \brief Get/Set the stiffness, how hard the constraint is
        ///@{
        void setStiffness(const double stiffness) { m_stiffness = stiffness; }
        double getStiffness() const { return m_stiffness; }
    ///@}

    protected:
        double m_stiffness = 0.0;
};
} // namespace imstk