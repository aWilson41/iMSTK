#include "imstkPbdInteractionPair.h"
#include "imstkSurfaceMesh.h"

namespace imstk {

bool isIntersect(const double& a, const double& b, const double& c, const double& d)
{
    if ( (a <= d && a >= c) || (c <= b && c >= a) ) return true;
    return false;
}

bool testAABBvsAABB(const double& min1_x, const double& max1_x,
                    const double& min1_y, const double& max1_y,
                    const double& min1_z, const double& max1_z,
                    const double& min2_x, const double& max2_x,
                    const double& min2_y, const double& max2_y,
                    const double& min2_z, const double& max2_z)
{
    return  isIntersect(min1_x, max1_x, min2_x, max2_x) &&
            isIntersect(min1_y, max1_y, min2_y, max2_y) &&
            isIntersect(min1_z, max1_z, min2_z, max2_z);
}

bool testLINEvsLINE(const double& x1, const double& y1, const double& z1,
                    const double& x2, const double& y2, const double& z2,
                    const double& x3, const double& y3, const double& z3,
                    const double& x4, const double& y4, const double& z4,
                    const double& prox1, const double& prox2)
{
    double min1_x, max1_x, min1_y, max1_y, min1_z, max1_z;
    if ( x1 < x2 ) { min1_x = x1; max1_x = x2;}
    else { min1_x = x2; max1_x = x1;}
    if ( y1 < y2 ) { min1_y = y1; max1_y = y2;}
    else { min1_y = y2; max1_y = y1;}
    if ( z1 < z2 ) { min1_z = z1; max1_z = z2;}
    else { min1_z = z2; max1_z = z1;}
    double min2_x, max2_x, min2_y, max2_y, min2_z, max2_z;
    if ( x3 < x4 ) { min2_x = x3; max2_x = x4;}
    else { min2_x = x4; max2_x = x3;}
    if ( y3 < y4 ) { min2_y = y3; max2_y = y4;}
    else { min2_y = y4; max2_y = y3;}
    if ( z3 < z4 ) { min2_z = z3; max2_z = z4;}
    else { min2_z = z4; max2_z = z3;}

    return  testAABBvsAABB(min1_x - prox1, max1_x + prox1, min1_y - prox1, max1_y + prox1, min1_z - prox1, max1_z + prox1,
                           min2_x - prox2, max2_x + prox2, min2_y - prox2, max2_y + prox2, min2_z - prox2, max2_z + prox2);

}

bool testPOINTvsTRIANGLE(const double& x1, const double& y1, const double& z1,
                    const double& x2, const double& y2, const double& z2,
                    const double& x3, const double& y3, const double& z3,
                    const double& x4, const double& y4, const double& z4,
                    const double& prox1, const double& prox2)
{
    double min_x, max_x, min_y, max_y, min_z, max_z;
    min_x = std::min(x2, std::min(x3,x4));
    max_x = std::max(x2, std::max(x3,x4));
    min_y = std::min(y2, std::min(y3,y4));
    max_y = std::max(y2, std::max(y3,y4));
    min_z = std::min(z2, std::min(z3,z4));
    max_z = std::max(z2, std::max(z3,z4));
    return  testAABBvsAABB(x1 - prox1, x1 + prox1, y1 - prox1, y1 + prox1, z1 - prox1, z1 + prox1,
                           min_x - prox2, max_x + prox2, min_y - prox2, max_y + prox2, min_z - prox2, max_z + prox2);
}

bool PbdInteractionPair::doBroadPhase()
{
    auto g1 = first->getCollidingGeometry();
    auto g2 = second->getCollidingGeometry();
    auto mesh1 = std::static_pointer_cast<Mesh>(g1);
    auto mesh2 = std::static_pointer_cast<Mesh>(g2);
    Vec3d min1, max1;
    mesh1->computeBoundingBox(min1, max1);
    Vec3d min2, max2;
    mesh2->computeBoundingBox(min2, max2);
    double prox1 = first->getProximity();
    double prox2 = second->getProximity();
    return  testAABBvsAABB(min1_x - prox1, max1_x + prox1, min1_y - prox1, max1_y + prox1, min1_z - prox1, max1_z + prox1,
                           min2_x - prox2, max2_x + prox2, min2_y - prox2, max2_y + prox2, min2_z - prox2, max2_z + prox2);

}

void PbdInteractionPair::doNarrowPhase()
{
    auto g1 = first->getCollidingGeometry();
    auto g2 = second->getCollidingGeometry();
    auto mesh1 = std::static_pointer_cast<SurfaceMesh>(g1);
    auto mesh2 = std::static_pointer_cast<SurfaceMesh>(g2);

    double prox1 = first->getProximity();
    double prox2 = second->getProximity();

    // brute force, use BVH or spatial grid woulbe be much better
    // point
    for (int i = 0; i < mesh1->getNumVertices(); ++i) {
        Vec3d& p = mesh1->getVertexPosition(i);
        std::vector<SurfaceMesh::TriangleArray> elements = mesh2->getTrianglesVertices();
        for (int j = 0; j < elements.size(); ++j) {
            Vec3i& e = elements[j];
            Vec3d& p0 = mesh2->getVertexPosition(e[0]);
            Vec3d& p1 = mesh2->getVertexPosition(e[1]);
            Vec3d& p2 = mesh2->getVertexPosition(e[2]);

            if (testPOINTvsTRIANGLE(p[0],p[1],p[2],
                                    p0[0],p0[1],p0[2],
                                    p1[0],p1[1],p1[2],
                                    p2[0],p2[1],p2[2], prox1, prox2))
            {
                PointTriangleConstraint* c = new PointTriangleConstraint;
                c->initConstraint(first, i, second, e[0], e[1], e[2]);
                m_collisionConstraints.push_back(c);
            }
        }
    }
    // edge
    // since we don't have edge structure, the following is not good
    int nV = mesh1->getNumVertices();
    std::vector<std::vector<bool>> E(nV, std::vector<bool>(nV, 1));
    std::vector<SurfaceMesh::TriangleArray> elements = mesh1->getTrianglesVertices();
    for (int k = 0; k < elements.size(); ++k) {
       Vec3i& tri = elements[k];
       unsigned int i1;
       unsigned int i2;
       i1 = tri[0];
       i2 = tri[1];
       if (E[i1][i2] && E[i2][i1]) {
           Vec3d& P = mesh1->getVertexPosition(i1);
           Vec3d& Q = mesh1->getVertexPosition(i2);
           std::vector<SurfaceMesh::TriangleArray> elements2 = mesh2->getTrianglesVertices();
           for (int j = 0; j < elements2.size(); ++j) {
               Vec3i& e = elements2[j];
               Vec3d& p0 = mesh2->getVertexPosition(e[0]);
               Vec3d& p1 = mesh2->getVertexPosition(e[1]);
               Vec3d& p2 = mesh2->getVertexPosition(e[2]);
               if (testLINEvsLINE(P[0],P[1],P[2],
                                  Q[0],Q[1],Q[2],
                                  p0[0],p0[1],p0[2],
                                  p1[0],p1[1],p1[2], prox1, prox2))
               {
                   EdgeEdgeConstraint* c = new EdgeEdgeConstraint;
                   c->initConstraint(first, i1, i2, second, e[0], e[1]);
                   m_collisionConstraints.push_back(c);
               }
               if (testLINEvsLINE(P[0],P[1],P[2],
                                  Q[0],Q[1],Q[2],
                                  p1[0],p1[1],p1[2],
                                  p2[0],p2[1],p2[2], prox1, prox2))
               {
                   EdgeEdgeConstraint* c = new EdgeEdgeConstraint;
                   c->initConstraint(first, i1, i2, second, e[1], e[2]);
                   m_collisionConstraints.push_back(c);
               }
               if (testLINEvsLINE(P[0],P[1],P[2],
                                  Q[0],Q[1],Q[2],
                                  p2[0],p2[1],p2[2],
                                  p0[0],p0[1],p0[2], prox1, prox2))
               {
                   EdgeEdgeConstraint* c = new EdgeEdgeConstraint;
                   c->initConstraint(first, i1, i2, second, e[2], e[0]);
                   m_collisionConstraints.push_back(c);
               }
           }
           E[i1][i2] = 0;
       }

       i1 = tri[1];
       i2 = tri[2];
       if (E[i1][i2] && E[i2][i1]) {
           Vec3d& P = mesh1->getVertexPosition(i1);
           Vec3d& Q = mesh1->getVertexPosition(i2);
           std::vector<SurfaceMesh::TriangleArray> elements2 = mesh2->getTrianglesVertices();
           for (int j = 0; j < elements2.size(); ++j) {
               Vec3i& e = elements2[j];
               Vec3d& p0 = mesh2->getVertexPosition(e[0]);
               Vec3d& p1 = mesh2->getVertexPosition(e[1]);
               Vec3d& p2 = mesh2->getVertexPosition(e[2]);
               if (testLINEvsLINE(P[0],P[1],P[2],
                                  Q[0],Q[1],Q[2],
                                  p0[0],p0[1],p0[2],
                                  p1[0],p1[1],p1[2], prox1, prox2))
               {
                   EdgeEdgeConstraint* c = new EdgeEdgeConstraint;
                   c->initConstraint(first, i1, i2, second, e[0], e[1]);
                   m_collisionConstraints.push_back(c);
               }
               if (testLINEvsLINE(P[0],P[1],P[2],
                                  Q[0],Q[1],Q[2],
                                  p1[0],p1[1],p1[2],
                                  p2[0],p2[1],p2[2], prox1, prox2))
               {
                   EdgeEdgeConstraint* c = new EdgeEdgeConstraint;
                   c->initConstraint(first, i1, i2, second, e[1], e[2]);
                   m_collisionConstraints.push_back(c);
               }
               if (testLINEvsLINE(P[0],P[1],P[2],
                                  Q[0],Q[1],Q[2],
                                  p2[0],p2[1],p2[2],
                                  p0[0],p0[1],p0[2], prox1, prox2))
               {
                   EdgeEdgeConstraint* c = new EdgeEdgeConstraint;
                   c->initConstraint(first, i1, i2, second, e[2], e[0]);
                   m_collisionConstraints.push_back(c);
               }
           }
           E[i1][i2] = 0;
       }
       i1 = tri[2];
       i2 = tri[0];
       if (E[i1][i2] && E[i2][i1]) {
           Vec3d& P = mesh1->getVertexPosition(i1);
           Vec3d& Q = mesh1->getVertexPosition(i2);
           std::vector<SurfaceMesh::TriangleArray> elements2 = mesh2->getTrianglesVertices();
           for (int j = 0; j < elements2.size(); ++j) {
               Vec3i& e = elements2[j];
               Vec3d& p0 = mesh2->getVertexPosition(e[0]);
               Vec3d& p1 = mesh2->getVertexPosition(e[1]);
               Vec3d& p2 = mesh2->getVertexPosition(e[2]);
               if (testLINEvsLINE(P[0],P[1],P[2],
                                  Q[0],Q[1],Q[2],
                                  p0[0],p0[1],p0[2],
                                  p1[0],p1[1],p1[2], prox1, prox2))
               {
                   EdgeEdgeConstraint* c = new EdgeEdgeConstraint;
                   c->initConstraint(first, i1, i2, second, e[0], e[1]);
                   m_collisionConstraints.push_back(c);
               }
               if (testLINEvsLINE(P[0],P[1],P[2],
                                  Q[0],Q[1],Q[2],
                                  p1[0],p1[1],p1[2],
                                  p2[0],p2[1],p2[2], prox1, prox2))
               {
                   EdgeEdgeConstraint* c = new EdgeEdgeConstraint;
                   c->initConstraint(first, i1, i2, second, e[1], e[2]);
                   m_collisionConstraints.push_back(c);
               }
               if (testLINEvsLINE(P[0],P[1],P[2],
                                  Q[0],Q[1],Q[2],
                                  p2[0],p2[1],p2[2],
                                  p0[0],p0[1],p0[2], prox1, prox2))
               {
                   EdgeEdgeConstraint* c = new EdgeEdgeConstraint;
                   c->initConstraint(first, i1, i2, second, e[2], e[0]);
                   m_collisionConstraints.push_back(c);
               }
           }
           E[i1][i2] = 0;
       }
    }
}

void PbdInteractionPair::doCollision()
{
    if (!m_collisionConstraints.empty()) {
        int i = 0;
        while (++i < maxIter){
            for (int k = 0; k < m_collisionConstraints.size(); ++k) {
                m_collisionConstraints[k]->solvePositionConstraint();
            }
        }
    }
}



}
