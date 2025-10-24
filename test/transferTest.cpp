#include <gtest/gtest.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h>
#include <geometrycentral/surface/signpost_intrinsic_triangulation.h>
#include <filesystem>

#include "polyscope/curve_network.h"
#include <chrono>

#include "eider/poisson.h"
#include <implot.h>

#include "geometrycentral/surface/transfer_functions.h"

using namespace geometrycentral::surface;
using namespace geometrycentral;
#include <eider/refine.h>

VertexData<Vector3> intrinsic_position(IntrinsicTriangulation& Tri, VertexPositionGeometry& inputG){
    auto& mesh = *Tri.intrinsicMesh;
    VertexData<Vector3> int_positions(mesh) ;
    for (Vertex v : mesh.vertices()) { int_positions[v] = Tri.vertexLocations[v].interpolate(inputG.vertexPositions); }
    return int_positions;
}

TEST(transfertTest,testL2Face) {
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" / "disk.stl";
    auto [parent_m,parent_g] = readManifoldSurfaceMesh(fds.string());

    IntegerCoordinatesIntrinsicTriangulation icit(*parent_m,*parent_g);
    icit.delaunayRefine(25,0.01);


    ManifoldSurfaceMesh& m = *icit.intrinsicMesh; m.compress();
    VertexData<Vector3> int_positions(m) ;
    for (Vertex v : m.vertices()) {
        int_positions[v] = icit.vertexLocations[v].interpolate(parent_g->vertexPositions);
    }

    FaceData<double> f_A(*parent_m), f_B(m);
    for (const Face f: parent_m->faces()) {
        f_A[f] = SurfacePoint(f,Vector3(1./3,1./3,1./3)).interpolate(parent_g->vertexPositions).x;
    }
    f_B = transferAtoB_L2(icit,f_A);
    polyscope::init();
    polyscope::SurfaceMesh* pm_A = polyscope::registerSurfaceMesh("original", parent_g->vertexPositions,parent_m->getFaceVertexList());
    polyscope::SurfaceMesh* pm_B = polyscope::registerSurfaceMesh("Intrinsic", int_positions,m.getFaceVertexList());
    pm_A->addFaceScalarQuantity("f",f_A);
    pm_B->addFaceScalarQuantity("f",f_B);
    pm_A->translate(glm::vec3(-1,0,0));
    pm_B->translate(glm::vec3( 1,0,0));


    polyscope::show();
}

TEST(transfertTest,deadVerticesRetainValues) {
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" / "quad.stl";
    auto [parent_m,parent_g] = readManifoldSurfaceMesh(fds.string());

    AdaptiveTriangulation atri(*parent_m,*parent_g);
    ManifoldSurfaceMesh& m = atri.mesh();
    VertexData<std::size_t> idx(m,0);
    std::vector<Face> faces { }; for (Face f: m.faces()) faces.push_back(f);
    atri.refine(faces);
    Vertex center_v;
    for (Vertex v: m.vertices()) if(!v.isBoundary()) center_v = v;
    ASSERT_NE(center_v,Vertex());
    idx[center_v] = 1;
    faces.clear();for (Face f: m.faces()) faces.push_back(f);
    Face f_back = m.face(m.nFaces()-2);
    atri.coarse(faces);
    ASSERT_TRUE(center_v.isDead());
    ASSERT_EQ(idx[center_v],1);
    ASSERT_TRUE(f_back.isDead());

    int i = 0;
    for (Halfedge he: f_back.adjacentHalfedges()) {
        ASSERT_TRUE(he.isDead()); i++;
    }
    ASSERT_EQ(i,3);
}



// We rely on the following undocumented behaviour of geometry central:
// The tangent space on each face f is build, s.t. halfedge(f) lies on the x  axis,
// with Tangent(he_f)[0] / length[he_f] = 1
TEST(transferTest,testAdaptiveTangentSpace) {
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" / "torus.stl";
    auto [parent_m,parent_g] = readManifoldSurfaceMesh(fds.string());
    AdaptiveTriangulation atri(*parent_m, *parent_g);

    ManifoldSurfaceMesh& m = atri.mesh();
    IntrinsicGeometryInterface& g = atri.geom();

    g.requireHalfedgeVectorsInFace();
    std::vector<Face> faces;
    for (Face f: m.faces()) {
        Vector2 vec_he = g.halfedgeVectorsInFace[f.halfedge()];
        ASSERT_EQ(vec_he[1],0);
        ASSERT_NEAR(vec_he[0]/g.edgeLengths[f.halfedge().edge()],1,0.0001);
        faces.push_back(f);
    }

    g.unrequireHalfedgeVectorsInFace();
    atri.refine(faces);
    g.refreshQuantities(); g.requireHalfedgeVectorsInFace();
    faces.clear();

    for (Face f: m.faces()) {
        Vector2 vec_he = g.halfedgeVectorsInFace[f.halfedge()];
        ASSERT_EQ(vec_he[1],0);
        ASSERT_NEAR(vec_he[0]/g.edgeLengths[f.halfedge().edge()],1,0.0001);
        faces.push_back(f);
    }

    g.unrequireHalfedgeVectorsInFace();
    atri.coarse(faces);
    g.refreshQuantities(); g.requireHalfedgeVectorsInFace();
    faces.clear();

    for (Face f: m.faces()) {
        Vector2 vec_he = g.halfedgeVectorsInFace[f.halfedge()];
        ASSERT_EQ(vec_he[1],0);
        ASSERT_NEAR(vec_he[0]/g.edgeLengths[f.halfedge().edge()],1,0.0001);
        faces.push_back(f);
    }
}

TEST(transferTest,testAdaptiveTangentSpaceSplitHe) {
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" / "triangle.stl";
    auto [parent_m,parent_g] = readManifoldSurfaceMesh(fds.string());
    AdaptiveTriangulation atri(*parent_m, *parent_g);

    ManifoldSurfaceMesh& m = atri.mesh();
    IntrinsicGeometryInterface& g = atri.geom();

    Face f = m.face(0); Halfedge he = f.halfedge();

    // mark opposite corner of f.he
    for (Halfedge h: f.adjacentHalfedges()) atri.marked_corner[h.corner()] = false;
    atri.marked_corner[he.corner()] = true;

    g.requireHalfedgeVectorsInFace();
    auto& g_tang = g.halfedgeVectorsInFace;
    Vector2 he_tang = g_tang[he];
    Vector2 (g_tang[he] + g_tang[he.next()]) - g_tang[he] * 0.5;
    atri.refine({f});
    g.refreshQuantities();
    ASSERT_EQ(m.nFaces(),2);
    Face f0 = m.face(0), f1 = m.face(1);
    // Split sets faces to orthogonal halfedges
    ASSERT_EQ(f0.halfedge().twin(), f1.halfedge());
}

TEST(transferTest,testVectorRotation) {
    Vector2 x = Vector2(1,0);
    Vector2 r = Vector2::fromAngle(PI/4);
    ASSERT_EQ(r*x, r);
    Vector2 y = r*x;
    ASSERT_EQ(y/x,r);
    std::cout << y/x <<std::endl;
}

TEST(transferTest,testSplitPreserveFace) {
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" / "quad.stl";
    auto [parent_m,parent_g] = readManifoldSurfaceMesh(fds.string());
    ManifoldSurfaceMesh& m = *parent_m;
    IntrinsicGeometryInterface& g = *parent_g;
    FaceData<bool> fd(m,false);

    Halfedge he;
    for (Edge e: m.edges() ) if (!e.isBoundary()) he = e.halfedge(); // should also hold for twin
    for (Face f: he.edge().adjacentFaces()) fd[f] = true;
    Vertex vi = he.tailVertex(), vj = he.tipVertex();
    Halfedge new_he = m.splitEdgeTriangular(he.edge());
    ASSERT_EQ(new_he.tipVertex(), vj);
    ASSERT_TRUE(!new_he.tailVertex().isBoundary());
    ASSERT_FALSE(fd[new_he.face()]);
    ASSERT_TRUE(fd[new_he.prevOrbitFace().twin().face()]);
    ASSERT_EQ(new_he.face().halfedge(),new_he.prevOrbitFace());
    ASSERT_EQ(new_he.prevOrbitFace().twin().face().halfedge(),new_he.prevOrbitFace().twin());
}

TEST(transfertTest,testL2FaceU) {
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" / "quad.stl";
    auto [parent_m,parent_g] = readManifoldSurfaceMesh(fds.string());
    ManifoldSurfaceMesh& m = *parent_m; VertexPositionGeometry& g = *parent_g;
    g.requireHalfedgeVectorsInFace();
    g.requireFaceTangentBasis();

    Edge e;
    for (Edge me: m.edges()) if(!me.isBoundary()) e= me;

    FaceData<Vector2> u(m);
    CornerData<bool> corner(m, false);
    for (Halfedge he: e.adjacentHalfedges()) {
        u[he.face()] = g.halfedgeVectorsInFace[he];
        corner[he.prevOrbitFace().corner()] = true;
    }

    Vertex vi = e.halfedge().tailVertex(), vj = e.halfedge().tipVertex();
    AdaptiveFaceTransfer transfer(m, g, u, corner);

    transfer.startRefine();
    Quad q (e.halfedge(),g.halfedgeVectorsInFace);
    Halfedge he = m.splitEdgeTriangular(e);
    g.refreshQuantities();
    Diamond d ( he,g.halfedgeVectorsInFace);
    transfer.refineEdge(q,d);
    transfer.endRefine();

    transfer.startCoarse();
     m.collapseEdgeTriangular(he);
     for (Halfedge ohe: vi.outgoingHalfedges()){ if(ohe.tipVertex() == vj) he = ohe; }
     g.refreshQuantities();
     transfer.coarseEdge(Quad(he,g.halfedgeVectorsInFace),d);
    transfer.endCoarse();


    HalfedgeData<double> frame_base(m,0);
    for (Face f: m.faces()) {
        frame_base[f.halfedge()] =1;
    }

    auto u_new = transfer.transfer();

    m.compress();
    g.refreshQuantities();

    polyscope::init();
    polyscope::SurfaceMesh* pm_A = polyscope::registerSurfaceMesh("original", g.vertexPositions,m.getFaceVertexList(),polyscopePermutations(m));
    FaceData<Vector3> e1(m),e2(m);
    for (Face f: m.faces()) { e1[f] = g.faceTangentBasis[f][0], e2[f] = g.faceTangentBasis[f][1];  }
    pm_A->addFaceTangentVectorQuantity("unrotated",u,e1,e2)->setEnabled(true);
    pm_A->addFaceTangentVectorQuantity("rotated",u_new,e1,e2)->setEnabled(true);
    pm_A->addHalfedgeScalarQuantity("frame",frame_base);

    polyscope::show();
}

TEST(transfertTest,testL2VerticesCoarse) {
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" / "grid.stl";
    auto [parent_m,parent_g] = readManifoldSurfaceMesh(fds.string());

    AdaptiveTriangulation atri(*parent_m, *parent_g);

    ManifoldSurfaceMesh& m = atri.mesh();
    IntrinsicGeometryInterface& g = atri.geom();



    std::vector<Face> faces;
    for (Face f: m.faces()){ faces.push_back(f); }
    atri.refine(faces);

    VertexData<Vector3> int_positions_A = intrinsic_position(atri.intrinsicTriangulation(),*parent_g);
    VertexData<double> fA(m,0);
    for(Vertex v: m.vertices()){
        fA[v] = int_positions_A[v].x;
    }

    std::vector<Face> facesR, facesC; for (Face f: m.faces()){
        Vector3 pos = SurfacePoint(f,Vector3::constant(1./3)).interpolate(int_positions_A);
        if(pos.x < 0) facesR.push_back(f);
        else facesC.push_back(f);
    }

    AdaptiveVertexTransfer transfer(atri.intrinsicTriangulation(),fA);
    atri.refine(facesR,&transfer);
    atri.coarse(facesC,&transfer);
    VertexData<double> fB = transfer.transfer();

    int_positions_A = intrinsic_position(atri.intrinsicTriangulation(),*parent_g);
    for(Vertex v: m.vertices()){
        ASSERT_NEAR(fB[v],int_positions_A[v].x,0.000001);
    }
}



