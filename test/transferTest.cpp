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

bool selectFace(polyscope::SurfaceMesh* myMesh, std::size_t* index) {
    myMesh->setSelectionMode(polyscope::MeshSelectionMode::FacesOnly);

    // get the mouse location from ImGui
    ImGuiIO &io = ImGui::GetIO();
    if (io.MouseClicked[0]) {
        // if clicked
        glm::vec2 screenCoords{io.MousePos.x, io.MousePos.y};
        polyscope::PickResult pickResult = polyscope::pickAtScreenCoords(screenCoords);

        // check out pickResult.isHit, pickResult.structureName, pickResult.depth, etc

        // get additional information if we clicked on a mesh
        if (pickResult.isHit && pickResult.structure == myMesh) {
            polyscope::SurfaceMeshPickResult meshPickResult =
                myMesh->interpretPickResult(pickResult);

            if (meshPickResult.elementType == polyscope::MeshElement::FACE) {
                *index = meshPickResult.index;
                return true;
            }
        }
    }
    return false;
}

std::array<FaceData<Vector3>,2> face_tangent_basis(ManifoldSurfaceMesh& m, VertexData<Vector3>& p) {
    VertexPositionGeometry g(m,p);
    g.requireFaceTangentBasis();
    FaceData<Vector3> e1(m),e2(m);
    for (Face f: m.faces()) { e1[f] = g.faceTangentBasis[f][0], e2[f] = g.faceTangentBasis[f][1];  }
    return {e1, e2};
}

std::vector<Face> selectedFaces(FaceData<int>& d) {
    std::vector<Face> faces;
    for (Face f: d.getMesh()->faces()) {
        if (d[f]) faces.push_back(f);
    }
    return faces;
}

TEST(transfertTest,testL2FaceU) {
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" / "quad.stl";
    auto [parent_m,parent_g] = readManifoldSurfaceMesh(fds.string());
    AdaptiveTriangulation atri(*parent_m, *parent_g);
    ManifoldSurfaceMesh& m = atri.mesh();
    IntrinsicGeometryInterface& g = atri.geom();

    Edge ce; for (Edge e: m.edges() ) if (!e.isBoundary()) ce = e;
    atri.marked_corner.fill(false);
    atri.marked_corner[ce.halfedge().prevOrbitFace().corner()] = true;
    atri.marked_corner[ce.halfedge().twin().next().corner()] = true;

    auto v_pos = intrinsic_position(atri.intrinsicTriangulation(),*parent_g);
    auto tang_basis  = face_tangent_basis(m,v_pos);
    FaceData<Vector2> u(m);
    for (Face f: m.faces()) {
        Vector3 e1 = tang_basis[0][f], e2 = tang_basis[1][f];
        Vector2 f1 = g.halfedgeVectorsInFace[f.halfedge()], f2 = f1.rotate90();
        Vector3 f1w = f1.x* e1 + f1.y*e2, f2w = f2.x*e1 + f2.y*e2;
        u[f] = Vector2(
            dot(f1w,Vector3(1,0,0)),
            dot(f2w,Vector3(1,0,0))
        );
    }

    g.requireHalfedgeVectorsInFace();

    polyscope::init();

    HalfedgeData<int> f_basis(m,0), marked_halfedges(m,0);
    auto vis = [&]()->polyscope::SurfaceMesh* {
        g.refreshQuantities();
        m.compress();
        v_pos = intrinsic_position(atri.intrinsicTriangulation(),*parent_g);
        tang_basis  = face_tangent_basis(m,v_pos);
        polyscope::SurfaceMesh* pm_A = polyscope::registerSurfaceMesh("original", v_pos,m.getFaceVertexList(),polyscopePermutations(m));
        pm_A->addFaceTangentVectorQuantity("vector field",u,tang_basis[0],tang_basis[1])->setEnabled(true);

        f_basis.fill(0);
        for (Face f:m.faces()){f_basis[f.halfedge()] = 1;}
        pm_A->addHalfedgeScalarQuantity("Base Halfedges",f_basis);

        marked_halfedges.fill(0);
        for (Halfedge he: m.halfedges()) if (atri.marked_corner[he.oppositeCorner()]) marked_halfedges[he] = 1;
        pm_A->addHalfedgeScalarQuantity("Marked Halfedges",marked_halfedges);

        return pm_A;
    };

    polyscope::SurfaceMesh* pm_A = vis();
    FaceData<int> selected(m,0);
    polyscope::state::userCallback = [&]() {
        // get the mouse location from ImGui
        std::size_t select_idx;
        if (selectFace(pm_A,&select_idx)) {
            int& i = selected[m.face(select_idx)];
            i = (i+1) % 2;
            pm_A->addFaceScalarQuantity("selected",selected,polyscope::DataType::CATEGORICAL)->setMapRange({0,1})->setEnabled(true);
        }
        if (ImGui::Button("Select all faces")) {
            selected.fill(1);
            pm_A = vis();
        }
        if (ImGui::Button("Refine selected faces")) {
            AdaptiveFaceTransfer transfer(m,g,u,atri.marked_corner);
            atri.refine(selectedFaces(selected),&transfer);
            atri.coarse({},&transfer);
            u = transfer.transfer();
            selected.fill(0);
            pm_A = vis();
        }
        if (ImGui::Button("Coarse selected faces")) {
            AdaptiveFaceTransfer transfer(m,g,u,atri.marked_corner);
            atri.refine({},&transfer);
            atri.coarse(selectedFaces(selected),&transfer);
            u = transfer.transfer();
            selected.fill(0);
            pm_A = vis();
        }
    };

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
        ASSERT_NEAR(fB[v],int_positions_A[v].x,0.00001);
    }
}



