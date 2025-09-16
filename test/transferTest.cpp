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
    atri.coarse(faces);
    ASSERT_TRUE(center_v.isDead());
    ASSERT_EQ(idx[center_v],1);
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

    AdaptiveTransfer transfer(atri.intrinsicTriangulation(),fA);
    atri.refine(facesR,&transfer);
    atri.coarse(facesC,&transfer);
    VertexData<double> fB = transfer.transfer();

    int_positions_A = intrinsic_position(atri.intrinsicTriangulation(),*parent_g);
    for(Vertex v: m.vertices()){
        ASSERT_NEAR(fB[v],int_positions_A[v].x,0.000001);
    }
}



