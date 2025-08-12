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

TEST(transfertTest,testCoarseRefine) {
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" / "quad.stl";
    auto [parent_m,parent_g] = readManifoldSurfaceMesh(fds.string());

    IntegerCoordinatesIntrinsicTriangulation icit(*parent_m,*parent_g);
    ManifoldSurfaceMesh& m = *icit.intrinsicMesh;
    IntrinsicGeometryInterface& g = icit;

    Edge c_edge; std::array<Edge,4> b_edges;
    int i = 0; for (Edge e: m.edges()) if (e.isBoundary()) b_edges[i++] = e; else c_edge = e;

    Halfedge he = icit.splitEdge(c_edge.halfedge(),0.5);
    icit.intrinsicMesh->collapseEdgeTriangular(he.edge());
    m.compress(); icit.refreshQuantities();

    ASSERT_EQ(m.nEdges(), parent_m->nEdges());
    // icit.joinEdge(he.vertex());



    VertexData<Vector3> int_positions(m) ;
    for (Vertex v : m.vertices()) {
        int_positions[v] = icit.vertexLocations[v].interpolate(parent_g->vertexPositions);
    }

    polyscope::init();
    polyscope::SurfaceMesh* pm_A = polyscope::registerSurfaceMesh("original", parent_g->vertexPositions,parent_m->getFaceVertexList());
    polyscope::SurfaceMesh* pm_B = polyscope::registerSurfaceMesh("Intrinsic", int_positions,m.getFaceVertexList());
    pm_A->setEdgeWidth(1);
    pm_B->setEdgeWidth(1);
    pm_A->translate(glm::vec3(-1,0,0));
    pm_B->translate(glm::vec3( 1,0,0));
    polyscope::show();
}
