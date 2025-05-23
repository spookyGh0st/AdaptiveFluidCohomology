#include <gtest/gtest.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h>
#include <filesystem>
#include <eider/homology.h>

#include <chrono>

TEST(homologyTest, perf_test)
{
    using namespace geometrycentral::surface;
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" /"torus_deformed.stl";
    auto [m,g] = readManifoldSurfaceMesh(fds.string());
    auto mst = computeMinimalSpanningTree(*m,*g);
    auto d_mst = computePrimalEdgesOfDualMaxST(*m,*g);
    EdgeData<int> eData1(*m,0);
    EdgeData<int> eData2(*m,0);
    for (Edge e: mst) { eData1[e] = 1; }
    for (Edge e: d_mst) {  eData2[e] = 1; }
    g->requireDualEdgeLengths();
    polyscope::init();
    polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh("M", g->vertexPositions,m->getFaceVertexList(), polyscopePermutations(*m));
    pm->addEdgeScalarQuantity("msp",eData1,polyscope::DataType::CATEGORICAL);
    pm->addEdgeScalarQuantity("co_msp",eData2,polyscope::DataType::CATEGORICAL);
    pm->addEdgeScalarQuantity("edge_lenghts",g->edgeLengths);
    pm->addEdgeScalarQuantity("dual edge lengths",g->dualEdgeLengths);
    polyscope::show();
}
