#include <gtest/gtest.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h>
#include <filesystem>
#include <eider/homology.h>

#include <chrono>

#include "polyscope/curve_network.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;
void updateTriagulationViz(SurfaceMesh& mesh, VertexPositionGeometry& geometry, EdgeData<EdgeType> t) {

    // Convert to 3D positions
    std::vector<std::vector<Vector3>> traces;
    std::vector<std::size_t> bridge_edges;
    std::vector<double> edge_color;
    for (Edge e : mesh.edges())
    {
        if (t[e] == EdgeType::minimal_st) {
            auto& line = traces.emplace_back();
            line.push_back(geometry.vertexPositions[e.firstVertex()]);
            line.push_back(geometry.vertexPositions[e.secondVertex()]);
            edge_color.push_back(0.0);
        } else if (t[e] == EdgeType::bridge) {
            auto& line = traces.emplace_back();
            line.push_back(geometry.vertexPositions[e.firstVertex()]);
            line.push_back(geometry.vertexPositions[e.secondVertex()]);
            edge_color.push_back(2.0);
        } else if (t[e] == EdgeType::maximal_co_st) {
            auto& line = traces.emplace_back();
            Vector3 v13 = Vector3(1./3,1./3,1./3);
            line.push_back(SurfacePoint(e.halfedge().face(),v13).interpolate(geometry.vertexPositions));
            line.push_back(SurfacePoint(e,0.5).interpolate(geometry.vertexPositions));
            line.push_back(SurfacePoint(e.halfedge().twin().face(),v13).interpolate(geometry.vertexPositions));
            edge_color.push_back(1.0);
            edge_color.push_back(1.0);
        }
    }
    std::vector<Vector3> tracesPts;
    std::vector<std::array<size_t, 2>> tracesEdgeInds;
    for (std::vector<Vector3>& line : traces) {
        if (line.size() < 2) continue;
        tracesPts.push_back(line[0]);
        for (size_t i = 0; i < line.size() - 1; i++) {
            tracesPts.push_back(line[i + 1]);
            tracesEdgeInds.push_back({tracesPts.size() - 2, tracesPts.size() - 1});
        }
    }

    // Register with polyscope
    auto psCurves = polyscope::registerCurveNetwork("intrinsic edges", tracesPts, tracesEdgeInds);
    psCurves->addEdgeScalarQuantity("Bridge Edge", edge_color)->setEnabled(true);
    psCurves->setRadius(0.01);
    psCurves->setEnabled(true);
}

TEST(homologyTest, TestHomotopyBasis)
{
    using namespace geometrycentral::surface;
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" /"two_torus.stl";
    auto [m,g] = readManifoldSurfaceMesh(fds.string());
    EdgeData<EdgeType> edge_data(*m, EdgeType::bridge);
    computeMinimalSpanningTree(*m,*g,edge_data);
    computePrimalEdgesOfDualMaxST(*m,*g, edge_data);
    auto d_edges = distinctEdges(*m,edge_data);
    Face x = m->face(0);
    auto dijkstra = co_dijkstra(*m,*g,edge_data,x);
    FaceData<int> prev(*m, 0); for (Face f: m->faces()) {
        Halfedge he = dijkstra.first[f];
        if (he != Halfedge()) { prev[f] = he.face().getIndex(); } else { prev[f] = -1; }
    }



    g->requireDualEdgeLengths();
    polyscope::init();
    polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh("M", g->vertexPositions,m->getFaceVertexList(), polyscopePermutations(*m));
    pm->addEdgeScalarQuantity("minimal - max_co - bridge",edge_data,polyscope::DataType::CATEGORICAL);
    pm->addEdgeScalarQuantity("edge_lenghts",g->edgeLengths);
    pm->addEdgeScalarQuantity("dual edge lengths",g->dualEdgeLengths);
    pm->addFaceScalarQuantity("dijkstra distances",dijkstra.second);
    pm->addFaceScalarQuantity("dijkstra prev",prev);

    updateTriagulationViz(*m, *g, edge_data);
    auto h_basis = homotopy_basis(*m,*g,x);
    int basis_i = 0;
    for (const auto& basis: h_basis) {
        EdgeData<int> b(*m, 0);
        for (Halfedge e: basis) { b[e.edge()] = 1;}
        pm->addEdgeScalarQuantity("homology basis " + std::to_string(basis_i),b,polyscope::DataType::CATEGORICAL);
        basis_i++;
    }
    polyscope::show();
}

TEST(homologyTest, TestDeRhamCohom)
{
    using namespace geometrycentral::surface;
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" /"torus.stl";
    auto [m,g] = readManifoldSurfaceMesh(fds.string());
    Face x = m->face(0);
    auto h_basis = homotopy_basis(*m,*g,x);

    polyscope::init();
    polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh("M", g->vertexPositions,m->getFaceVertexList(), polyscopePermutations(*m));
    int basis_i = 0;
    g->requireDECOperators();

    g->requireFaceTangentBasis();
    FaceData<Vector3> e1(*m),e2(*m);
    for (Face f: m->faces()) { e1[f] = g->faceTangentBasis[f][0], e2[f] = g->faceTangentBasis[f][1]; }

    for (const auto& basis: h_basis) {
        EdgeData<int> b(*m, 0);
        for (Halfedge e: basis) { b[e.edge()] = 1;}
        auto df =delta_form(*m, basis);
        pm->addEdgeScalarQuantity("delta form " + std::to_string(basis_i), df);
        Eigen::VectorXd l  = g->d1 * df.raw();
        pm->addFaceScalarQuantity("delta form ex der " + std::to_string(basis_i),FaceData<double>(*m,l)) ;
        EdgeData<double> pf = pressure_project(*m, df, *g);
        FaceData<double> dpf = FaceData<double>(*m, g->d1 * pf.raw());
        pm->addEdgeScalarQuantity("pressure projection " + std::to_string(basis_i), pf);
        pm->addFaceScalarQuantity("pressure projection ex der " + std::to_string(basis_i), dpf);
        pm->addEdgeScalarQuantity("homology basis " + std::to_string(basis_i),b,polyscope::DataType::CATEGORICAL);
        FaceData<Vector2> wi = whitney_interpolation(*m,*g,pf);
        pm->addFaceTangentVectorQuantity("homology basis Whitney " + std::to_string(basis_i), wi, e1,e2);
        basis_i++;
    }
    polyscope::show();
}

TEST(homologyTest, TestWhitney)
{
    using namespace geometrycentral::surface;
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" /"torus.stl";
    auto [m,g] = readManifoldSurfaceMesh(fds.string());
    g->requireHalfedgeVectorsInFace();
    g->requireFaceTangentBasis();
    auto basis = orthonormal_hom_basis(*m,*g);

    g->requireFaceTangentBasis();
    FaceData<Vector3> e1(*m),e2(*m);
    for (Face f: m->faces()) { e1[f] = g->faceTangentBasis[f][0], e2[f] = g->faceTangentBasis[f][1]; }

    polyscope::init();
    polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh("M", g->vertexPositions,m->getFaceVertexList(), polyscopePermutations(*m));
    std::size_t i = 0;
    for (const auto& b: basis) {
        pm->addFaceTangentVectorQuantity("Hom basis" + std::to_string(i),b,e1,e2);
        i++;
    }
    polyscope::show();
}
