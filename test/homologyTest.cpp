#include <gtest/gtest.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h>
#include <filesystem>
#include <eider/homology.h>

#include <chrono>

#include "Stopwatch.h"
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
            if (e.halfedge().isInterior())
                line.push_back(SurfacePoint(e.halfedge().face(),v13).interpolate(geometry.vertexPositions));
            line.push_back(0.5* (geometry.vertexPositions[e.firstVertex()] + geometry.vertexPositions[e.secondVertex()]));
            if (e.halfedge().twin().isInterior())
                line.push_back(SurfacePoint(e.halfedge().twin().face(),v13).interpolate(geometry.vertexPositions));
            for (int i = 0; i < line.size()-1; ++i) {
                edge_color.push_back(1.0);
            }
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
    fds = fds.parent_path()/ "models" /"grid_hole.stl";
    auto [m,g] = readManifoldSurfaceMesh(fds.string());
    EdgeData<EdgeType> edge_data(*m, EdgeType::bridge);
    g->requireDualEdgeLengths();
    computePrimalEdgesOfDualMaxST(*m,edge_data, [&l = g->dualEdgeLengths] (Edge a, Edge b) {return l[a] < l[b]; });
    computeMinimalSpanningTree(*m,edge_data,[&l = g->edgeLengths] (Edge a, Edge b) {return l[a] > l[b]; });
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
        EdgeData<int> b(*m, 0), b_reduced(*m,0);
        auto red_basis = reduce_co_loop(*m,basis);
        for (Halfedge e: basis) { b[e.edge()] = 1;}
        for (Halfedge e: red_basis) { b_reduced[e.edge()] = 1;}
        pm->addEdgeScalarQuantity("homology basis " + std::to_string(basis_i),b,polyscope::DataType::CATEGORICAL);
        pm->addEdgeScalarQuantity("reduced homology basis " + std::to_string(basis_i),b_reduced,polyscope::DataType::CATEGORICAL);
        basis_i++;
    }
    polyscope::show();
}

TEST(homologyTest, TestDeRhamCohom)
{
    using namespace geometrycentral::surface;
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" /"grid_hole.stl";
    auto [m,g] = readManifoldSurfaceMesh(fds.string());
    Face x = m->face(0);
    auto h_basis = homotopy_basis(*m,*g,m->face(0));

    polyscope::init();
    polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh("M", g->vertexPositions,m->getFaceVertexList(), polyscopePermutations(*m));
    int basis_i = 0;
    g->requireDECOperators();

    g->requireFaceTangentBasis();
    FaceData<Vector3> e1(*m),e2(*m);
    for (Face f: m->faces()) { e1[f] = g->faceTangentBasis[f][0], e2[f] = g->faceTangentBasis[f][1]; }
    PressureProjectionSolver pp_solver {};
    pp_solver.compute(*g);

    for (const auto& basis: h_basis) {
        EdgeData<int> b(*m, 0);
        for (Halfedge e: basis) { b[e.edge()] = 1;}
        auto df =delta_form(*m, reduce_co_loop(*m, basis));
        // auto df =delta_form(*m, basis);
        pm->addEdgeScalarQuantity("delta form " + std::to_string(basis_i), df);
        Eigen::VectorXd l  = g->d1 * df.toVector();
        pm->addFaceScalarQuantity("delta form ex der " + std::to_string(basis_i),FaceData<double>(*m,l)) ;
        EdgeData<double> pf = pp_solver.solve(*m, df);
        FaceData<double> dpf = FaceData<double>(*m, g->d1 * pf.toVector());
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
    fds = fds.parent_path()/ "models" /"torus.obj";
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

TEST(homologyTest, Intrinsic)
{
    using namespace geometrycentral::surface;
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" /"grid_hole.stl";
    auto [m,g] = readManifoldSurfaceMesh(fds.string());
    IntegerCoordinatesIntrinsicTriangulation ig (*m,*g);
    ig.flipToDelaunay();
    g->requireFaceTangentBasis();

    ManifoldSurfaceMesh& im = *ig.intrinsicMesh;
    ig.requireHalfedgeVectorsInFace();
    auto basis = orthonormal_hom_basis(im,ig);

    VertexData<Vector3> vd(im);
    for (Vertex v: m->vertices()) { vd[im.vertex(v.getIndex())] = g->vertexPositions[v]; }
    VertexPositionGeometry vpg(im,vd);

    vpg.requireFaceTangentBasis();
    FaceData<Vector3> e1(im),e2(im);
    for (Face f: im.faces()) { e1[f] = vpg.faceTangentBasis[f][0], e2[f] = vpg.faceTangentBasis[f][1]; }

    polyscope::init();
    polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh("M", vpg.vertexPositions,im.getFaceVertexList(), polyscopePermutations(im));
    std::size_t i = 0;
    for (const auto& b: basis) {
        pm->addFaceTangentVectorQuantity("Hom basis" + std::to_string(i),b,e1,e2);
        i++;
    }
    polyscope::show();
}

TEST(homologyTest, TestPerformance)
{
    using namespace geometrycentral::surface;
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" /"torus_max.stl";
    auto [m,g] = readManifoldSurfaceMesh(fds.string());
    g->requireHalfedgeVectorsInFace();
    g->requireFaceTangentBasis();
    double time = 0;
    std::vector<FaceData<Vector2>> basis;
    {
        auto s = Stopwatch(time);
        basis = orthonormal_hom_basis(*m, *g);
    }
    std::cout<< "#faces:  " << m->nFaces() << std::endl;
    std::cout<< "#edges:  " << m->nEdges() << std::endl;
    std::cout<< "#vertices:  " << m->nVertices() << std::endl;
    std::cout<< "1betti:  " << basis.size() << std::endl;
    std::cout<< "runtime: " << time << " ms" << std::endl;
}

TEST(homologyTest, testProjection)
{
    using namespace geometrycentral::surface;
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" /"torus_bounded_max.stl";
    auto [m,g] = readManifoldSurfaceMesh(fds.string());
    PressureProjectionSolver P {};
    P.compute(*g);
    EdgeData<double> x(*m, Eigen::VectorXd::Random(m->nEdges()));
    EdgeData<double> px = P.solve(*m, x);
    Eigen::VectorXd dTx = (g->d0.transpose() * px.toVector());
    ASSERT_LE(dTx.norm(), 1e-10) << "projection does not project on kernel of A^T";
    ASSERT_LE(dTx.maxCoeff(), 1e-10) << "projection does not project on kernel of A^T ";

    Eigen::VectorXd xPx = px.toVector() - x.toVector();
    double d = px.toVector().transpose() * xPx;
    ASSERT_LE(std::abs(d), 1e-10) << "projection is not orthogonal";
}
