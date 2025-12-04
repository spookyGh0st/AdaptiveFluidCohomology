#include <gtest/gtest.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h>
#include <filesystem>
#include <eider/homology.h>

#include <chrono>

#include "Stopwatch.h"
#include "eider/homotopy.h"
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

Halfedge add_boundary_edge(Face x, EdgeData<EdgeType> &data) {
    Halfedge h;
    for (Halfedge he : x.adjacentHalfedges()) {
        if (he.edge().isBoundary()) {
            h = he;
            break;
        }
    }
    if (h == Halfedge())
        return Halfedge();
    data[h.edge()] = EdgeType::maximal_co_st;
    return h;
}

// Compute for each edge e in G\T^*  the length of the shortest loop in T^* that contains e
EdgeData<double> edge_coloop_lengths(ManifoldSurfaceMesh &mesh, FaceData<double> &dist, const EdgeData<EdgeType> &data) {
    EdgeData<double> loop_length(mesh, NAN);
    for (Edge e : mesh.edges()) {
        if (data[e] == EdgeType::maximal_co_st)
            continue;
        double d = 0;
        for (Halfedge he : e.adjacentHalfedges())
            if (he.isInterior())
                d += dist[he.face()];
        loop_length[e] = d;
    }
    return loop_length;
}

TEST(homologyTest, TestHomotopyBasis)
{
    using namespace geometrycentral::surface;
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" /"cylinder.stl";
    auto [m,g] = readManifoldSurfaceMesh(fds.string());
    ManifoldSurfaceMesh& mesh = *m;
    VertexPositionGeometry& geom = *g;
    g->requireDualEdgeLengths();
    Face x = arbitrary_base_face(*m);

    EdgeData<EdgeType> edge_data(mesh, EdgeType::bridge);
    geom.requireEdgeLengths();
    geom.requireDualEdgeLengths();
    auto prev_dist = co_dijkstra(mesh, geom, edge_data, x, false);
    for (Face f : mesh.faces()) {
        if (f == x)
            continue;
        Edge e = prev_dist.first[f].edge();
        assert(!e.isBoundary());
        edge_data[e] = EdgeType::maximal_co_st;
    }
    Halfedge x_he = add_boundary_edge(x, edge_data);
    if (mesh.hasBoundary() && x_he == Halfedge()) throw std::runtime_error("Called with non-boundary face for mesh with boundary");
    auto coloop_lengths =
        edge_coloop_lengths(mesh, prev_dist.second, edge_data);
    computeMinimalSpanningTree( mesh, edge_data, [&l = coloop_lengths](Edge a, const Edge &b) { return l[a] > l[b]; });

    auto dist_edges = distinctEdges(mesh, edge_data);
    std::vector<std::vector<Halfedge>> co_loops;
    for (Edge e : dist_edges) {
        co_loops.push_back(homotopy_co_loop(prev_dist.first, x, e, x_he));
    }
    geom.unrequireDualEdgeLengths();

    auto dijkstra = co_dijkstra(*m, *g, edge_data, x, true);
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
    pm->addEdgeScalarQuantity("coloop lengths",coloop_lengths);
    pm->addFaceScalarQuantity("dijkstra distances",dijkstra.second);
    pm->addFaceScalarQuantity("dijkstra prev",prev);

    updateTriagulationViz(*m, *g, edge_data);
    int basis_i = 0;
    for (const auto& basis: co_loops) {
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
    fds = fds.parent_path()/ "models" /"cylinder.stl";
    auto [m,g] = readManifoldSurfaceMesh(fds.string());
    Face x = m->face(0);
    auto h_basis = greedy_homotopy_basis(*m,*g, arbitrary_base_face(*m));

    polyscope::init();
    polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh("M", g->vertexPositions,m->getFaceVertexList(), polyscopePermutations(*m));
    int basis_i = 0;
    g->requireDECOperators();

    g->requireFaceTangentBasis();
    FaceData<Vector3> e1(*m),e2(*m);
    for (Face f: m->faces()) { e1[f] = g->faceTangentBasis[f][0], e2[f] = g->faceTangentBasis[f][1]; }
    PressureProjectionSolver pp_solver {};
    pp_solver.compute(*m,*g);

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
    fds = fds.parent_path()/ "models" /"cylinder.stl";
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
    fds = fds.parent_path()/ "models" /"halftorus.stl";
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
    P.compute(*m,*g);
    EdgeData<double> x(*m, Eigen::VectorXd::Random(m->nEdges()));
    EdgeData<double> px = P.solve(*m, x);
    Eigen::VectorXd dTx = (g->d0.transpose() * px.toVector());
    ASSERT_LE(dTx.norm(), 1e-10) << "projection does not project on kernel of A^T";
    ASSERT_LE(dTx.maxCoeff(), 1e-10) << "projection does not project on kernel of A^T ";

    Eigen::VectorXd xPx = px.toVector() - x.toVector();
    double d = px.toVector().transpose() * xPx;
    ASSERT_LE(std::abs(d), 1e-10) << "projection is not orthogonal";
}

template <typename DerivedA, typename DerivedB>
::testing::AssertionResult EigenNear(
    const char* expr1,
    const char* expr2,
    const char* tol_expr,
    const Eigen::MatrixBase<DerivedA>& a,
    const Eigen::MatrixBase<DerivedB>& b,
    double tol)
{
    if (a.rows() != b.rows() || a.cols() != b.cols()) {
        return ::testing::AssertionFailure() << "Size mismatch: "
                                             << expr1 << " is " << a.rows() << "x" << a.cols()
                                             << ", " << expr2 << " is " << b.rows() << "x" << b.cols();
    }

    double max_diff = (a - b).cwiseAbs().maxCoeff();
    if (max_diff <= tol) {
        return ::testing::AssertionSuccess();
    } else {
        std::stringstream ss;
        ss << "Maximum absolute difference " << max_diff
           << " exceeds tolerance " << tol
           << "\n" << expr1 << " =\n" << a
           << "\n" << expr2 << " =\n" << b;
        return ::testing::AssertionFailure() << ss.str();
    }
}

// Macro for nicer syntax
#define EXPECT_EIGEN_NEAR(val1, val2, tol) \
EXPECT_TRUE(EigenNear(#val1, #val2, #tol, val1, val2, tol))

#define ASSERT_EIGEN_NEAR(val1, val2, tol) \
ASSERT_TRUE(EigenNear(#val1, #val2, #tol, val1, val2, tol))


void reportMatrix(std::string name, Eigen::MatrixXd A) {
    std::cout << name <<": " << A.rows() << "x" << A.cols() << std::endl;
    Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
    std::cout << "  Rank:       " << lu.rank() << std::endl;
    std::cout << "  Invertible: " << lu.isInvertible() << std::endl;
}

/// This test shows, That (I-AA^+) x is equivalent to solving the least square system, even if the matrixes are not invertible
TEST(homologyTest, d0Test) {
    using namespace geometrycentral::surface;
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" /"band.stl";
    auto [mesh,geom] = readManifoldSurfaceMesh(fds.string());
    geom->required0();
    geom->requireDECOperators();


    // Eigen::MatrixXd A = hodgeStar1Galerkin3Point(*mesh,*geom,geom->edgeIndices).transpose() * geom->d0 * geom->hodge0Inverse.transpose();
    Eigen::MatrixXd A = geom->d0;
    Eigen::MatrixXd AT = A.transpose();
    auto m = A.rows(), n = A.cols();
    ASSERT_EQ(m,mesh->nEdges());
    ASSERT_EQ(n,mesh->nVertices());


    std::cout <<std::fixed << std::setprecision(2);
    std::cout << "*_0^{-1}\n" << geom->hodge0Inverse.toDense() << std::endl;
    std::cout << "*_1\n" << hodgeStar1Galerkin3Point(*mesh,*geom,geom->edgeIndices).toDense() << std::endl;


    reportMatrix("I", Eigen::MatrixXd::Identity(m,m));
    reportMatrix("A", A);
    std::cout << "d0\n" << geom->d0.toDense() << std::endl;
    std::cout << "New d0\n" << A << std::endl;
    reportMatrix("A A^T", A * AT);
    reportMatrix("A^T A", AT * A);

    auto homotopy_b = greedy_homotopy_basis(*mesh,*geom,arbitrary_base_face(*mesh));
    auto cycle = Singular_Circle(*mesh, homotopy_b[0]);
    auto df = delta_form(*mesh,cycle);

    EXPECT_EIGEN_NEAR(geom->d1 * df.toVector(), Eigen::VectorXd::Zero(mesh->nFaces()),0.00001);

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const Eigen::VectorXd& S = svd.singularValues();
    Eigen::MatrixXd S_inv = Eigen::MatrixXd::Zero(S.size(),S.size());

    double tol = 1e-9;

    // Invert only non-zero singular values
    for (int i = 0; i < S.size(); ++i) {
        if (S(i) > tol) {
            S_inv(i, i) = 1.0 / S(i);
        }
    }
    Eigen::MatrixXd Ap = svd.matrixV() * S_inv * svd.matrixU().transpose();
    reportMatrix("A^+", Ap);
    int r = (S.array() > tol).count();
    Eigen::MatrixXd U_r = svd.matrixU().leftCols(r);  // keep only non-zero singular vectors
    Eigen::MatrixXd AAp = U_r * U_r.transpose();
    reportMatrix("A^+",AAp);
    ASSERT_EIGEN_NEAR(A*Ap,AAp,0.00001);

    Eigen::VectorXd b = df.toVector();

    // solve argmin ||Ax - b||
    Eigen::FullPivLU<Eigen::MatrixXd> solver(AT*A);
    Eigen::VectorXd x_ls = solver.solve(AT*b);

    Eigen::VectorXd b_hat = AAp * b;
    EXPECT_EIGEN_NEAR(A * x_ls, b_hat,0.00001);
    EXPECT_EIGEN_NEAR(AT *(b-A*x_ls), Eigen::VectorXd::Zero(n),0.00001);

    Eigen::VectorXd h = b - A * x_ls;
    Eigen::VectorXd h_hat = b - b_hat;

    EXPECT_EIGEN_NEAR(AT *h, Eigen::VectorXd::Zero(n),0.00001);
    EXPECT_EIGEN_NEAR(AT *h_hat, Eigen::VectorXd::Zero(n),0.00001);

    EXPECT_EIGEN_NEAR(geom->d1 * h, Eigen::VectorXd::Zero(mesh->nFaces()),0.00001);
    EXPECT_EIGEN_NEAR(geom->d1 *h_hat, Eigen::VectorXd::Zero(mesh->nFaces()),0.00001);
}
