#include <gtest/gtest.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h>
#include <filesystem>
#include <eider/refine.h>

#include "polyscope/curve_network.h"
#include <chrono>

#include "Stopwatch.h"
#include "eider/poisson.h"
#include <implot.h>

#include "eider/util.h"

using namespace geometrycentral::surface;
using namespace geometrycentral;

void updateTriagulationViz(IntrinsicTriangulation& signpostTri, VertexPositionGeometry& geometry) {
    using namespace geometrycentral::surface;
    using namespace geometrycentral;

    // Get the edge traces
    EdgeData<std::vector<SurfacePoint>> traces = signpostTri.traceAllIntrinsicEdgesAlongInput();

    // Convert to 3D positions
    std::vector<std::vector<Vector3>> traces3D(traces.size());
    size_t i = 0;
    for (Edge e : signpostTri.mesh.edges()) {
        for (SurfacePoint& p : traces[e]) {
            traces3D[i].push_back(p.interpolate(geometry.inputVertexPositions));
        }
        i++;
    }
    std::vector<Vector3> tracesPts;
    std::vector<std::array<size_t, 2>> tracesEdgeInds;
    for (std::vector<Vector3>& line : traces3D) {
        if (line.size() < 2) continue;
        tracesPts.push_back(line[0]);
        for (size_t i = 0; i < line.size() - 1; i++) {
            tracesPts.push_back(line[i + 1]);
            tracesEdgeInds.push_back({tracesPts.size() - 2, tracesPts.size() - 1});
        }
    }


    // Register with polyscope
    auto psCurves = polyscope::registerCurveNetwork("intrinsic edges", tracesPts, tracesEdgeInds);
    psCurves->setRadius(0.001);
    psCurves->setEnabled(true);
}


void uniform_refine(AdaptiveTriangulation& intrT) {
    auto v = std::vector<Face>();
    v.reserve(intrT.mesh().nFaces());
    for (auto f: intrT.mesh().faces()) {
        v.push_back(f);
    }
    intrT.refine(v);

}

double refinement_ms = 0;

TEST(refineTest, testSplit)
{
    using namespace geometrycentral;
    using namespace geometrycentral::surface;
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" / "quad.stl";
    auto [parent_m,parent_g] = readManifoldSurfaceMesh(fds.string());

    IntegerCoordinatesIntrinsicTriangulation icit(*parent_m,*parent_g);
    // icit.flipToDelaunay();

    AdaptiveTriangulation atri(icit);

    ManifoldSurfaceMesh& m = atri.mesh();
    FaceData<int> selection(m,0);
    FaceData<double> residual(m,0);

    polyscope::SurfaceMesh* pm;
    std::vector<double> n_tri, error, h1, h2;
    auto vis_mg = [&m,&icit,&parent_g,&selection, &pm,&n_tri, &error, &h1, &h2,&residual, &atri]() {
        m.compress(); icit.refreshQuantities();
        VertexData<Vector3> int_positions(m) ;
        for (Vertex v : m.vertices()) {
            int_positions[v] = icit.vertexLocations[v].interpolate(parent_g->vertexPositions);
        }
        pm = polyscope::registerSurfaceMesh("M", int_positions,m.getFaceVertexList());
        pm->setAllPermutations(polyscopePermutations(m));
        pm->setSelectionMode(polyscope::MeshSelectionMode::FacesOnly);
        pm->addFaceScalarQuantity("selected faces", selection);
        HalfedgeData<int> he_d(atri.mesh(),0);
        for (Face f: atri.mesh().faces()) {
            he_d[atri.getRefinementEdge(f)] = 1;
        }
        pm->addHalfedgeScalarQuantity("refinement edges", he_d);


        VertexData<double> f(m,1);
        VertexData<double> u(m,0);

        StreamFunctionSolver S {};
        S.compute(m, icit);
        S.solve(m,icit,u,f);
        residual = poisson_residual_error_sqr(m, icit, u, f);
        double s = 0; for (Face f: m.faces()) { s += residual[f]; }
        pm->addVertexScalarQuantity("u",u);
        pm->addVertexScalarQuantity("f",f);
        pm->addFaceScalarQuantity("residual error",residual);
        n_tri.push_back(m.nFaces());
        error.push_back(s);
        h1.push_back(max_diameter(m,icit));
        h2.push_back(std::cbrt(h1.back()*h1.back()));
    };


    polyscope::init();
    ImPlot::CreateContext();
    vis_mg();

    float theta = 0.6;


    polyscope::state::userCallback = [&]()
    {

        ImGuiIO& io = ImGui::GetIO();
        if (io.MouseClicked[0] && !io.WantCaptureMouse)
        {
            glm::vec2 screenCoords{io.MousePos.x, io.MousePos.y};
            polyscope::PickResult pickResult = polyscope::pickAtScreenCoords(screenCoords);

            // get additional information if we clicked on a mesh
            if(pickResult.isHit && pickResult.structure == pm) {
                polyscope::SurfaceMeshPickResult meshPickResult =
                  pm->interpretPickResult(pickResult);

                if(meshPickResult.elementType == polyscope::MeshElement::FACE) {
                    std::cout << "clicked face " << meshPickResult.index << std::endl;
                    Face f = m.face(meshPickResult.index);
                    if (selection[f] == 1) selection[f] =0;
                    else selection[f] = 1;
                    pm->addFaceScalarQuantity("selected faces", selection)->setEnabled(true);
                }
            }
        }

        ImGui::DragFloat("etha",&theta,0.025,0,1);
        ImGui::SameLine();
        if (ImGui::Button("Select doerfler"))
        {
            auto faces = select_doerfler(m, residual, theta, 1);
            selection.fill(0);
            for (Face f:faces) {
                selection[f] = 1;
            }
            vis_mg();
        }
        // Build a UI element to edit a parameter, which will
        // appear in the onscreen panel
        if (ImGui::Button("Refine Selected"))
        {
            std::vector<Face> faces;
            for (Face f:m.faces()) {
                if (selection[f] == 1) faces.push_back(f);
            }
            atri.refine(faces);
            selection.fill(0);
            vis_mg();
        }
        ImGui::SameLine();
        if (ImGui::Button("Uniform Refine"))
        {
            {
                Stopwatch sw(refinement_ms);
                uniform_refine(atri);
            }
            vis_mg();
        }
        ImGui::Text(("duration: "+ std::to_string(refinement_ms) + " ms").c_str());

        if (ImPlot::BeginPlot("My Plot")) {
            ImPlot::SetupAxis(ImAxis_X1, "#Triangles",ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxis(ImAxis_Y1, "residual error estimate", ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxis(ImAxis_Y2, "Optimal", ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_Opposite);
            // ImPlot::SetupAxisScale(ImAxis_X1);

            std::size_t n = n_tri.size();
            std::vector<double> bound1(n), bound2(n);
            for (int i = 0; i < n_tri.size(); ++i) {
                bound1[i] = 1. / std::sqrt(double(n_tri[i]));
                bound2[i] = 1. / std::cbrt(double(n_tri[i]));
            }

            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::SetAxes(ImAxis_X1, ImAxis_Y1);
            ImPlot::PlotLine("residual error estimate", n_tri.data(), error.data(),n_tri.size());
            ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
            ImPlot::PlotLine("#T^-(1/2)", n_tri.data(), bound1.data(),n_tri.size());
            ImPlot::PlotLine("#T^-(1/3)", n_tri.data(), bound2.data(),n_tri.size());
            ImPlot::PlotLine("h", n_tri.data(), h1.data(),n_tri.size());
            ImPlot::PlotLine("h^(2/3)", n_tri.data(), h2.data(),n_tri.size());
            ImPlot::EndPlot();
        }
        ImPlot::ShowDemoWindow();

    };// specify the callback
    polyscope::show();
    ImPlot::DestroyContext();
}

TEST(refineTest, perf_test)
{
    using namespace geometrycentral::surface;
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/"models"/"L.stl";
    auto [m,g] = readManifoldSurfaceMesh(fds.string());

    IntegerCoordinatesIntrinsicTriangulation icit(*m,*g);
    AdaptiveTriangulation atri(icit);
    for (int i = 0; i< 10; i++)
    {
        uniform_refine(atri);
    }
}

TEST(afemTest, estimate)
{
    using namespace geometrycentral::surface;
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" / "L_fine.stl";
    auto [m,g] = readManifoldSurfaceMesh(fds.string());

    VertexData<double> f(*m,1);
    VertexData<double> u(*m,0);

    StreamFunctionSolver S {};
    S.compute(*m, *g);
    S.solve(*m,*g,u,f);
    auto e = poisson_residual_error_sqr(*m, *g, u, f);

    polyscope::init();
    polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh("M", g->vertexPositions,m->getFaceVertexList());
    pm->addFaceScalarQuantity("residual error",e)->setEnabled(true);
    pm->addVertexScalarQuantity("u",u);
    pm->addVertexScalarQuantity("f",f);
    polyscope::show();
}

double discrete_energy(IntrinsicGeometryInterface& T, VertexData<double> U)
{
    T.requireCotanLaplacian();
    // double d=0; for (Face f: T.mesh.faces()) d+=T.faceAreas[f]*grad(T,f,U).norm2();
    // return d;
    Eigen::VectorXd x = U.toVector();
    double d = x.transpose() * (T.cotanLaplacian) * x;
    T.unrequireCotanLaplacian();
    return d;
}

double du_norm_sqr() { return 1.06422; }

double abs_error(IntrinsicGeometryInterface& T, VertexData<double> U)
{
    return std::sqrt(du_norm_sqr() - discrete_energy(T,U));
}
double rel_error(IntrinsicGeometryInterface& T, VertexData<double> U)
{
    return std::sqrt(du_norm_sqr() - discrete_energy(T,U))/std::sqrt(du_norm_sqr());
}

TEST(refineTest, uniform_vs_adaptive_refinement) {
    float theta = 0.5;
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" / "L_coarse.stl";
    auto [parent_m,parent_g] = readManifoldSurfaceMesh(fds.string());

    IntegerCoordinatesIntrinsicTriangulation uniform_g(*parent_m,*parent_g), adaptive_g(*parent_m, *parent_g);
    AdaptiveTriangulation uniform_a(uniform_g), adaptive_a(adaptive_g);
    ManifoldSurfaceMesh &uniform_m = *uniform_g.intrinsicMesh, &adaptive_m = *adaptive_g.intrinsicMesh;
    while (uniform_m.nFaces() < 3072) {uniform_refine(uniform_a); uniform_refine(adaptive_a); }
    FaceData<double> res_u(uniform_m,0), res_a(adaptive_m,0);
    struct Data
    {
        std::vector<double> abs_error, rel_error, residual, nTri;
    };
    Data data_uniform, data_adaptive;


    auto vis_mg = [&](std::string name, IntegerCoordinatesIntrinsicTriangulation& icit, Data& data, FaceData<double>& residual) {
        ManifoldSurfaceMesh& m = *icit.intrinsicMesh;
        m.compress();
        VertexData<Vector3> int_positions(m) ;
        for (Vertex v : m.vertices()) {
            int_positions[v] = icit.vertexLocations[v].interpolate(parent_g->vertexPositions);
        }
        polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh(name, int_positions,m.getFaceVertexList());

        VertexData<double> f(m,1);
        VertexData<double> u(m,0);

        StreamFunctionSolver S {};
        S.compute(m, icit);
        S.solve(m,icit,u,f);
        pm->addVertexScalarQuantity(name + " - u",u);
        pm->addVertexScalarQuantity(name + " - f",f);
        residual = poisson_residual_error_sqr(m, icit, u, f);
        double s = 0; for (Face f: m.faces()) { s += residual[f]; }
        auto* sq = pm->addFaceScalarQuantity(name + " - residual error",residual);
        data.nTri.push_back(m.nFaces());
        data.abs_error.push_back(abs_error(icit,u));
        data.residual.push_back(s);
        data.rel_error.push_back(rel_error(icit,u));
    };
    auto vis_all = [&]() { vis_mg("uniform", uniform_g,data_uniform,res_u); vis_mg("adaptive", adaptive_g,data_adaptive,res_a); };
    auto adaptive_refine = [&]() {
        VertexData<double> f(adaptive_m,1); VertexData<double> u(adaptive_m,0);
        StreamFunctionSolver S {}; S.compute(adaptive_m, adaptive_g); S.solve(adaptive_m,adaptive_g,u,f);

        res_a = poisson_residual_error_sqr(adaptive_m, adaptive_g, u, f);
        auto faces = select_doerfler(adaptive_m, res_a, theta, 1);
        adaptive_a.refine(faces);
    };

    polyscope::init();
    vis_all();
    ImPlot::CreateContext();

    polyscope::state::userCallback = [&]()
    {
        ImGui::DragFloat("theta",&theta,0.01,0,1);
        if (ImGui::Button("Refine until 2xT"))
        {
            adaptive_refine();
            while (uniform_m.nFaces()*2 < adaptive_m.nFaces()) uniform_refine(uniform_a);
            vis_all();
        }

        if (ImPlot::BeginPlot("Absolut error")) {
            ImPlot::SetupAxis(ImAxis_X1, "#Triangles",ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxis(ImAxis_Y1, "absolut error", ImPlotAxisFlags_AutoFit);
            // ImPlot::SetupAxis(ImAxis_Y2, "Optimal", ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_Opposite);
            ImPlot::SetupAxisScale(ImAxis_Y1,ImPlotScale_Log10);
            ImPlot::SetupAxisScale(ImAxis_X1,ImPlotScale_Log10);

            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::SetAxes(ImAxis_X1, ImAxis_Y1);
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::PlotLine("uniform", data_uniform.nTri.data(), data_uniform.abs_error.data(),data_uniform.nTri.size());
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::PlotLine("adaptive", data_adaptive.nTri.data(), data_adaptive.abs_error.data(),data_adaptive.nTri.size());
            ImPlot::EndPlot();
        }
        if (ImPlot::BeginPlot("Relative Error")) {
            ImPlot::SetupAxis(ImAxis_X1, "#Triangles",ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxis(ImAxis_Y1, "relative Error", ImPlotAxisFlags_AutoFit);
            // ImPlot::SetupAxis(ImAxis_Y2, "Optimal", ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_Opposite);
            ImPlot::SetupAxisScale(ImAxis_Y1,ImPlotScale_Log10);
            ImPlot::SetupAxisScale(ImAxis_X1,ImPlotScale_Log10);

            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::SetAxes(ImAxis_X1, ImAxis_Y1);
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::PlotLine("uniform", data_uniform.nTri.data(), data_uniform.rel_error.data(),data_uniform.nTri.size());
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::PlotLine("adaptive", data_adaptive.nTri.data(), data_adaptive.rel_error.data(),data_adaptive.nTri.size());
            ImPlot::EndPlot();
        }
        if (ImPlot::BeginPlot("Residual")) {
            ImPlot::SetupAxis(ImAxis_X1, "#Triangles",ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxis(ImAxis_Y1, "residual error estimate", ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1,ImPlotScale_Log10);
            ImPlot::SetupAxisScale(ImAxis_X1,ImPlotScale_Log10);

            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::PlotLine("uniform residual", data_uniform.nTri.data(), data_uniform.residual.data(),data_uniform.nTri.size());
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::PlotLine("adaptive residual", data_adaptive.nTri.data(), data_adaptive.residual.data(),data_adaptive.nTri.size());
            ImPlot::EndPlot();
        }
    };// specify the callback
    polyscope::show();
    ImPlot::DestroyContext();


}

TEST(refineTest,testCollapsingSameElement) {
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" / "quad.stl";
    auto [parent_m,parent_g] = readManifoldSurfaceMesh(fds.string());
    IntegerCoordinatesIntrinsicTriangulation icit(*parent_m,*parent_g);
    AdaptiveTriangulation atri(icit);
    ManifoldSurfaceMesh& m = atri.mesh();

    Edge c_edge; std::array<Edge,4> b_edges;
    int i = 0; for (Edge e: m.edges()) if (e.isBoundary()) b_edges[i++] = e; else c_edge = e;

    Halfedge he = atri.vertex_bisection(c_edge.halfedge());
    Halfedge phe = he.prevOrbitFace().twin().prevOrbitFace();
    i = 0; for (Edge e: he.vertex().adjacentEdges()) b_edges[i++] = e;
    Vertex v = atri.vertex_biunion(he).vertex();
    Edge a_edge;
    i = 0; for (Edge e: b_edges) if (!e.isDead()) a_edge = e;

    ASSERT_TRUE(he.isDead());
    ASSERT_FALSE(phe.edge().isDead());
    ASSERT_FALSE(a_edge.isDead());
    ASSERT_EQ(phe.edge(),a_edge);

    EdgeData<double> d(m,0);
    d[a_edge] = 1;

    icit.refreshQuantities(); m.compress();

    VertexData<Vector3> int_positions(m) ;
    for (Vertex v : m.vertices()) {
        int_positions[v] = icit.vertexLocations[v].interpolate(parent_g->vertexPositions);
    }

    polyscope::init();
    polyscope::SurfaceMesh* pm_A = polyscope::registerSurfaceMesh("original", parent_g->vertexPositions,parent_m->getFaceVertexList());
    polyscope::SurfaceMesh* pm_B = polyscope::registerSurfaceMesh("Intrinsic", int_positions,m.getFaceVertexList());
    pm_B->setAllPermutations(polyscopePermutations(m));
    pm_A->addFaceScalarQuantity("index",parent_m->getFaceIndices())->setEnabled(true);
    pm_B->addFaceScalarQuantity("index",m.getFaceIndices())->setEnabled(true);
    pm_B->addEdgeScalarQuantity("d",d)->setEnabled(true);
    pm_A->setEdgeWidth(1);
    pm_B->setEdgeWidth(1);
    pm_A->translate(glm::vec3(-1,0,0));
    pm_B->translate(glm::vec3( 1,0,0));
    polyscope::show();
}

TEST(refineTest,testCoarsingCondition) {
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" / "quad.stl";
    auto [parent_m,parent_g] = readManifoldSurfaceMesh(fds.string());
    IntegerCoordinatesIntrinsicTriangulation icit(*parent_m,*parent_g);
    AdaptiveTriangulation atri(icit);
    ManifoldSurfaceMesh& m = *icit.intrinsicMesh;

    Edge c_edge; std::array<Edge,4> b_edges;
    int i = 0; for (Edge e: m.edges()) if (e.isBoundary()) b_edges[i++] = e; else c_edge = e;

    Halfedge he = atri.vertex_bisection(c_edge.halfedge());
    ASSERT_EQ(atri.mesh().nFaces(),4);
    ASSERT_FALSE(he.vertex().isBoundary());
    Face f1 = he.face(), f2 = he.prevOrbitFace().twin().face();
    Face f3 = he.twin().next().twin().face(), f4 = he.twin().face();
    ASSERT_GT(atri.idx[f1],atri.idx[f2]);
    ASSERT_GT(atri.idx[f3],atri.idx[f4]);
    he = atri.vertex_bisection(he.next());
    ASSERT_TRUE(he.isInterior());
    ASSERT_NE(he,Halfedge());
    auto nhe = atri.coarse_halfedge(he.vertex());
    ASSERT_EQ(nhe, he);
    Halfedge v = atri.vertex_biunion(he);
    ASSERT_NE(v, Halfedge());

    icit.refreshQuantities(); m.compress();


    VertexData<Vector3> int_positions(m) ;
    for (Vertex v : m.vertices()) {
        int_positions[v] = icit.vertexLocations[v].interpolate(parent_g->vertexPositions);
    }

    polyscope::init();
    polyscope::SurfaceMesh* pm_A = polyscope::registerSurfaceMesh("original", parent_g->vertexPositions,parent_m->getFaceVertexList());
    polyscope::SurfaceMesh* pm_B = polyscope::registerSurfaceMesh("Intrinsic", int_positions,m.getFaceVertexList());
    pm_B->setAllPermutations(polyscopePermutations(m));
    FaceData<std::size_t> face_idx(m);
    for (Face f:m.faces()) { face_idx[f] = atri.idx[f];}
    pm_B->addFaceScalarQuantity("index",face_idx)->setEnabled(true);
    pm_A->setEdgeWidth(1);
    pm_B->setEdgeWidth(1);
    pm_A->translate(glm::vec3(-1,0,0));
    pm_B->translate(glm::vec3( 1,0,0));
    polyscope::show();
}

TEST(refineTest,testCoarseRefine) {
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" / "nvb-coarsening.stl";
    auto [parent_m,parent_g] = readManifoldSurfaceMesh(fds.string());

    IntegerCoordinatesIntrinsicTriangulation icit(*parent_m,*parent_g);
    AdaptiveTriangulation atri(icit);
    ManifoldSurfaceMesh& m = atri.mesh();
    IntrinsicGeometryInterface& g = atri.geom();

    Edge c_edge; std::array<Edge,4> b_edges;
    // int i = 0; for (Edge e: m.edges()) if (e.isBoundary()) b_edges[i++] = e; else c_edge = e;
    c_edge = m.edge(0);

    Halfedge he1 = atri.vertex_bisection(c_edge.halfedge());
    Halfedge he2 = atri.vertex_bisection(he1.next());
    std::vector<Face> faces { }; for (Face f: m.faces()) faces.push_back(f);
    atri.coarse(faces);
    icit.refreshQuantities(); m.compress();

    VertexData<Vector3> int_positions(m) ;
    for (Vertex v : m.vertices()) {
        int_positions[v] = icit.vertexLocations[v].interpolate(parent_g->vertexPositions);
    }

    polyscope::init();
    polyscope::SurfaceMesh* pm_A = polyscope::registerSurfaceMesh("original", parent_g->vertexPositions,parent_m->getFaceVertexList());
    polyscope::SurfaceMesh* pm_B = polyscope::registerSurfaceMesh("Intrinsic", int_positions,m.getFaceVertexList());
    pm_A->addFaceScalarQuantity("index",parent_m->getFaceIndices())->setEnabled(true);
    pm_B->addFaceScalarQuantity("index",m.getFaceIndices())->setEnabled(true);
    pm_A->setEdgeWidth(1);
    pm_B->setEdgeWidth(1);
    pm_A->translate(glm::vec3(-1,0,0));
    pm_B->translate(glm::vec3( 1,0,0));
    polyscope::show();
}
