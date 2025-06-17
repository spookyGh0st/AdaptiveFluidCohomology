#include <gtest/gtest.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h>
#include <geometrycentral/surface/signpost_intrinsic_triangulation.h>
#include <filesystem>
#include <eider/refine.h>

#include "polyscope/curve_network.h"
#include <chrono>

#include "Stopwatch.h"
#include "eider/poisson.h"
#include <implot.h>

void updateTriagulationViz(gcs::IntrinsicTriangulation& signpostTri, gcs::VertexPositionGeometry& geometry) {
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


void uniform_refine(gcs::IntrinsicTriangulation& intrT) {
    auto v = std::vector<gcs::Face>();
    v.reserve(intrT.intrinsicMesh->nFaces());
    for (auto f: intrT.intrinsicMesh->faces()) {
        v.push_back(f);
    }
    refine(intrT,v);

}

double refinement_ms = 0;

TEST(refineTest, testSplit)
{
    using namespace geometrycentral;
    using namespace geometrycentral::surface;
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" / "L_coarse.stl";
    auto [parent_m,parent_g] = readManifoldSurfaceMesh(fds.string());

    IntegerCoordinatesIntrinsicTriangulation icit(*parent_m,*parent_g);
    icit.flipToDelaunay();
    // icit.delaunayRefine();

    ManifoldSurfaceMesh& m = *icit.intrinsicMesh;
    FaceData<int> selection(m,0);

    polyscope::SurfaceMesh* pm;
    std::vector<double> n_tri, error;
    auto vis_mg = [&m,&icit,&parent_g,&selection, &pm,&n_tri, &error]() {
        m.compress();
        VertexData<Vector3> int_positions(m) ;
        for (Vertex v : m.vertices()) {
            int_positions[v] = icit.vertexLocations[v].interpolate(parent_g->vertexPositions);
        }
        pm = polyscope::registerSurfaceMesh("M", int_positions,m.getFaceVertexList());
        pm->setSelectionMode(polyscope::MeshSelectionMode::FacesOnly);
        pm->addFaceScalarQuantity("selected faces", selection)->setEnabled(true);

        VertexData<double> f(m,1);
        VertexData<double> u(m,0);

        StreamFunctionSolver S {};
        S.compute(m, icit);
        S.solve(m,icit,u,f);
        auto e = poisson_residual_error(m, icit, f,u);
        double s = 0; for (Face f: m.faces()) { s += e[f]; }
        pm->addFaceScalarQuantity("residual error",e)->setEnabled(true);
        n_tri.push_back(m.nFaces());
        error.push_back(s);
    };


    polyscope::init();
    ImPlot::CreateContext();
    vis_mg();


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

        // Build a UI element to edit a parameter, which will
        // appear in the onscreen panel
        if (ImGui::Button("Refine Selected"))
        {
            std::vector<Face> faces;
            for (Face f:m.faces()) {
                if (selection[f] == 1) faces.push_back(f);
            }
            refine(icit, faces);
            selection.fill(0);
            vis_mg();
        }
        ImGui::SameLine();
        if (ImGui::Button("Uniform Refine"))
        {
            {
                Stopwatch sw(refinement_ms);
                uniform_refine(icit);
            }
            vis_mg();
        }
        ImGui::Text(("duration: "+ std::to_string(refinement_ms) + " ms").c_str());

        if (ImPlot::BeginPlot("My Plot")) {
            ImPlot::SetupAxis(ImAxis_X1, "#Triangles",ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxis(ImAxis_Y1, "residual error estimate", ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxis(ImAxis_Y2, "Optimal", ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_Opposite);
            ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);

            std::vector<double> bound1(n_tri.size()), bound2(n_tri.size());
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
    icit.requireHalfedgeVectorsInFace();
    for (int i = 0; i< 10; i++)
    {
        uniform_refine(icit);
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
    auto e = poisson_residual_error(*m, *g, f,u);

    polyscope::init();
    polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh("M", g->vertexPositions,m->getFaceVertexList());
    pm->addFaceScalarQuantity("residual error",e)->setEnabled(true);
    pm->addVertexScalarQuantity("u",u);
    pm->addVertexScalarQuantity("f",f);
    polyscope::show();
}

TEST(refineTest, laplacian_test)
{
    using namespace geometrycentral::surface;
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" / "L.stl";
    auto [m,g] = readManifoldSurfaceMesh(fds.string());
    g->requireDECOperators();
    Eigen::MatrixXd L = g->L1;
}
