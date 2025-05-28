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
void refine_callback(gcs::IntrinsicTriangulation& tri, gcs::VertexPositionGeometry& geometry) { // gets executed per-frame

    // Build a UI element to edit a parameter, which will
    // appear in the onscreen panel
    if (ImGui::Button("Uniform Refine"))
    {
        {
            Stopwatch sw(refinement_ms);
            uniform_refine(tri);
        }
        updateTriagulationViz(tri,geometry);
    }
    ImGui::Text(("duration: "+ std::to_string(refinement_ms) + " ms").c_str());
}

TEST(refineTest, testSplit)
{
    using namespace geometrycentral::surface;
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" / "torus.stl";
    auto [m,g] = readManifoldSurfaceMesh(fds.string());

    IntegerCoordinatesIntrinsicTriangulation icit(*m,*g);
    icit.flipToDelaunay();
    // uniform_refine(icit);
    // uniform_refine(icit);
    // uniform_refine(icit);

    polyscope::init();
    polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh("M", g->vertexPositions,m->getFaceVertexList());
    updateTriagulationViz(icit, *g);
    polyscope::state::userCallback = [&]() { refine_callback(icit, *g); };// specify the callback
    polyscope::show();
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

TEST(refineTest, laplacian_test)
{
    using namespace geometrycentral::surface;
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" / "L.stl";
    auto [m,g] = readManifoldSurfaceMesh(fds.string());
    g->requireDECOperators();
    Eigen::MatrixXd L = g->L1;
}
