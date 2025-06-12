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
void refine_callback(gcs::IntrinsicTriangulation& tri, gcs::SurfaceMesh& mesh, gcs::VertexPositionGeometry& geometry, gcs::FaceData<int>& selection) { // gets executed per-frame

    polyscope::SurfaceMesh* pm = polyscope::getSurfaceMesh("M");
    ImGuiIO& io = ImGui::GetIO();
    if (io.MouseClicked[0])
    {
        // if clicked
        glm::vec2 screenCoords{io.MousePos.x, io.MousePos.y};
        polyscope::PickResult pickResult = polyscope::pickAtScreenCoords(screenCoords);

        // check out pickResult.isHit, pickResult.structureName, pickResult.depth, etc

        // get additional information if we clicked on a mesh
        if(pickResult.isHit && pickResult.structure == pm) {
            polyscope::SurfaceMeshPickResult meshPickResult =
              pm->interpretPickResult(pickResult);

            if(meshPickResult.elementType == polyscope::MeshElement::FACE) {
                std::cout << "clicked face " << meshPickResult.index << std::endl;
                gcs::Face f = mesh.face(meshPickResult.index);
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
        std::vector<gcs::Face> faces;
        for (gcs::Face f:mesh.faces()) {
            if (selection[f] == 1) faces.push_back(f);
        }
        refine(tri, faces);
        updateTriagulationViz(tri,geometry);
        for (gcs::Face f:mesh.faces()) { selection[f] = 0; }
    }
    ImGui::SameLine();
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
    FaceData<int> selection(*m,0);

    polyscope::init();
    polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh("M", g->vertexPositions,m->getFaceVertexList());
    pm->addFaceScalarQuantity("selected faces", selection)->setEnabled(true);

    pm->setSelectionMode(polyscope::MeshSelectionMode::FacesOnly);
    // get the mouse location from ImGui




    updateTriagulationViz(icit, *g);
    polyscope::state::userCallback = [&]() { refine_callback(icit, *m, *g, selection); };// specify the callback
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
