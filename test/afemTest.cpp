#include "Stopwatch.h"

#include <eider/AdaptiveFluidSolver.h>

#include "eider/homotopy.h"
#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include "geometrycentral/surface/meshio.h"
#include <filesystem>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include <gtest/gtest.h>

using namespace geometrycentral;
using namespace geometrycentral::surface;

TEST(afemTest, testDefaultBehauvior)
{
  std::filesystem::path fds(__FILE__);
  fds = fds.parent_path()/ "models" /"L.stl";
  auto [pm,pg] = readManifoldSurfaceMesh(fds.string());
  auto intrTri = IntegerCoordinatesIntrinsicTriangulation(*pm,*pg);
  ManifoldSurfaceMesh& m = *intrTri.intrinsicMesh;
  VertexData<double> f(m,0);
  Face face = m.face(0);
  f[*face.adjacentVertices().begin()] = 1;
  intrTri.faceInsertionCallbackList.push_back([&](Face face, Vertex v) {f[v] = 2; });
  Vertex v = intrTri.splitFace(face,Vector3(0.33,0.333,0.33));
  ASSERT_EQ(f[v], 2);
}

wc_wrapper init_taylor(SurfaceMesh& mesh, VertexData<Vector3> geo, double vorticity_distance, Vector2 offset, Vector2 cuttof, Vector2 cuttof_offset) {
    VertexData<double> w(mesh,0);
    double k = 2*PI / vorticity_distance;
    double A = 1;
    for (Vertex v : mesh.vertices())
    {
        Vector3 p = geo[v];
        double x = p.x  + offset.x, y = p.y + offset.y;
        if (abs(p.x+cuttof_offset.x) > cuttof.x || abs(p.y+cuttof_offset.y) > cuttof.y) { continue;}
        w[v] = 2 * A *k * cos(k*x) * cos(k*y);
    }
    wc_wrapper wc;
    wc.w = w;
    return wc;
}
VertexData<Vector3> intrinsic_geom(IntrinsicTriangulation& Tri, VertexPositionGeometry& inputG){
    auto& mesh = *Tri.intrinsicMesh;
    VertexData<Vector3> int_positions(mesh) ;
    for (Vertex v : mesh.vertices()) { int_positions[v] = Tri.vertexLocations[v].interpolate(inputG.vertexPositions); }
    return int_positions;
}

struct TaylorInitializer {
    float v_dist  = 0.5, x_add = v_dist/4, y_add = 0, x_cuttof = v_dist/2, y_cuttof = v_dist/4, x_cutt_offset = 0, y_cutt_offset = 0.5;
    wc_wrapper wc(IntrinsicTriangulation& intTri, VertexPositionGeometry& pg) {
        return init_taylor(*intTri.intrinsicMesh, intrinsic_geom(intTri,pg), v_dist, Vector2(x_add, y_add), Vector2(x_cuttof, y_cuttof), Vector2(x_cutt_offset, y_cutt_offset));
    }
    void callback() {
        if (ImGui::TreeNode("Initial Taylor Vorticises")) {
            ImGui::InputFloat("vorticity distance",&v_dist,0.125,0.5);
            ImGui::InputFloat("x_add",&x_add,0.125,0.5); ImGui::InputFloat("y_add",&y_add,0.125,0.5);
            ImGui::InputFloat("x_cuttof",&x_cuttof,0.125,0.5); ImGui::InputFloat("y_cuttoff",&y_cuttof,0.125,0.5);
            ImGui::InputFloat("x_cut offset",&x_cutt_offset,0.125,0.5); ImGui::InputFloat("y_cut offset",&y_cutt_offset,0.125,0.5);
            ImGui::TreePop();
        }
    }
};

struct MeshSelecter {
    std::vector<std::string> names;
    std::vector<std::filesystem::path> paths;
    int selected_mesh = -1;
    std::vector<bool> toggles;
    MeshSelecter() {
        namespace fs = std::filesystem;
        std::filesystem::path dir(__FILE__);
        dir = dir.parent_path()/ "models";
        for (const auto &entry : fs::recursive_directory_iterator(dir)) {
            if (fs::is_regular_file(entry.path()) && entry.path().extension() == ".stl") {
                names.emplace_back(entry.path().filename());
                paths.emplace_back(entry.path());
            }
        }
        toggles = std::vector<bool>(names.size(), false);
    }
    bool select(std::unique_ptr<ManifoldSurfaceMesh>& mesh, std::unique_ptr<VertexPositionGeometry>& geom) {
        // Simple selection popup (if you want to show the current selection inside the Button itself,
        // you may want to build a string using the "###" operator to preserve a constant ID with a variable label)
        bool changed = false;
        if (ImGui::Button("Select Mesh"))
            ImGui::OpenPopup("my_select_popup");
        ImGui::SameLine();
        ImGui::TextUnformatted(selected_mesh == -1 ? "<None>" : names[selected_mesh].c_str());
        if (ImGui::BeginPopup("my_select_popup"))
        {
            ImGui::SeparatorText("Select Mesh");
            for (int i = 0; i < names.size(); i++) {
                if (ImGui::Selectable(names[i].c_str())) {
                    changed = true;
                    selected_mesh = i;
                    std::tie(mesh,geom) = readManifoldSurfaceMesh(paths[selected_mesh].string());
                    std::cout << "Selected " << names[selected_mesh] << std::endl;
                }
            }
            ImGui::EndPopup();
        }
        return changed;
    }
};

void visDopriConf(DOPRI5_conf &conf) {
    if (ImGui::TreeNode("Dopri56 Configuration")) {
        ImGui::InputDouble("Absolut Error", &conf.Atol_i, 0, 0, "%.10f");
        ImGui::InputDouble("relative Error", &conf.Rtol_i, 0, 0, "%.10f");
        ImGui::InputDouble("fac Min", &conf.facmin);
        ImGui::InputDouble("fac max", &conf.faxmax);
        ImGui::TreePop();
    }
}
void visDoerflerConf(DoeflerConf &conf) {
    if (ImGui::TreeNode("Doerfler Configuration")) {
        ImGui::InputDouble("Coarse: theta",&conf.theta_coarse);
        ImGui::InputDouble("Coarse: threshold",&conf.threshold_coarse);
        ImGui::InputDouble("Refine: theta",&conf.theta_refine);
        ImGui::InputDouble("Refine: threshold",&conf.threshold_refine);
        ImGui::TreePop();
    }
}

struct AdaptiveFluidVisualization {
    MeshSelecter selecter;
    TaylorInitializer taylor;
    std::unique_ptr<ManifoldSurfaceMesh> pMesh;
    std::unique_ptr<VertexPositionGeometry> pGeom;
    std::unique_ptr<IntrinsicTriangulation> intrT;
    std::unique_ptr<AdaptiveTriangulation> aTri;
    std::unique_ptr<AdaptiveFluidSolver> solver;
    DOPRI5_conf dopri5;
    DoeflerConf doefler;
    void load() {
        intrT = std::make_unique<IntegerCoordinatesIntrinsicTriangulation>(*pMesh, *pGeom);
        intrT->flipToDelaunay();
        aTri = std::make_unique<AdaptiveTriangulation>(*intrT);
        wc_wrapper wc = taylor.wc(*intrT,*pGeom);
        solver = std::make_unique<AdaptiveFluidSolver>(*aTri,wc,dopri5,doefler);
    }

    struct SolverState {
        bool running = false;
        bool adaptive_run = false;
        int adapt_after = 1;
        int step= 0;
    } state;

    void visualize() {
        if (!aTri) return;

        ManifoldSurfaceMesh& mesh = aTri->mesh();
        mesh.compress();
        VertexData<Vector3> int_positions = intrinsic_geom(*intrT,*pGeom);
        VertexPositionGeometry g(mesh,int_positions);
        polyscope::SurfaceMesh* polym = polyscope::registerSurfaceMesh("M", int_positions,mesh.getFaceVertexList(), polyscopePermutations(mesh));
        // polym->setAllPermutations(polyscopePermutations(mesh));

        g.requireFaceTangentBasis();
        FaceData<Vector3> e1(mesh),e2(mesh);
        for (Face f: mesh.faces()) { e1[f] = g.faceTangentBasis[f][0], e2[f] = g.faceTangentBasis[f][1];  }

        auto v = solver->velocity();
        polym->addFaceTangentVectorQuantity("velocity",v.u,e1,e2);
        polym->addVertexScalarQuantity("vorticity",solver->wc.w);

        for (int i = 0; i < solver->h.size(); ++i) {
            polym->addFaceTangentVectorQuantity("harmonic form " + std::to_string(i),solver->h[i],e1,e2);
        }

        HalfedgeData<int> d(mesh,0);
        for (int i = 0; i < solver->hom.homologyB.size(); ++i) {
            for (Halfedge e: mesh.halfedges()){
                if(solver->hom.homologyB[i].nextLeft[e].has_value()) {
                    if (*solver->hom.homologyB[i].nextLeft[e]) d[e] = i*2+1; else d[e] = i*2+2;
                }
            }
        }
        polym->addHalfedgeScalarQuantity("Hom", d);
    }

    void callBack() {
        if (selecter.select(pMesh,pGeom) ||ImGui::Button("reset")) { load(); visualize(); }
        taylor.callback();

        if (!aTri) return;

        if (ImGui::Button("adapt")) { solver->adapt(); visualize(); }
        ImGui::Checkbox("Run", &state.running);
        ImGui::Checkbox("Refine while running", &state.adaptive_run);
        ImGui::InputInt("Adapt after:", &state.adapt_after);
        if (state.running || ImGui::Button("Advance")) {
            solver->step(); state.step++;
            if (state.adaptive_run && state.step % state.adapt_after == 0) solver->adapt();
            visualize();
        }
        ImGui::Text("dt: %.6f", solver->dt);
        for (auto& c: solver->wc.c) {
            ImGui::Text("c: %.6f", c);
        }

        visDoerflerConf(solver->doerflerConf);
        visDopriConf(solver->conf);
    }
};

TEST(afemTest, AdaptiveFluidCohomology)
{
    AdaptiveFluidVisualization vis;
    polyscope::init();

    polyscope::state::userCallback = [&]()
    {
      vis.callBack();
    };

    polyscope::show();
}

TEST(afemTest, testPathConsistency){
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" / "disk.stl";
    auto [parent_m,parent_g] = readManifoldSurfaceMesh(fds.string());

    IntegerCoordinatesIntrinsicTriangulation icit(*parent_m,*parent_g);
    AdaptiveTriangulation atri(icit);
    ManifoldSurfaceMesh& m = atri.mesh();
    IntrinsicGeometryInterface& g = atri.geom();

    auto homotopy_b= greedy_homotopy_basis(m,g, arbitrary_base_face(m));
    auto homologyBasis = singular_homology_basis(m,homotopy_b);

    auto& h = homologyBasis;
    for (int h_idx = 0; h_idx < h.size(); ++h_idx) {
        icit.edgeSplitCallbackList.emplace_back([&,h_idx](Edge e, Halfedge he1, Halfedge he2) {
          onSplit(e, he1, he2, h[h_idx].nextLeft); });
    }
    for (int h_idx = 0; h_idx < h.size(); ++h_idx) {
        icit.edgeCollapseCallbackList.push_back([&,h_idx](Halfedge he) {
          onCollapse(he, h[h_idx].nextLeft); });
    }

    bool first = true;
    auto vis = [& ](){
      m.compress();
      VertexData<Vector3> int_positions(m) ;
      for (Vertex v : m.vertices()) {
          int_positions[v] = icit.vertexLocations[v].interpolate(parent_g->vertexPositions);
      }



      polyscope::SurfaceMesh* pm_A = polyscope::registerSurfaceMesh("original", parent_g->vertexPositions,parent_m->getFaceVertexList());
      polyscope::SurfaceMesh* pm_B = polyscope::registerSurfaceMesh("Intrinsic", int_positions,m.getFaceVertexList());
      pm_B->setAllPermutations(polyscopePermutations(m));
      pm_A->addFaceScalarQuantity("index",parent_m->getFaceIndices());
      pm_B->addFaceScalarQuantity("index",m.getFaceIndices());
      pm_B->addEdgeScalarQuantity("edgeCoords",icit.normalCoordinates.edgeCoords);
      pm_B->addHalfedgeScalarQuantity("roundabouts",icit.normalCoordinates.roundabouts);
      pm_B->addVertexScalarQuantity("roundabouts degree",icit.normalCoordinates.roundaboutDegrees);
      pm_B->addEdgeScalarQuantity("lenghts",icit.edgeLengths);
      pm_A->setEdgeWidth(1);
      pm_B->setEdgeWidth(1);
      if(first){
          pm_A->translate(glm::vec3(-1,0,0));
          pm_B->translate(glm::vec3( 1,0,0)); first = false;
      }
      HalfedgeData<int> d(m,0);
      for (int i = 0; i < homologyBasis.size(); ++i) {
          for (Halfedge e: m.halfedges()){
              if(homologyBasis[i].nextLeft[e].has_value()) {
                  if (*homologyBasis[i].nextLeft[e]) d[e] = i*2+1; else d[e] = i*2+2;
                  if (!e.isInterior()) d[e.twin()] = i*2+3;
              }
          }
      }
      pm_B->addHalfedgeScalarQuantity("Hom", d);
    };


    polyscope::init();
    std::vector<Face> faces {  };
    for (Face f: m.faces()) faces.push_back(f);
    vis();

    polyscope::state::userCallback = [&]()
    {
      if (ImGui::Button("Refine"))
      {
          std::vector<Face> faces {  };
          for (Face f: m.faces()) faces.push_back(f);
          atri.refine(faces);
          vis();
      }
      if (ImGui::Button("Coarse")) {
          std::vector<Face> faces {  };
          for (Face f: m.faces()) faces.push_back(f);
          atri.coarse(faces);
          vis();
      }
        ImGui::ShowDemoWindow();
    };

    polyscope::show();
}

