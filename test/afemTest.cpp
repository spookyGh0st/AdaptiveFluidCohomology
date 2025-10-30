#include "Stopwatch.h"

#include <eider/AdaptiveFluidSolver.h>

#include "Evaluator.h"
#include "TaylorInitializer.h"
#include "eider/homotopy.h"
#include "eider/util.h"
#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include "geometrycentral/surface/meshio.h"
#include <algorithm>
#include <filesystem>
#include <implot.h>
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

VertexData<Vector3> intrinsic_geom(IntrinsicTriangulation& Tri, VertexPositionGeometry& inputG){
    auto& mesh = *Tri.intrinsicMesh;
    VertexData<Vector3> int_positions(mesh) ;
    for (Vertex v : mesh.vertices()) { int_positions[v] = Tri.vertexLocations[v].interpolate(inputG.vertexPositions); }
    return int_positions;
}


#include <cstddef> // for size_t
#include <utility>

void sortByNames(std::vector<std::string>& names,
                 std::vector<std::filesystem::path>& paths) {
    if (names.size() != paths.size()) {
        throw std::runtime_error("Vectors must have the same size");
    }

    // Create index vector
    std::vector<size_t> indices(names.size());
    for (size_t i = 0; i < indices.size(); ++i) indices[i] = i;

    // Sort indices based on names
    std::sort(indices.begin(), indices.end(),
              [&](size_t a, size_t b) { return names[a] < names[b]; });

    // Apply sorted order to names and paths
    std::vector<std::string> sortedNames;
    std::vector<std::filesystem::path> sortedPaths;
    sortedNames.reserve(names.size());
    sortedPaths.reserve(paths.size());

    for (size_t i : indices) {
        sortedNames.push_back(std::move(names[i]));
        sortedPaths.push_back(std::move(paths[i]));
    }

    // Replace original vectors
    names = std::move(sortedNames);
    paths = std::move(sortedPaths);
}


struct MeshSelecter {
    std::vector<std::string> names;
    std::vector<std::filesystem::path> paths;
    int selected_mesh = -1;
    std::vector<bool> toggles;
    MeshSelecter(std::string name, std::unique_ptr<ManifoldSurfaceMesh>& mesh, std::unique_ptr<VertexPositionGeometry>& geom) {
        namespace fs = std::filesystem;
        std::filesystem::path dir(__FILE__);
        dir = dir.parent_path()/ "models";
        for (const auto &entry : fs::recursive_directory_iterator(dir)) {
            if (fs::is_regular_file(entry.path()) && entry.path().extension() == ".stl") {
                names.emplace_back(entry.path().filename());
                paths.emplace_back(entry.path());
            }
        }
        sortByNames(names,paths);
        toggles = std::vector<bool>(names.size(), false);

        for (int i = 0; i < names.size(); i++) {
            if(name == names[i]) selected_mesh = i;
        }
        if(name != "" && selected_mesh == -1) throw  std::runtime_error("Mesh does not exist");
        if (selected_mesh!=-1){
            std::tie(mesh,geom) = readManifoldSurfaceMesh(paths[selected_mesh].string());
        }
    }
    bool select(std::unique_ptr<ManifoldSurfaceMesh>& mesh, std::unique_ptr<VertexPositionGeometry>& geom) {
        // Simple selection popup (if you want to show the current selection inside the Button itself,
        // you may want to build a string using the "###" operator to preserve a constant ID with a variable label)
        ImGui::PushID(this);
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
        ImGui::PopID();
        return changed;
    }
};

void visDopriConf(DOPRI5_conf &conf) {
    // Precision presets
    const char* precision_presets[] = { "Low", "Medium", "High" };
    static int current_preset = 1; // default = Medium

    if (ImGui::TreeNode("Dopri56 Configuration")) {

        if (ImGui::Combo("Precision Preset", &current_preset, precision_presets, IM_ARRAYSIZE(precision_presets))) {
            switch (current_preset) {
            case 0: // Low precision (fast)
                conf = DOPRI5Preset(DOPRI5PresetConf::LOW); break;
            case 1: // Medium precision (balanced)
                conf = DOPRI5Preset(DOPRI5PresetConf::MEDIUM); break;
            case 2: // High precision (slow, stable)
                conf = DOPRI5Preset(DOPRI5PresetConf::HIGH); break;
            }
        }

        ImGui::InputDouble("Absolute Tolerance", &conf.Atol_i, 0.0, 0.0, "%.2e");
        ImGui::InputDouble("Relative Tolerance", &conf.Rtol_i, 0.0, 0.0, "%.2e");
        ImGui::InputDouble("Step Factor Min",    &conf.facmin, 0.0, 0.0, "%.2f");
        ImGui::InputDouble("Step Factor Max",    &conf.faxmax, 0.0, 0.0, "%.2f");
        ImGui::TreePop();
    }
}

void visDoerflerConf(DoeflerConf &conf) {
    // Preset names
    const char* doerfler_presets[] = {
        "Low",
        "Medium",
        "High",
        "Uniform Refine",
        "Uniform Coarse"
    };
    static int current_preset = 0; // default = Conservative

    if (ImGui::TreeNode("Doerfler Configuration")) {

        if (ImGui::Combo("Preset", &current_preset, doerfler_presets, IM_ARRAYSIZE(doerfler_presets))) {
            switch (current_preset) {
            case 0: // Conservative refinement
                conf = DoerflerPreset(DoerflerPresetConf::LOW); break;
            case 1: // Aggressive refinement
                conf = DoerflerPreset(DoerflerPresetConf::MEDIUM); break;
            case 2: // Aggressive coarsening
                conf = DoerflerPreset(DoerflerPresetConf::HIGH); break;
            case 3: // Uniform refine
                conf = DoerflerPreset(DoerflerPresetConf::UNIFORM_REFINE); break;
            case 4: // Uniform coarse
                conf = DoerflerPreset(DoerflerPresetConf::UNIFORM_COARSE); break;
            }
        }

        ImGui::InputDouble("Refine: theta",     &conf.theta_refine);
        ImGui::InputDouble("Refine: threshold", &conf.threshold_refine, 0.0, 0.0, "%.2e");
        ImGui::InputDouble("Coarse: theta",     &conf.theta_coarse);
        ImGui::InputDouble("Coarse: threshold", &conf.threshold_coarse, 0.0,0.0, "%.2e");

        ImGui::TreePop();
    }
}


void timePlot(const std::string& name, const std::function<void()>& plotFunc) {
    ImGui::PushID(name.c_str());
    if (ImGui::TreeNode(name.c_str())) {
        if (ImPlot::BeginPlot(name.c_str())) {
            ImPlot::SetupAxis(ImAxis_X1, "Time",ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxis(ImAxis_Y1, name.c_str(), ImPlotAxisFlags_AutoFit);
            plotFunc();
            ImPlot::EndPlot();
        }
        ImGui::TreePop();
    }
    ImGui::PopID();
}



FaceData<Vector2> lamb_form(ManifoldSurfaceMesh& mesh, const VertexData<double>& w, const FaceData<Vector2>& u){
    FaceData<Vector2> lamb(mesh);
    for (Face f: mesh.faces()){ lamb[f] =Lamb(f,w,u); }
    return lamb;
}


struct AdaptiveFluidPlotter {
    std::vector<float> t_data, fluid_velocity, vorticity_data, stream_function_data, dt, lamb, dw;
    std::vector<std::vector<float>> c_data, dc;
    std::vector<std::vector<float>> hom_sum;
    std::vector<float> adapt_data;
    int adapt = 0;
    std::vector<std::vector<float>> adapt_dc;
    std::vector<float> adapt_fluid_velocity, adapt_vorticity_data, adapt_stream_function_data, adapt_lamb, adapt_dw;
    Evaluator evt, eva;
    void registerAll(Evaluator& ev,int h_size){
        ev.reg("dopri - dt", [](EvData d) { return d.dp5s.t_past; });
        ev.reg("dopri - attempts", [](EvData d) { return d.dp5s.attempts; });
        ev.reg("velocity", [](EvData d) { return L2Norm(d.vel.u,d.geom); });
        ev.reg("streamfunction", [](EvData d) { return integral(d.vel.stream_function,d.geom); });
        ev.reg("w", [](EvData d) { return integral(d.wc.w,d.geom); });
        ev.reg("dw", [](EvData d) { return integral(d.rhs.w,d.geom); });
        for (int i = 0; i < h_size; ++i) {
            ev.reg("c"+std::to_string(i), [i](EvData d) { return d.wc.c[i]; });
            ev.reg("dc"+std::to_string(i), [i](EvData d) { return d.rhs.c[i]; });
        }

    }
    AdaptiveFluidPlotter(AdaptiveFluidSolver& solver){
        registerAll(evt,solver.h.size());
        registerAll(eva,solver.h.size());
    }

    void onAdapt(AdaptiveFluidSolver& solver){
        ManifoldSurfaceMesh& mesh = solver.tri.mesh();
        IntrinsicGeometryInterface& geom = solver.tri.geom();
        wc_wrapper rhs = evalRHS(mesh,geom,solver.wc,solver.h,solver.S);
        EvData data { mesh, geom, solver.velocity(), rhs, solver.wc, solver.h, DOPRI5_sample() };
        eva.onStep(data,mesh.nFaces());
    }

    void onStep(AdaptiveFluidSolver& solver, DOPRI5_sample sample){
        ManifoldSurfaceMesh& mesh = solver.tri.mesh();
        IntrinsicGeometryInterface& geom = solver.tri.geom();
        wc_wrapper rhs = evalRHS(mesh,geom,solver.wc,solver.h,solver.S);
        EvData data { mesh, geom, solver.velocity(), rhs, solver.wc, solver.h, sample };
        evt.onStep(data,solver.elapsed_time);
    }

    void plote(std::string name, Evaluator& e){
        if(ImGui::TreeNode(name.c_str())){
            auto groups = evt.group_columns();
            for (auto const& [prefix, vecs] : groups) {
                timePlot(prefix,[&]() {
                  for(const auto i: vecs){
                      ImPlot::PlotLine(evt.y[i].name.c_str(),  e.x.data(), e.y[i].data.data(),int(e.x.size()));
                  }
                });
            }
            ImGui::TreePop();
        }

    }

    void callBack(){
        plote("Time plot", evt);
        plote("Adapt plot", eva);
    }
};

struct AdaptiveFluidPlotterComparator{
    std::vector<std::pair<std::string, AdaptiveFluidPlotter*>> plotters;
    void callBack(){
        ImGui::PushID(this);
        if(ImGui::TreeNode("Comparator")){
            timePlot("Delta Time",[&]() {
                for(const auto p: plotters){
                    if(p.second == nullptr) continue;
                    ImPlot::PlotLine((p.first + " dt").c_str(),  p.second->t_data.data(), p.second->dt.data(),int(p.second->t_data.size()));
                }
            });
            timePlot("Homology Coefficients",[&]() {
              for(const auto p: plotters){
                  if(p.second == nullptr) continue;
                  for (int j = 0; j < p.second->c_data.size() ; ++j) {
                      ImPlot::PlotLine((p.first + " c" + std::to_string(j)).c_str(), p.second->t_data.data(), p.second->c_data[j].data(),int(p.second->t_data.size())); }
              }
            });
            timePlot("velocity norm",[&]() {
              for(const auto p: plotters){
                  if(p.second == nullptr) continue;
                  ImPlot::PlotLine((p.first +" fluid").c_str(), p.second->t_data.data(), p.second->fluid_velocity.data(), int(p.second->t_data.size()));
              }
            });
            timePlot("vorticity/streamfunction norm",[&]() {
              for(const auto p: plotters){
                  if(p.second == nullptr) continue;
                  ImPlot::PlotLine((p.first + " vorticity").c_str(), p.second->t_data.data(), p.second->vorticity_data.data(), int(p.second->t_data.size()));
                  ImPlot::PlotLine((p.first + " streamfunction").c_str(), p.second->t_data.data(), p.second->stream_function_data.data(), int(p.second->t_data.size()));
              }
            });
            ImGui::TreePop();
        }
        ImGui::PopID();

    }
};

struct AdaptiveFluidVisualization {
    std::string name;
    std::unique_ptr<ManifoldSurfaceMesh> pMesh = nullptr;
    std::unique_ptr<VertexPositionGeometry> pGeom = nullptr;
    std::unique_ptr<AdaptiveFluidSolver> solver;
    std::unique_ptr<AdaptiveFluidPlotter> plotter;
    TaylorInitializer taylor;
    AdaptiveFluidSolverData solverData;
    MeshSelecter selecter;
    double delauny_angle = 25, delauny_circ = std::numeric_limits<double>::infinity();
    bool fix_c = false, delauny_refine = false;
    int start_refine = 0;

    void load() {
        auto intrT = std::make_unique<IntegerCoordinatesIntrinsicTriangulation>(*pMesh, *pGeom);
        if (delauny_refine) {
            intrT->delaunayRefine(delauny_angle,delauny_circ);
            intrT->flipToDelaunay();
            intrT->intrinsicMesh->compress();
            intrT->refreshQuantities();
        }
        auto nMesh = intrT->intrinsicMesh->copy();
        pGeom = std::make_unique<VertexPositionGeometry>(*nMesh,intrinsic_geom(*intrT,*pGeom).reinterpretTo(*nMesh));
        pMesh = std::move(nMesh);

        AdaptiveFluidSolverData data = solver? solver->data() : AdaptiveFluidSolverData();
        solver = std::make_unique<AdaptiveFluidSolver>(*pMesh,*pGeom,data);
        solver->wc.w = TaylorInitializer().wc(solver->tri.intrinsicTriangulation(),*pGeom).w;
        for (int i = 0; i < start_refine; ++i) {
            solver->adapt();
        }
        solver->wc.w = TaylorInitializer().wc(solver->tri.intrinsicTriangulation(),*pGeom).w;
        plotter = std::make_unique<AdaptiveFluidPlotter>(*solver);
    }

    AdaptiveFluidVisualization(std::string mesh = "") :selecter(mesh,pMesh,pGeom) {
        if(pMesh != nullptr) load();
    }

    bool init = true;

    struct SolverState {
        bool running = false;
        bool adaptive_run = false;
        int adapt_after = 1;
        int step= 0;
    } state;

    void visualize() {
        if (!solver) return;

        ManifoldSurfaceMesh& mesh = solver->tri.mesh();
        IntrinsicTriangulation& intrT = solver->tri.intrinsicTriangulation();
        mesh.compress();
        VertexData<Vector3> int_positions = intrinsic_geom(intrT,*pGeom);
        VertexPositionGeometry g(mesh,int_positions);
        polyscope::SurfaceMesh* polym = polyscope::registerSurfaceMesh(name, int_positions,mesh.getFaceVertexList(), polyscopePermutations(mesh));
        g.requireFaceTangentBasis();
        FaceData<Vector3> e1(mesh),e2(mesh);
        for (Face f: mesh.faces()) { e1[f] = g.faceTangentBasis[f][0], e2[f] = g.faceTangentBasis[f][1];  }

        std::vector<FaceData<double>> dw_face(solver->h.size(), FaceData<double>(mesh));
        auto v = solver->velocity();
        auto dwdc = evalRHS(mesh,solver->tri.geom(),solver->wc,solver->h,solver->S,&dw_face);
        polym->addFaceTangentVectorQuantity("velocity",v.u,e1,e2);
        auto* vorticityV = polym->addVertexScalarQuantity("vorticity",solver->wc.w);
        if(init) vorticityV->setEnabled(true);
        vorticityV->setColorMap("coolwarm");
        polym->addVertexScalarQuantity("vorticity - change",dwdc.w);
        polym->addVertexScalarQuantity("stream_function",v.stream_function);


        auto full_h =solver->hom.fullHarmonicBasis();
        for (int i = 0; i < solver->h.size(); ++i) {
            auto &basis = solver->hom.homologyB[i];
            auto df = delta_form(mesh, basis);
            polym->addEdgeScalarQuantity("harmonic form Whitney" + std::to_string(i),df);

            FaceData<Vector2> b = solver->h[i];
            polym->addFaceTangentVectorQuantity("harmonic form " + std::to_string(i),b,e1,e2);
            auto dwf = dw_face[i];
            polym->addFaceScalarQuantity("dw_face - " + std::to_string(i),dwf);

            polym->addEdgeScalarQuantity("h -df" + std::to_string(i),full_h.df[i]);
            polym->addEdgeScalarQuantity("h -proj df" + std::to_string(i),full_h.proj_df[i]);
            polym->addFaceTangentVectorQuantity("h -unorthorgonal" + std::to_string(i),full_h.h_unorth[i],e1,e2);
            polym->addFaceTangentVectorQuantity("h -orthorgonal" + std::to_string(i),full_h.h_orth[i],e1,e2);
        }
        for (int i = 0; i < solver->h_interpolated.size(); ++i) {
            polym->addFaceTangentVectorQuantity("h -interpolated" + std::to_string(i),solver->h_interpolated[i],e1,e2);
        }

        polym->addFaceTangentVectorQuantity("Lamb Form ", lamb_form(mesh,solver->wc.w,v.u),e1,e2);

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

    void callBackImpl() {
        taylor.callback();
        if (ImGui::TreeNode("Delauny Settings")){
            ImGui::Checkbox("Delauny Refine", &delauny_refine);
            ImGui::InputDouble("Delauny Angle", &delauny_angle);
            ImGui::InputDouble("Delauny circ", &delauny_circ);
            ImGui::InputInt("Start Refinement", &start_refine);
            ImGui::TreePop();
        }
        if (selecter.select(pMesh, pGeom) || ImGui::Button("reset")) { load(); visualize(); }
        if (!solver) { return; }
        if (ImGui::TreeNode("Control")) {

            ImGui::SameLine();
            if (ImGui::Button("adapt")) {
                solver->adapt();
                plotter->onAdapt(*solver);
                visualize();
            }

            ImGui::Checkbox("Run", &state.running);
            ImGui::Checkbox("Adapt Space", &solver->adapte_space);
            ImGui::Checkbox("Adapt Time", &solver->adapt_time);
            ImGui::Checkbox("Fix C", &fix_c);
            ImGui::Checkbox("Use interpolated h", &solver->use_interpolated_h);
            ImGui::InputDouble("Delta Time", &solver->dt);
            if (ImGui::Button("Set h as interpolated h"))
                solver->h_interpolated = solver->h;

            if (state.running || ImGui::Button("Advance")) {
                wc_wrapper wc = solver->wc;
                DOPRI5_sample dopri5_sample = solver->step();
                state.step++;
                plotter->onStep(*solver,dopri5_sample);
                if (fix_c)
                    solver->wc.c = wc.c;
                visualize();
            }
            ImGui::TreePop();
        }

        visDoerflerConf(solver->doerflerConf);
        visDopriConf(solver->conf);

        plotter->callBack();
    }

    void callBack() {
        if (init) {visualize(); init = false;}
        ImGui::PushID(this);
        if (ImGui::TreeNode("Mesh")){
            callBackImpl();
            ImGui::TreePop();
        }
        ImGui::PopID();
    }
};

TEST(afemTest, AdaptiveFluidCohomology)
{
    AdaptiveFluidVisualization visA("cheese_min.stl");
    visA.name = "Normal";
    AdaptiveFluidVisualization visB;
    visB.name = "Refined";
    AdaptiveFluidPlotterComparator cmp;
    polyscope::init();
    ImPlot::CreateContext();

    polyscope::state::userCallback = [&]()
    {
      visA.callBack();
      visB.callBack();
      cmp.plotters = { {"Normal", visA.plotter.get() }, { "Refined", visB.plotter.get() }};
      cmp.callBack();
    };

    polyscope::show();
}

TEST(afemTest, testPathConsistency){
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" / "disk.stl";
    auto [parent_m,parent_g] = readManifoldSurfaceMesh(fds.string());

    AdaptiveTriangulation atri(*parent_m,*parent_g);
    auto& icit = atri.intrinsicTriangulation();
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
