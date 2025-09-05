#include "Stopwatch.h"

#include <eider/AdaptiveFluidSolver.h>

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

#include <cstddef> // for size_t

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
        sortByNames(names,paths);
        toggles = std::vector<bool>(names.size(), false);
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
                conf.Rtol_i  = 1e-3;
                conf.Atol_i  = 1e-6;
                conf.facmin  = 0.2;
                conf.faxmax  = 10.0;
                break;

            case 1: // Medium precision (balanced)
                conf.Rtol_i  = 1e-6;
                conf.Atol_i  = 1e-9;
                conf.facmin  = 0.2;
                conf.faxmax  = 5.0;
                break;

            case 2: // High precision (slow, stable)
                conf.Rtol_i  = 1e-8;
                conf.Atol_i  = 1e-11;
                conf.facmin  = 0.25;
                conf.faxmax  = 2.0;
                break;
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
        "Conservative",
        "Aggressive Refine",
        "Aggressive Coarse",
        "Uniform Refine",
        "Uniform Coarse"
    };
    static int current_preset = 0; // default = Conservative

    if (ImGui::TreeNode("Doerfler Configuration")) {

        if (ImGui::Combo("Preset", &current_preset, doerfler_presets, IM_ARRAYSIZE(doerfler_presets))) {
            switch (current_preset) {
            case 0: // Conservative refinement
                conf.theta_refine     = 0.1;
                conf.threshold_refine = 1e-6;
                conf.theta_coarse     = 0.1;
                conf.threshold_coarse = 1e-8;
                break;

            case 1: // Aggressive refinement
                conf.theta_refine     = 0.5;
                conf.threshold_refine = 1e-8;
                conf.theta_coarse     = 0.05;
                conf.threshold_coarse = 1e-9;
                break;

            case 2: // Aggressive coarsening
                conf.theta_refine     = 0.3;
                conf.threshold_refine = 1e-6;
                conf.theta_coarse     = 0.3;
                conf.threshold_coarse = 1e-6;
                break;

            case 3: // Uniform refine
                conf.theta_refine     = 1.0;   // force refine
                conf.threshold_refine = 0.0;   // no threshold
                conf.theta_coarse     = 0.0;   // disable coarsening
                conf.threshold_coarse = 0.0;   // disable coarsening
                break;

            case 4: // Uniform coarse
                conf.theta_refine     = 0.0;   // disable refining
                conf.threshold_refine = 1e9;   // effectively prevent refining
                conf.theta_coarse     = 1.0;   // force coarsening
                conf.threshold_coarse = 1e9;   // no threshold
                break;
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

double L2Norm(FaceData<Vector2> d, IntrinsicGeometryInterface& geom){
    double s = 0;
    for (Face f: geom.mesh.faces()) { s += d[f].norm2() * geom.faceAreas[f]; }
    return std::sqrt(s);
}
double L2Norm(VertexData<double> d, IntrinsicGeometryInterface& geom){
    double s = 0;
    for (Face f: geom.mesh.faces()){
        double fs = 0;
        for(Vertex v: f.adjacentVertices()){
            // TODO: abs only debug
            fs += std::abs(d[v]);
        }
        s += fs/3*geom.faceAreas[f];
    }
    return s;
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
    AdaptiveFluidPlotter(AdaptiveFluidSolver& solver){
        c_data.resize(solver.h.size());
        hom_sum.resize(solver.h.size());
        dc.resize(solver.h.size());
        adapt_dc.resize(solver.h.size());
    }

    void onAdapt(AdaptiveFluidSolver& solver){
        adapt_data.push_back(adapt++);

        ManifoldSurfaceMesh& mesh = solver.tri.mesh();
        IntrinsicGeometryInterface& geom = solver.tri.geom();

        auto vel = solver.velocity();
        adapt_fluid_velocity.push_back(L2Norm(vel.u,geom));
        adapt_vorticity_data.push_back(L2Norm(solver.wc.w,geom));
        adapt_stream_function_data.push_back(L2Norm(vel.stream_function,geom));
        adapt_lamb.push_back(L2Norm(lamb_form(mesh,solver.wc.w,vel.u),geom));

        wc_wrapper rhs = evalRHS(mesh,geom,solver.wc,solver.h,solver.S);
        for (int i = 0; i < solver.h.size(); ++i) {
            adapt_dc[i].push_back(rhs.c[i]);
        }
        adapt_dw.push_back(L2Norm(rhs.w,geom));
    }

    void onStep(AdaptiveFluidSolver& solver){
        t_data.push_back(solver.elapsed_time);
        dt.push_back(solver.dt);

        ManifoldSurfaceMesh& mesh = solver.tri.mesh();
        IntrinsicGeometryInterface& geom = solver.tri.geom();
        wc_wrapper rhs = evalRHS(mesh,geom,solver.wc,solver.h,solver.S);

        for (int i = 0; i < solver.h.size(); ++i) {
            c_data[i].push_back(solver.wc.c[i]);
            hom_sum[i].push_back(L2Norm(solver.h[i],geom));
            dc[i].push_back(rhs.c[i]);
        }

        auto vel = solver.velocity();

        fluid_velocity.push_back(L2Norm(vel.u,geom));
        vorticity_data.push_back(L2Norm(solver.wc.w,geom));
        vorticity_data.push_back(L2Norm(vel.stream_function,geom));
        lamb.push_back(L2Norm(lamb_form(mesh,solver.wc.w,vel.u),geom));
        dw.push_back(L2Norm(rhs.w,geom));
    }

    void callBack(){
        ImGui::PushID(this);
        if(ImGui::TreeNode("Time Plots")){
            timePlot("Delta Time",[&]() {
              ImPlot::PlotLine("dt",  t_data.data(), dt.data(),int(t_data.size()));
            });
            timePlot("Homology Coefficients",[&]() {
              for (int j = 0; j < c_data.size() ; ++j) { ImPlot::PlotLine(("c" + std::to_string(j)).c_str(), t_data.data(), c_data[j].data(),int(t_data.size())); }
            });
            timePlot("Harmonic basis norm",[&]() {
              for (int j = 0; j < hom_sum.size(); ++j) { ImPlot::PlotLine(("c" + std::to_string(j)).c_str(), t_data.data(), hom_sum[j].data(), int(t_data.size())); }
            });
            timePlot("velocity norm",[&]() {
              ImPlot::PlotLine("fluid", t_data.data(), fluid_velocity.data(), int(t_data.size()));
            });
            timePlot("vorticity/streamfunction norm",[&]() {
              ImPlot::PlotLine("vorticity", t_data.data(), vorticity_data.data(), int(t_data.size()));
              ImPlot::PlotLine("streamfunction", t_data.data(), stream_function_data.data(), int(t_data.size()));
            });
            timePlot("dc",[&]() {
              for (int j = 0; j < c_data.size() ; ++j) { ImPlot::PlotLine(("dc" + std::to_string(j)).c_str(), t_data.data(), dc[j].data(),int(t_data.size())); }
            });
            timePlot("dw",[&]() {
              ImPlot::PlotLine("dw", t_data.data(), dw.data(),int(t_data.size()));
            });
            ImGui::TreePop();
        }
        if(ImGui::TreeNode("Adapt_plots")){
            timePlot("velocity norm",[&]() {
              ImPlot::PlotLine("velocity", adapt_data.data(), adapt_fluid_velocity.data(), int(adapt_data.size()));
            });
            timePlot("vorticity/streamfunction norm",[&]() {
              ImPlot::PlotLine("vorticity", adapt_data.data(), adapt_vorticity_data.data(), int(adapt_data.size()));
              ImPlot::PlotLine("streamfunction", adapt_data.data(), adapt_stream_function_data.data(), int(adapt_data.size()));
            });
            timePlot("Lamb norm",[&]() {
              ImPlot::PlotLine("Lamb form", adapt_data.data(), adapt_lamb.data(), int(adapt_data.size()));
            });
            timePlot("dc",[&]() {
              for (int j = 0; j < c_data.size() ; ++j) { ImPlot::PlotLine(("dc" + std::to_string(j)).c_str(), adapt_data.data(), adapt_dc[j].data(),int(adapt_data.size())); }
            });
            timePlot("dw",[&]() {
              ImPlot::PlotLine("dw", adapt_data.data(), adapt_dw.data(),int(adapt_data.size()));
            });
            ImGui::TreePop();
        }
        ImGui::PopID();
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
    MeshSelecter selecter;
    TaylorInitializer taylor;
    std::unique_ptr<ManifoldSurfaceMesh> pMesh;
    std::unique_ptr<VertexPositionGeometry> pGeom;
    std::unique_ptr<IntrinsicTriangulation> intrT;
    std::unique_ptr<AdaptiveTriangulation> aTri;
    std::unique_ptr<AdaptiveFluidSolver> solver;
    std::unique_ptr<AdaptiveFluidPlotter> plotter;
    DOPRI5_conf dopri5;
    DoeflerConf doefler;
    double delauny_angle = 25, delauny_circ = std::numeric_limits<double>::infinity();
    bool normalize_basis = true, fix_c = false;

    void load() {
        if(solver){ dopri5= solver->conf; doefler = solver->doerflerConf; }
        intrT = std::make_unique<IntegerCoordinatesIntrinsicTriangulation>(*pMesh, *pGeom);
        intrT->delaunayRefine(delauny_angle,delauny_circ);
        intrT->intrinsicMesh->compress();
        intrT->refreshQuantities();
        auto nMesh = intrT->intrinsicMesh->copy();
        pGeom = std::make_unique<VertexPositionGeometry>(*nMesh,intrinsic_geom(*intrT,*pGeom).reinterpretTo(*nMesh));
        pMesh = std::move(nMesh);
        intrT = std::make_unique<IntegerCoordinatesIntrinsicTriangulation>(*pMesh, *pGeom);

        intrT->mesh.compress(); intrT->refreshQuantities();
        aTri = std::make_unique<AdaptiveTriangulation>(*intrT);
        wc_wrapper wc = taylor.wc(*intrT,*pGeom);
        solver = std::make_unique<AdaptiveFluidSolver>(*aTri,wc,dopri5,doefler);
        plotter = std::make_unique<AdaptiveFluidPlotter>(*solver);
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
        polyscope::SurfaceMesh* polym = polyscope::registerSurfaceMesh(name, int_positions,mesh.getFaceVertexList(), polyscopePermutations(mesh));
        // polym->setAllPermutations(polyscopePermutations(mesh));

        g.requireFaceTangentBasis();
        FaceData<Vector3> e1(mesh),e2(mesh);
        for (Face f: mesh.faces()) { e1[f] = g.faceTangentBasis[f][0], e2[f] = g.faceTangentBasis[f][1];  }

        std::vector<FaceData<double>> dw_face(solver->h.size(), FaceData<double>(mesh));
        auto v = solver->velocity();
        auto dwdc = evalRHS(mesh,solver->tri.geom(),solver->wc,solver->h,solver->S,&dw_face);
        polym->addFaceTangentVectorQuantity("velocity",v.u,e1,e2);
        polym->addVertexScalarQuantity("vorticity",solver->wc.w);
        polym->addVertexScalarQuantity("vorticity - change",dwdc.w);
        polym->addVertexScalarQuantity("stream_function",v.stream_function);


        AdaptivePressureProjectionSolver pp_solver;
        pp_solver.compute(g);
        for (int i = 0; i < solver->h.size(); ++i) {
            auto &basis = solver->hom.homologyB[i];
            auto df = delta_form(mesh, basis);
            EdgeData<double> pf = pp_solver.solveWithGuess(mesh, df, &solver->hom.pf_guess[i]);
            FaceData<Vector2> h = whitney_interpolation(mesh, g, pf);
            polym->addFaceTangentVectorQuantity("harmonic form Whitney" + std::to_string(i),h,e1,e2);

            FaceData<Vector2> b = solver->h[i];
            if (normalize_basis) for (Face f: mesh.faces()) { b[f] /= g.faceArea(f); }
            polym->addFaceTangentVectorQuantity("harmonic form " + std::to_string(i),b,e1,e2);
            auto dwf = normalize_basis ? dw_face[i] / g.faceAreas : dw_face[i];
            polym->addFaceScalarQuantity("dw_face - " + std::to_string(i),dwf);
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
            ImGui::InputDouble("Delauny Angle", &delauny_angle);
            ImGui::InputDouble("Delauny circ", &delauny_circ);
            ImGui::TreePop();
        }
        if (selecter.select(pMesh, pGeom) || ImGui::Button("reset")) { load(); visualize(); }
        if (!aTri) { return; }
        if (ImGui::TreeNode("Control")) {

            ImGui::SameLine();
            if (ImGui::Button("adapt")) {
                solver->adapt();
                plotter->onAdapt(*solver);
                visualize();
            }

            ImGui::Checkbox("Run", &state.running);
            ImGui::SameLine();
            ImGui::Checkbox("Adaptive", &state.adaptive_run);

            ImGui::InputInt("Adapt after:", &state.adapt_after);
            if (state.running || ImGui::Button("Advance")) {
                wc_wrapper wc = solver->wc;
                solver->step();
                state.step++;
                if (state.adaptive_run && state.step % state.adapt_after == 0){
                    solver->adapt();
                    plotter->onAdapt(*solver);
                }
                plotter->onStep(*solver);
                if (fix_c)
                    solver->wc.c = wc.c;
                visualize();
            }
            ImGui::TreePop();
        }

        if(ImGui::TreeNode("Debug Configuration")){
            if(ImGui::Checkbox("Normalize Basis", &normalize_basis)) visualize();
            ImGui::Checkbox("Fix C", &fix_c);
            ImGui::TreePop();
        }
        visDoerflerConf(solver->doerflerConf);
        visDopriConf(solver->conf);

        plotter->callBack();
    }

    void callBack() {
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
    AdaptiveFluidVisualization visA;
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

