#include "Evaluator.h"
#include "TaylorInitializer.h"
#include "eider/util.h"
#include "eider/Stopwatch.h"
#include <algorithm>
#include <eider/AdaptiveFluidSolver.h>
#include <filesystem>
#include <gtest/gtest.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

namespace fs = std::filesystem;

fs::path results_folder() {
    fs::path file_path = __FILE__;
    fs::path project_root = file_path.parent_path().parent_path().parent_path();
    fs::path results_f = project_root/ "tex" / "thesis"/"figures"/"results";
    return results_f;
}

fs::path run_folder(std::filesystem::path results_f) {
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&t), "%Y-%m-%dT%H-%M-%S");
    fs::path test_case = results_f / ("run_" + ss.str());
    fs::create_directory(test_case);
    return test_case;
}

void copyFolder(const fs::path& source, const fs::path& destination) {
    try {
        if (!fs::exists(source) || !fs::is_directory(source)) {
            throw std::runtime_error("Source folder does not exist or is not a directory.");
        }

        // Create the destination folder if it does not exist
        if(fs::exists(destination)){
            if (!destination.filename().string().starts_with("run_tc")) { throw std::runtime_error("Careful, we delete it here."); }

            fs::remove_all(destination);
        }
        fs::create_directories(destination);

        // Iterate through all files and directories in source
        for (const auto& entry : fs::recursive_directory_iterator(source)) {
            const auto& path = entry.path();
            auto relativePath = fs::relative(path, source);
            fs::path destPath = destination / relativePath;

            if (fs::is_directory(path)) {
                fs::create_directories(destPath);
            } else if (fs::is_regular_file(path)) {
                fs::copy_file(path, destPath, fs::copy_options::overwrite_existing);
            } else {
                std::cout << "Skipping non-regular file: " << path << "\n";
            }
        }

        std::cout << "Folder copied successfully.\n";
    }
    catch (const fs::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << '\n';
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
    }
}

namespace geometrycentral::surface{

VertexData<Vector3> intrinsic_geom(IntrinsicTriangulation& Tri, VertexPositionGeometry& inputG){
    auto& mesh = *Tri.intrinsicMesh;
    VertexData<Vector3> int_positions(mesh) ;
    for (Vertex v : mesh.vertices()) { int_positions[v] = Tri.vertexLocations[v].interpolate(inputG.vertexPositions); }
    return int_positions;
}


class TestCase {
  public:
    virtual ~TestCase() = default;
    std::chrono::duration<double> total_time;
    std::string name;
    std::string label;
    TestCase(std::string name, std::string label) :name(name), label(label), total_time(0) {}
    virtual void runUntil(double t) = 0;
    virtual void write(fs::path) = 0;
    virtual polyscope::SurfaceMesh* visualize(fs::path f_screenshots = fs::path()) = 0;
};

using MeshP = std::unique_ptr<ManifoldSurfaceMesh>;
using GeomP = std::unique_ptr<VertexPositionGeometry>;

std::pair<MeshP, GeomP > uniform_refine(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, int refine){
    AdaptiveTriangulation atri(mesh,geom);
    for (int i = 0; i < refine; ++i) {
        std::vector<Face> faces; for (Face f: atri.mesh().faces()) faces.push_back(f);
        atri.refine(faces);
    }
    atri.mesh().compress();
    atri.intrinsicTriangulation().refreshQuantities();
    auto nMesh = atri.mesh().copy();
    auto nGeom = std::make_unique<VertexPositionGeometry>(*nMesh,intrinsic_geom(atri.intrinsicTriangulation(),geom).reinterpretTo(*nMesh));
    return { std::move(nMesh), std::move(nGeom)};
}

class TaylorVorticesCase : public TestCase{
  public:
    std::unique_ptr<ManifoldSurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geom;
    std::unique_ptr<AdaptiveFluidSolver> solver;
    Evaluator eval_a, eval_t;

    DOPRI5_sample dp5s;

    DoeflerConf doerflerConf = DoerflerPreset(DoerflerPresetConf::LOW);
    DOPRI5_conf dopri5Conf = DOPRI5Preset(DOPRI5PresetConf::LOW);;
    void registerAll(Evaluator& ev,int h_size){
        ev.reg("dt (s)", [](EvData d) { return d.dp5s.t_past; });
        // ev.reg("attempts", [](EvData d) { return d.dp5s.attempts; });
        ev.reg(R"($\|u\|_2$)", [](EvData d) { return L2Norm(d.vel.u,d.geom); });
        ev.reg(R"($\int_M \psi$)", [](EvData d) { return integral(d.vel.stream_function,d.geom); });
        ev.reg(R"($\int_M w$)", [](EvData d) { return integral(d.wc.w,d.geom); });
        ev.reg(R"($\int_M \frac{d}{dt} w$)", [](EvData d) { return integral(d.rhs.w,d.geom); });
        ev.reg("wall time per simulation time (s)", [](EvData d) { return d.time_per_sim_sec; });
        ev.reg(R"($\eta_R$)", [](EvData d) { return d.poison_residual_error; });
        for (int i = 0; i < h_size; ++i) {
            ev.reg("$c_IDX("+std::to_string(i) +")$", [i](EvData d) {
                return d.wc.c[i];
            });
            ev.reg(R"($\frac{d}{dt} c_IDX()"+std::to_string(i) +")$", [i](EvData d) { return d.rhs.c[i]; });
        }
    }
    void registerStep(Evaluator& ev,int h_size){
        ev.reg(R"($\|u\|_2$)", [](EvData d) { return L2Norm(d.vel.u,d.geom); });
        ev.reg(R"($\int_M \psi$)", [](EvData d) { return integral(d.vel.stream_function,d.geom); });
        ev.reg(R"($\int_M w$)", [](EvData d) { return integral(d.wc.w,d.geom); });
        ev.reg(R"($\int_M \frac{d}{dt} w$)", [](EvData d) { return integral(d.rhs.w,d.geom); });
        // ev.reg(R"($\#F$)", [](EvData d) { return d.mesh.nFaces(); });
        ev.reg(R"($\eta_R$)", [](EvData d) { return d.poison_residual_error; });
        for (int i = 0; i < h_size; ++i) {
            ev.reg("$c_IDX("+std::to_string(i) +")$", [i](EvData d) {
              return d.wc.c[i];
            });
            ev.reg(R"($\frac{d}{dt} c_IDX()"+std::to_string(i) +")$", [i](EvData d) { return d.rhs.c[i]; });
        }
    }

    EvData evData(double time_per_sim_sec){
        ManifoldSurfaceMesh& mesh = solver->tri.mesh();
        IntrinsicGeometryInterface& geom = solver->tri.geom();
        std::vector<FaceData<double>> face_dc (solver->h.size(),FaceData<double>(mesh));
        auto rhs = evalRHS(mesh,geom,solver->wc,solver->h, solver->S,&face_dc);
        auto u = solver->velocity();
        auto errorsqr = poisson_residual_error_sqr(mesh,geom,solver->wc.w,u.stream_function);
        double prsq = 0;
        for (Face f: mesh.faces()) prsq += std::sqrt(errorsqr[f]);
        EvData data(mesh,geom,u, rhs,solver->wc,solver->h,dp5s,face_dc,time_per_sim_sec,prsq);
        return data;
    }

    TaylorVorticesCase(std::string name, std::string label, MeshP&& pmesh, GeomP&& pgeom, AdaptiveFluidSolverData data)
        : TestCase(name,label), mesh(std::move(pmesh)), geom(std::move(pgeom)), solver(std::make_unique<AdaptiveFluidSolver>(*mesh,*geom,data)),
          dopri5Conf(data.dopri5Conf), doerflerConf(data.doerflerConf)
    {
        solver->wc.w = TaylorInitializer().wc(solver->tri.intrinsicTriangulation(),*geom).w;
        registerAll(eval_t, solver->h.size());
        registerAll(eval_a, solver->h.size());
        dp5s = DOPRI5_sample(solver->wc,solver->dt,solver->dt,0);
    }

    void runUntil(double t) override{
        while (solver->elapsed_time < t){
            Stopwatch s;
            dp5s = solver->step();
            auto wall_time = s.stop();
            double tpss = (s.stop() /std::chrono::duration<double>(dp5s.t_past));
            EvData ev  = evData(tpss);
            eval_t.onStep(ev,solver->elapsed_time);
            total_time+= wall_time;
        }
    }

    polyscope::SurfaceMesh* visualize(fs::path f_screenshots) override{
        AdaptiveTriangulation& atri = solver->tri;
        auto vp = intrinsic_geom(atri.intrinsicTriangulation(),*geom);
        auto& m = atri.mesh();
        polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh(name, vp,atri.mesh().getFaceVertexList(), polyscopePermutations(atri.mesh()));
        pm->setEdgeColor(glm::vec3(0,0,0));
        pm->setEdgeWidth(0.3);
        auto g =VertexPositionGeometry(m,vp);

        auto* vsq = pm->addVertexScalarQuantity("vorticity",solver->wc.w);
        vsq->setColorMap("coolwarm");
        vsq->setMapRange({-30,30})->setEnabled(true);
        colormap = { vsq->getColorMap(), vsq->getMapRange() };
        if(!f_screenshots.empty()){
            polyscope::screenshot(f_screenshots/"vorticity.png",true);

            g.requireFaceTangentBasis();
            FaceData<Vector3> e1(m),e2(m);
            for (Face f: m.faces()) { e1[f] = g.faceTangentBasis[f][0], e2[f] = g.faceTangentBasis[f][1];  }

            auto* ftvq = pm->addFaceTangentVectorQuantity("u",solver->velocity().u, e1,e2)->setEnabled(true);
            polyscope::screenshot(f_screenshots/"u.png",true);
            ftvq->setEnabled(false);


            int i = 0;
            for (auto& h:solver->h) {
                std::string name = "h_{"+std::to_string(i) + "}";
                auto* ftq = pm->addFaceTangentVectorQuantity("h"+std::to_string(i++),h, e1,e2);
                ftq->setEnabled(true);
                ftq->setVectorLengthRange(2);
                ftq->setVectorLengthScale(0.04);
                polyscope::screenshot(f_screenshots/(name+".png"),true);
                ftq->setEnabled(false);
            }
            pm->remove();
        }
        return pm;
    }

    std::pair<std::string,std::pair<double,double>> colormap;
    void writeCM(fs::path p){
        std::ofstream ss(p);
        ss << colormap.first << "," << colormap.second.first << "," << colormap.second.second;
        ss.close();
    }

    void write(fs::path p) override {
        to_csv(dopri5Conf,p/ (name+"dopri5.csv"));
        to_csv(doerflerConf,p/ (name+"doerfler.csv"));
        eval_t.saveCSV_T(p/ (name + "-measurements.csv"));
        writeCM(p/"colormap.csv");
    }
};

void refineAndReset(TaylorVorticesCase& aCase, GeomP& geom) {
    //aCase.eval_a = Evaluator();
    //aCase.registerAll(aCase.eval_a, aCase.solver->h.size());
    for (int i = 0; i < 16; ++i) {
        aCase.solver->adapt();
        // aCase.eval_a.onStep(aCase.evData(), 1);
    }
    aCase.solver->wc.w = TaylorInitializer().wc(aCase.solver->tri.intrinsicTriangulation(), *geom).w;
}
std::shared_ptr<TaylorVorticesCase> makeTaylorCase(
    std::string name, std::string label,
    ManifoldSurfaceMesh& pmesh, VertexPositionGeometry& pgeom,
    bool adapt_time, bool adapt_space,
    DOPRI5PresetConf dp5conf = DOPRI5PresetConf::MEDIUM,
    DoerflerPresetConf doerflerConf = DoerflerPresetConf::LOW,
    MARKING_STRATEGY marking_strategy = MARKING_STRATEGY::PATTERN
    )
{

    // deep copy mesh and geom
    std::unique_ptr<ManifoldSurfaceMesh> mesh = pmesh.copy();
    GeomP  geom= std::make_unique<VertexPositionGeometry>(*mesh,pgeom.vertexPositions.reinterpretTo(*mesh));
    // init data
    AdaptiveFluidSolverData data { DOPRI5Preset(dp5conf), DoerflerPreset(doerflerConf), 0.01, adapt_time, adapt_space, marking_strategy };

    auto aCase = std::make_shared<TaylorVorticesCase>(
        name, label, std::move(mesh), std::move(geom), data
    );

    if (adapt_space) {
        refineAndReset(*aCase, aCase->geom);
    }

    return aCase;
}

std::shared_ptr<TaylorVorticesCase> taylorVortices_ASAT(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom) {
    return makeTaylorCase("ASAT", "Adaptive Space Adaptive Time (ASAT)", mesh, geom, true, true);
}

std::shared_ptr<TaylorVorticesCase> taylorVortices_SSST(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom) {
    return makeTaylorCase("SSST","Static Space Static Time (SSST)", mesh, geom, false, false);
}

std::shared_ptr<TaylorVorticesCase> taylorVortices_ASST(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom) {
    return makeTaylorCase("ASST", "Adaptive Space Static Time", mesh, geom, false, true);
}

std::shared_ptr<TaylorVorticesCase> taylorVortices_SSAT(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom) {
    return makeTaylorCase("SSAT","Static Space Adaptive Time (SSAT)", mesh, geom, true, false);
}


std::shared_ptr<TaylorVorticesCase> taylorVortices_OR(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom) {
    return makeTaylorCase("OR", "Original (Yin et al., 2023)", mesh, geom, false, false);
}

class Comparator {
  public:
    std::vector<std::shared_ptr<TestCase>> testcases;

    void runUntil(double t) {
        std::cout << "Running until " << t << std::endl;
        for(auto& tc: testcases){
            std::cout << " - " << tc->name << std::endl;
            tc->runUntil(t);
        }
    }

    int step = 0;
    void visualize(fs::path f_screenshots = fs::path()) {
        polyscope::removeAllStructures();
        for(auto& tc: testcases){
            if (f_screenshots.empty()) {
                tc->visualize();
            }else {
                fs::path f_tc = f_screenshots / (tc->name+std::to_string(step));
                fs::create_directories(f_tc);
                tc->visualize(f_tc);
            }
        }
        step++;
    };

    void csv_total_times(fs::path filename){
        std::ofstream file(filename);
        if (!file.is_open()) { throw std::runtime_error("Cannot open file for writing: " + filename.string()); }

        file << ",wall time (s)" << std::endl;

        for (int i = 0; i < testcases.size(); ++i) {
            file << testcases[i]->name << "," << testcases[i]->total_time.count() << std::endl;
        }
        file.close();
    };
    void csv_comparator(fs::path filename){
        std::ofstream file(filename);
        if (!file.is_open()) { throw std::runtime_error("Cannot open file for writing: " + filename.string()); }
        file << "name,label" << std::endl;
        for (int i = 0; i < testcases.size(); ++i) {
            file << testcases[i]->name << ",\"" << testcases[i]->label << "\"" << std::endl;
        }
        file.close();
    }

    void write(fs::path p) {
        fs::path f_data = p / "data";
        fs::create_directories(f_data);
        for(auto& tc: testcases){
            tc->write(f_data);
        }
        csv_total_times(f_data / "total_time.csv");
        csv_comparator(p / "config.csv");
    }
};
struct CaseFolder{
    fs::path fsrc;
    fs::path fmodels;
    fs::path fev;
    fs::path f_screenshots;
    fs::path flatest;
    CaseFolder(std::string name)
        : fsrc(fs::path(__FILE__).parent_path()),
          fmodels(fsrc / "models" ),
          fev(run_folder(results_folder())),
          f_screenshots(fev / "snapshots"),
          flatest(fev.parent_path() / ("run_"+name)) {
        fs::create_directories(f_screenshots);
    }
};

void init_ps(Comparator& cmp){
    polyscope::init();
    cmp.visualize();
    polyscope::view::ensureViewValid();
    polyscope::view::fov = 30;
    polyscope::view::projectionMode = polyscope::ProjectionMode::Orthographic;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    polyscope::view::setWindowSize(1080,1080);
    polyscope::options::ssaaFactor = 3;
}

TEST(EvaluatorTest, Evaluate)
{
    CaseFolder cf("tc1");

    auto [mesh,geom] = readManifoldSurfaceMesh(cf.fmodels /"cheese_min.stl");
    auto [meshO,geomO] = readManifoldSurfaceMesh(cf.fmodels/"cheese_oriented.stl");

    Comparator cpm;
    cpm.testcases = {
        taylorVortices_OR(*meshO, *geomO),
        // taylorVortices_SSAT(*meshO,*geomO),
        // taylorVortices_SSST_REFINED(*mesh,*geom),
        taylorVortices_ASAT(*mesh,*geom),
    };

    init_ps(cpm);


    cpm.visualize(cf.f_screenshots);
    cpm.runUntil(1.5); cpm.visualize(cf.f_screenshots);
    cpm.runUntil(3); cpm.visualize(cf.f_screenshots);
    cpm.runUntil(4.5); cpm.visualize(cf.f_screenshots);
    cpm.runUntil(6); cpm.visualize(cf.f_screenshots);

    cpm.write(cf.fev);
    copyFolder(cf.fev,cf.flatest);
    cpm.visualize();
    polyscope::show();
}

TEST(EvaluatorTest, EvaluateDiffDopri)
{
    CaseFolder cf ("tc2");
    auto [mesh,geom] = readManifoldSurfaceMesh(cf.fmodels /"cheese_min.stl");
    auto [meshO,geomO] = readManifoldSurfaceMesh(cf.fmodels/"cheese_oriented.stl");
    Comparator cpm;
    cpm.testcases = {
        taylorVortices_OR(*meshO, *geomO),
        makeTaylorCase("low", "Low Precision", *mesh, *geom, true, true,DOPRI5PresetConf::VERY_LOW,DoerflerPresetConf::LOW),
        makeTaylorCase("vhigh", "Very High Precision", *mesh, *geom, true, true,DOPRI5PresetConf::MEDIUM,DoerflerPresetConf::MEDIUM),
    };

    init_ps(cpm);
    cpm.visualize(cf.f_screenshots);
    cpm.runUntil(3); cpm.visualize(cf.f_screenshots);
    cpm.runUntil(6); cpm.visualize(cf.f_screenshots);

    cpm.write(cf.fev);
    copyFolder(cf.fev,cf.flatest);
    cpm.visualize();
    polyscope::show();
}
TEST(EvaluatorTest, EvaluateAdapt)
{
    CaseFolder cf ("tc3");
    auto [mesh,geom] = readManifoldSurfaceMesh(cf.fmodels /"cheese_min.stl");
    Comparator cpm;
    cpm.testcases = {
        makeTaylorCase("adaptive", "Adaptive Refinement", *mesh, *geom, false, false,DOPRI5PresetConf::HIGH,DoerflerPresetConf::VERY_HIGH),
        makeTaylorCase("uniform", "Uniform Refinement", *mesh, *geom, false, false,DOPRI5PresetConf::HIGH,DoerflerPresetConf::UNIFORM_REFINE),
    };
    auto* tc1 = dynamic_cast<TaylorVorticesCase*>(cpm.testcases[0].get());
    auto* tc2 = dynamic_cast<TaylorVorticesCase*>(cpm.testcases[1].get());
    init_ps(cpm);


    Evaluator ev1, ev2;
    tc1->registerStep(ev1,tc1->solver->h.size());
    tc2->registerStep(ev2,tc2->solver->h.size());

    // Step 0
    ev1.onStep(tc1->evData(1),tc1->solver->tri.mesh().nFaces());
    ev2.onStep(tc2->evData(1),tc2->solver->tri.mesh().nFaces());
    cpm.visualize(cf.f_screenshots); std::cout << " - vis" << std::endl;

    auto n = 32768; int i =0;
    std::cout << "adaptive" << std::endl;
    std::vector<int> indices;
    for (int j = 1; j < 5; j++) {
        int idx = round((double)(j * (33 - 1)) / (5 - 1));
        indices.push_back(idx);
    }
    while (tc1->solver->tri.mesh().nFaces() <n){
        // TODO: plot against #faces
        std::cout << " - Refining step " << i <<  "faces: " << tc1->solver->tri.mesh().nFaces() << std::endl;
        tc1->solver->adapt();
        ev1.onStep(tc1->evData(1),tc1->solver->tri.mesh().nFaces());
        bool vis = std::ranges::contains(indices,i);
        if(vis) {
            tc2->solver->adapt();
            ev2.onStep(tc2->evData(1),tc2->solver->tri.mesh().nFaces());
            cpm.visualize(cf.f_screenshots); std::cout << " - vis" << std::endl;
        };
        i++;
    }
    cpm.write(cf.fev);
    ev1.saveCSV_T(cf.fev/"data"/(tc1->name+"-measurements.csv"));
    ev2.saveCSV_T(cf.fev/"data"/(tc2->name+"-measurements.csv"));
    copyFolder(cf.fev,cf.flatest);
    cpm.visualize();
    polyscope::show();
}

TEST(EvaluatorTest,evaluateInitialMarkings) {
    CaseFolder cf ("tc4");
    auto [mesh,geom] = readManifoldSurfaceMesh(cf.fmodels /"cheese_min.stl");
    Comparator cpm;
    cpm.testcases = {
        makeTaylorCase("pattern", "pattern-initialized", *mesh, *geom, false, false,DOPRI5PresetConf::VERY_LOW,DoerflerPresetConf::UNIFORM_REFINE,MARKING_STRATEGY::PATTERN),
        makeTaylorCase("random", "random-initialized", *mesh, *geom, false, false,DOPRI5PresetConf::VERY_LOW,DoerflerPresetConf::UNIFORM_REFINE,MARKING_STRATEGY::RANDOM),
    };
    auto* tc1 = dynamic_cast<TaylorVorticesCase*>(cpm.testcases[0].get());
    auto* tc2 = dynamic_cast<TaylorVorticesCase*>(cpm.testcases[1].get());
    int n = 9000;
    while (tc1->solver->tri.mesh().nFaces()< n) {
        tc1->solver->adapt();
    }
    while (tc2->solver->tri.mesh().nFaces()< n) {
        tc2->solver->adapt();
    }


    init_ps(cpm);
    cpm.visualize(cf.f_screenshots);
    cpm.runUntil(3); cpm.visualize(cf.f_screenshots);
    cpm.runUntil(6); cpm.visualize(cf.f_screenshots);

    cpm.write(cf.fev);
    copyFolder(cf.fev,cf.flatest);
    cpm.visualize();
    polyscope::show();

}

TEST(EvaluatorTest, Clean)
{
    fs::path folder = results_folder();
    if (!fs::exists(folder) || !fs::is_directory(folder)) {
        std::cerr << "Folder does not exist: " << folder << std::endl;
        return;
    }

    for (const auto& entry : fs::directory_iterator(folder)) {
        if (entry.is_directory()) {
            const std::string name = entry.path().filename().string();
            if (name.rfind("run_", 0) == 0) {  // name starts with "run_"
                try {
                    fs::remove_all(entry.path());
                    std::cout << "Removed folder: " << entry.path() << std::endl;
                } catch (const fs::filesystem_error& e) {
                    std::cerr << "Error removing folder " << entry.path() << ": " << e.what() << std::endl;
                }
            }
        }
    }
}

}
