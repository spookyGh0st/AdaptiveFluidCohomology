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
    virtual ManifoldSurfaceMesh& getMesh() = 0;
};

using MeshP = std::unique_ptr<ManifoldSurfaceMesh>;
using GeomP = std::unique_ptr<VertexPositionGeometry>;

std::pair<MeshP, GeomP > uniform_refine(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, int refine,MARKING_STRATEGY strategy = MARKING_STRATEGY::PATTERN){
    AdaptiveTriangulation atri(mesh,geom,strategy);
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
std::pair<MeshP, GeomP > delauny_refine(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, double degree = 25, double h_size = std::numeric_limits<double>::infinity()){
    IntegerCoordinatesIntrinsicTriangulation atri(mesh,geom);
    atri.delaunayRefine(degree);
    atri.intrinsicMesh->compress();
    atri.refreshQuantities();
    auto nMesh = atri.intrinsicMesh->copy();
    auto nGeom = std::make_unique<VertexPositionGeometry>(*nMesh,intrinsic_geom(atri,geom).reinterpretTo(*nMesh));
    return { std::move(nMesh), std::move(nGeom)};
}

class TaylorVorticesCase : public TestCase{
  public:
    std::unique_ptr<ManifoldSurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geom;
    std::unique_ptr<AdaptiveFluidSolver> solver;
    Evaluator eval_a, eval_t;
    bool fix_c = false;

    DOPRI5_sample dp5s;

    DoeflerConf doerflerConf = DoerflerPreset(DoerflerPresetConf::LOW);
    DOPRI5_conf dopri5Conf = DOPRI5Preset(DOPRI5PresetConf::LOW);;
    void registerAll(Evaluator& ev,int h_size){
        registerProperties(ev,defaultTimeProperties(),h_size);
    }
    void registerStep(Evaluator& ev,int h_size){
        ev.reg(R"($\|u\|_2$)", [](EvData d) { return L2Norm(d.vel.u,d.geom); });
        // ev.reg(R"($\int_M \psi$)", [](EvData d) { return integral(d.vel.stream_function,d.geom); });
        // ev.reg(R"($\int_M w$)", [](EvData d) { return integral(d.wc.w,d.geom); });
        ev.reg(R"($\int_M \frac{d}{dt} w$)", [](EvData d) { return integral(d.rhs.w,d.geom); });
        // ev.reg(R"($\#F$)", [](EvData d) { return d.mesh.nFaces(); });
        ev.reg(R"($\eta_R$)", [](EvData d) { return d.poison_residual_error; });
        for (int i = 0; i < h_size; ++i) {
            // # ev.reg("$c_IDX("+std::to_string(i) +")$", [i](EvData d) {
            // #   return d.wc.c[i];
            // # });
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
        for (Face f: mesh.faces()) prsq += errorsqr[f];
        prsq = std::sqrt(prsq);
        EvData data(mesh,geom,u, rhs,solver->wc,solver->h,dp5s,face_dc,time_per_sim_sec,prsq,solver->h_interpolated,intrinsic_geom(solver->tri.intrinsicTriangulation(),*this->geom));
        return data;
    }

    TaylorVorticesCase(std::string name, std::string label, MeshP&& pmesh, GeomP&& pgeom, AdaptiveFluidSolverData data, bool fix_c = false)
        : TestCase(name,label), mesh(std::move(pmesh)), geom(std::move(pgeom)), solver(std::make_unique<AdaptiveFluidSolver>(*mesh,*geom,data)),
          dopri5Conf(data.dopri5Conf), doerflerConf(data.doerflerConf), fix_c(fix_c)
    {
        solver->wc.w = TaylorInitializer().wc(solver->tri.intrinsicTriangulation(),*geom).w;
        registerAll(eval_t, solver->h.size());
        registerAll(eval_a, solver->h.size());
        dp5s = DOPRI5_sample(solver->wc,solver->dt,solver->dt,0);
    }

    void runUntil(double t) override{
        while (solver->elapsed_time < t){
            auto c_old = solver->wc.c;
            Stopwatch s;
            dp5s = solver->step();
            auto wall_time = s.stop();
            if (fix_c) solver->wc.c = c_old;
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
            auto fullh = solver->hom.fullHarmonicBasis();
            for (auto& h:solver->h) {
                std::string name = "h_{"+std::to_string(i) + "}";
                auto* ftq = pm->addFaceTangentVectorQuantity("h"+std::to_string(i),h, e1,e2);
                ftq->setEnabled(true);
                ftq->setVectorLengthRange(2);
                ftq->setVectorLengthScale(0.04);
                polyscope::screenshot(f_screenshots/(name+".png"),true);
                ftq->setEnabled(false);

                name = "h-df_"+std::to_string(i);
                auto* esq = pm->addEdgeScalarQuantity(name,fullh.df[i]);
                esq->setEnabled(true);
                polyscope::screenshot(f_screenshots/(name+".png"),true);
                esq->setEnabled(false);

                name = "h-proj_df_"+std::to_string(i);
                esq = pm->addEdgeScalarQuantity(name,fullh.proj_df[i]);
                esq->setEnabled(true);
                polyscope::screenshot(f_screenshots/(name+".png"),true);
                esq->setEnabled(false);

                i++;

            }
            pm->remove();
        }
        return pm;
    }

    ManifoldSurfaceMesh& getMesh() override { return solver->tri.mesh(); }

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
    AdaptiveFluidSolverData d, bool fix_c = false
    )
{

    // deep copy mesh and geom
    std::unique_ptr<ManifoldSurfaceMesh> mesh = pmesh.copy();
    GeomP  geom= std::make_unique<VertexPositionGeometry>(*mesh,pgeom.vertexPositions.reinterpretTo(*mesh));
    // init data

    auto aCase = std::make_shared<TaylorVorticesCase>( name, label, std::move(mesh), std::move(geom), d, fix_c);

    if (d.adaptive_space) { refineAndReset(*aCase, aCase->geom); }

    return aCase;
}


std::shared_ptr<TaylorVorticesCase> taylorVortices_ASAT(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom) {
    AdaptiveFluidSolverData data(DOPRI5PresetConf::MEDIUM,DoerflerPresetConf::LOW,0.01,true,true,MARKING_STRATEGY::PATTERN,false);
    return makeTaylorCase("ASAT", "Adaptive Space Adaptive Time (ASAT)", mesh, geom, data);
}

std::shared_ptr<TaylorVorticesCase> taylorVortices_OR(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom) {
    AdaptiveFluidSolverData data(DOPRI5PresetConf::MEDIUM,DoerflerPresetConf::LOW,0.01,false,false,MARKING_STRATEGY::PATTERN,false);
    return makeTaylorCase("OR", "Original (Yin et al., 2023)", mesh, geom, data);
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

    void csv_faces(fs::path filename){
        std::ofstream file(filename);
        if (!file.is_open()) { throw std::runtime_error("Cannot open file for writing: " + filename.string()); }

        file << ",nfaces" << std::endl;

        for (int i = 0; i < testcases.size(); ++i) {
            file << testcases[i]->name << "," << testcases[i]->getMesh().nFaces() << std::endl;
        }
        file.close();
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
        csv_faces(f_data / "nfaces.csv");
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
    auto [meshO, geomO] = uniform_refine(*mesh,*geom,4);
    // auto [meshO,geomO] = readManifoldSurfaceMesh(cf.fmodels/"cheese_oriented.stl");

    Comparator cpm;
    AdaptiveFluidSolverData data_comp_h(DOPRI5PresetConf::LOW,DoerflerPresetConf::LOW,0.01,true,true,MARKING_STRATEGY::PATTERN,false,false);
    AdaptiveFluidSolverData data_interp_ha= data_comp_h; data_interp_ha.interpolate_harmonic_basis = true; data_interp_ha.use_interpolated_harmonic_basis =true;
    cpm.testcases = {
        taylorVortices_OR(*meshO, *geomO),
        makeTaylorCase("AR", "Adaptive, recomputed h (AR)", *mesh, *geom, data_comp_h),
        makeTaylorCase("AI", "Adaptive, interpolated h (AI)", *mesh, *geom, data_interp_ha),
    };
    init_ps(cpm);


    cpm.visualize(cf.f_screenshots);
    cpm.runUntil(1.5); cpm.visualize(cf.f_screenshots);
    cpm.runUntil(3.0); cpm.visualize(cf.f_screenshots);
    cpm.runUntil(4.5); cpm.visualize(cf.f_screenshots);
     cpm.runUntil(6); cpm.visualize(cf.f_screenshots);

    cpm.write(cf.fev);
    copyFolder(cf.fev,cf.flatest);
    cpm.visualize();
    polyscope::show();
}

TEST(EvaluatorTest, EvaluateLongTerm)
{
    CaseFolder cf("tc8");
    auto [mesh,geom] = readManifoldSurfaceMesh(cf.fmodels /"cheese_min.stl");

    Comparator cpm;
    AdaptiveFluidSolverData data_comp_h(DOPRI5PresetConf::MEDIUM,DoerflerPresetConf::LOW,0.01,false,true,MARKING_STRATEGY::PATTERN,false,false);
    AdaptiveFluidSolverData data_interp_ha= data_comp_h; data_interp_ha.interpolate_harmonic_basis = true; data_interp_ha.use_interpolated_harmonic_basis =true;
    auto tc1 = makeTaylorCase("AR", "Adaptive, recomputed h (AR)", *mesh, *geom, data_comp_h);
    auto tc2 = makeTaylorCase("AI", "Adaptive, interpolated h (AI)", *mesh, *geom, data_interp_ha);
    ExportProperty p = EXPORT_int_w | EXPORT_C| EXPORT_velocity | EXPORT_Vort_Y;
    tc1->eval_t.y.clear(); registerProperties(tc1->eval_t,p,tc1->solver->h.size());
    tc2->eval_t.y.clear(); registerProperties(tc2->eval_t,p,tc2->solver->h.size());
    cpm.testcases = { taylorVortices_OR(*mesh, *geom), tc1,tc2  }; // TODO: Uniform Refine Original
    init_ps(cpm);

    cpm.visualize(cf.f_screenshots);
    cpm.runUntil(10); cpm.visualize(cf.f_screenshots);

    cpm.write(cf.fev);
    copyFolder(cf.fev,cf.flatest);
    cpm.visualize();
    polyscope::show();
}

TEST(EvaluatorTest, EvaluateFixC)
{
    CaseFolder cf("tc6");
    auto [mesh,geom] = readManifoldSurfaceMesh(cf.fmodels/"cheese_oriented.stl");

    Comparator cpm;
    AdaptiveFluidSolverData data(DOPRI5PresetConf::MEDIUM,DoerflerPresetConf::LOW,0.01,false,false,MARKING_STRATEGY::PATTERN,false);
    cpm.testcases = {
        makeTaylorCase("constant", "Constant harmonic coefficient", *mesh, *geom, data,true),
        makeTaylorCase("dynamic", "dynamic harmonic coefficient", *mesh, *geom, data,false),
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
    auto [meshO, geomO] = uniform_refine(*mesh,*geom,4);
    Comparator cpm;
    cpm.testcases = {
        taylorVortices_OR(*meshO, *geomO),
        makeTaylorCase("low", "Low Precision", *mesh, *geom, AdaptiveFluidSolverData(DOPRI5PresetConf::LOW,DoerflerPresetConf::LOW,0.01,true,true,MARKING_STRATEGY::PATTERN,false,false)),
        makeTaylorCase("vhigh", "Very High Precision", *mesh, *geom, AdaptiveFluidSolverData(DOPRI5PresetConf::HIGH,DoerflerPresetConf::MEDIUM,0.01,true,true,MARKING_STRATEGY::PATTERN,false,false)),
    };

    init_ps(cpm);
    cpm.visualize(cf.f_screenshots);
    cpm.runUntil(1);
    cpm.runUntil(2);
    cpm.runUntil(3); cpm.visualize(cf.f_screenshots);
    cpm.runUntil(6); cpm.visualize(cf.f_screenshots);

    cpm.write(cf.fev);
    copyFolder(cf.fev,cf.flatest);
    cpm.visualize();
    polyscope::show();
}
TEST(EvaluatorTest, EvaluateAdapt)
{
    CaseFolder cf ("tc3.1");
    CaseFolder cf2 ("tc3.2");
    auto [mesh,geom] = readManifoldSurfaceMesh(cf.fmodels /"cheese_min.stl");
    Comparator cpm, cmp2;
    auto tc1 = makeTaylorCase("AR", "Adaptive Refinement (recomputed harmonic forms)", *mesh, *geom, AdaptiveFluidSolverData(DOPRI5PresetConf::HIGH,DoerflerPresetConf::VERY_HIGH,0.01,false,false,MARKING_STRATEGY::PATTERN,false,false));
    auto tc2 = makeTaylorCase("UN", "Uniform Refinement", *mesh, *geom, AdaptiveFluidSolverData(DOPRI5PresetConf::HIGH,DoerflerPresetConf::UNIFORM_REFINE,0.01,false,false));
    auto tc3 = makeTaylorCase("AI", "Adaptive Refinement (interpolated harmonic forms)", *mesh, *geom, AdaptiveFluidSolverData(DOPRI5PresetConf::HIGH,DoerflerPresetConf::VERY_HIGH,0.01,false,false,MARKING_STRATEGY::PATTERN,true,true));
    cpm.testcases = { tc2, tc1 }; cmp2.testcases = {tc2,tc3 };
    init_ps(cpm); init_ps(cmp2);

    Evaluator ev1, ev2, ev3;
    tc1->registerStep(ev1,tc1->solver->h.size());
    tc2->registerStep(ev2,tc2->solver->h.size());
    tc3->registerStep(ev3,tc3->solver->h.size());

    // Step 0
    ev1.onStep(tc1->evData(1),tc1->solver->tri.mesh().nFaces());
    ev2.onStep(tc2->evData(1),tc2->solver->tri.mesh().nFaces());
    ev3.onStep(tc3->evData(1),tc3->solver->tri.mesh().nFaces());
    cpm.visualize(cf.f_screenshots);
    cmp2.visualize(cf.f_screenshots);

    auto n = 35000; int i =0;
    while (tc2->solver->tri.mesh().nFaces() < n) {
        std::cout << " - Refining step " << i <<  "faces: " << tc2->solver->tri.mesh().nFaces() << std::endl;
        tc2->solver->adapt();
        ev2.onStep(tc2->evData(1), tc2->solver->tri.mesh().nFaces());

        while (tc1->solver->tri.mesh().nFaces() < tc2->solver->tri.mesh().nFaces()) {
            tc1->solver->adapt();
            ev1.onStep(tc1->evData(1), tc1->solver->tri.mesh().nFaces());
        }
        while (tc3->solver->tri.mesh().nFaces() < tc2->solver->tri.mesh().nFaces()) {
            tc3->solver->adapt();
            ev3.onStep(tc3->evData(1), tc3->solver->tri.mesh().nFaces());
        }
        if(i%2 == 0){
            cpm.visualize(cf.f_screenshots);
            cmp2.visualize(cf2.f_screenshots);
        }
        i++;
    }
    cpm.write(cf.fev);
    ev2.saveCSV_T(cf.fev/"data"/(tc2->name+"-measurements.csv"));
    ev1.saveCSV_T(cf.fev/"data"/(tc1->name+"-measurements.csv"));
    copyFolder(cf.fev,cf.flatest);

    cmp2.write(cf2.fev);
    ev2.saveCSV_T(cf2.fev/"data"/(tc2->name+"-measurements.csv"));
    ev3.saveCSV_T(cf2.fev/"data"/(tc3->name+"-measurements.csv"));
    copyFolder(cf2.fev,cf2.flatest);

    cpm.visualize();
    polyscope::show();
}
TEST(EvaluatorTest, EvaluatePerformance1)
{
    CaseFolder cf("tc8");

    auto [mesh,geom] = readManifoldSurfaceMesh(cf.fmodels /"cheese_min.stl");

    AdaptiveFluidSolverData data_comp_h(DOPRI5PresetConf::LOW,DoerflerPresetConf::LOW,0.01,true,true,MARKING_STRATEGY::PATTERN,false,false);
    AdaptiveFluidSolver solver(*mesh,*geom, data_comp_h);
    solver.wc.w = TaylorInitializer().wc(solver.tri.intrinsicTriangulation(),*geom).w;
    while (solver.elapsed_time < 1) {
        solver.step();
    }
}
TEST(EvaluatorTest, EvaluatePerformance2)
{
    CaseFolder cf("tc8");

    auto [mesh,geom] = readManifoldSurfaceMesh(cf.fmodels /"cheese_min.stl");

    AdaptiveFluidSolverData data_comp_h(DOPRI5PresetConf::LOW,DoerflerPresetConf::LOW,0.01,true,true,MARKING_STRATEGY::PATTERN,true,true);
    AdaptiveFluidSolver solver(*mesh,*geom, data_comp_h);
    solver.wc.w = TaylorInitializer().wc(solver.tri.intrinsicTriangulation(),*geom).w;
    while (solver.elapsed_time < 1) {
        solver.step();
    }
}


// x, y in [0, 1]
double double_shear_layer(double x, double y) {
    const double delta = 0.05; // Perturbation magnitude
    const double rho   = 100.0; // Shear layer thickness/sharpness

    auto sech2 = [](double z) { double t = 1.0 / std::cosh(z); return t * t; };

    // Two opposing vorticity strips perturbed by a cosine wave
    return delta * std::cos(2.0 * M_PI * x)
           + rho * sech2(rho * (y - 0.25))
           - rho * sech2(rho * (y - 0.75));
}

TEST(VideoTest,TorusAnimation) {
    CaseFolder cf("tc9");
    auto [mesh,geom,param] = readParameterizedManifoldSurfaceMesh(cf.fmodels /"torus.obj");
    VertexData<Vector2> uv(*mesh); for (Vertex v: mesh->vertices()) { uv[v] = (*param)[v.corner()]; }

    AdaptiveFluidSolverData data_comp_h(DOPRI5PresetConf::HIGH,DoerflerPresetConf::LOW,0.001,true, true,MARKING_STRATEGY::LONGEST_EDGE,false,false);
    data_comp_h.doerflerConf.threshold_refine = 5;
    data_comp_h.doerflerConf.threshold_coarse = 1;
    AdaptiveFluidSolver solver(*mesh,*geom, data_comp_h);

    for (int i = 0; i < 16; ++i) {
        for (Vertex v: solver.tri.mesh().vertices()){
            auto vec = solver.tri.intrinsicTriangulation().vertexLocations[v].interpolate(uv);
            solver.wc.w[v] = double_shear_layer(vec.y,vec.x);
        }
        solver.adapt();
    }


    polyscope::init();

    bool running = false;
    polyscope::state::userCallback = [&]()
    {
        ImGui::Checkbox("run", &running);
        if(running) {
            solver.step();
            solver.tri.mesh().compress();
            solver.tri.geom().refreshQuantities();
            auto *pm = polyscope::registerSurfaceMesh("mesh", intrinsic_geom(solver.tri.intrinsicTriangulation(), *geom), solver.tri.mesh().getFaceVertexList(), polyscopePermutations(solver.tri.mesh()));
            pm->addVertexScalarQuantity("vorticity", solver.wc.w)->setEnabled(true);
            auto u = solver.velocity();
            auto errorsqr = poisson_residual_error_sqr(solver.tri.mesh(), solver.tri.geom(), solver.wc.w, u.stream_function);
            pm->addFaceScalarQuantity("errorsqr", errorsqr);;
        }

    };
    polyscope::show();

}

TEST(EvaluatorTest,evaluateInitialMarkings) {
    CaseFolder cf ("tc4");
    auto [mesh,geom] = readManifoldSurfaceMesh(cf.fmodels /"cheese_min.stl");
    Comparator cpm;
    cpm.testcases = {
        makeTaylorCase("pattern", "pattern-initialized", *mesh, *geom, AdaptiveFluidSolverData(DOPRI5PresetConf::VERY_LOW,DoerflerPresetConf::UNIFORM_REFINE,0.01,false,false,MARKING_STRATEGY::PATTERN,true)),
        makeTaylorCase("random", "random-initialized", *mesh, *geom, AdaptiveFluidSolverData(DOPRI5PresetConf::VERY_LOW,DoerflerPresetConf::UNIFORM_REFINE,0.01,false,false,MARKING_STRATEGY::RANDOM,true)),
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
    tc1->solver->wc.w = TaylorInitializer().wc(tc1->solver->tri.intrinsicTriangulation(),*geom).w;
    tc2->solver->wc.w = TaylorInitializer().wc(tc2->solver->tri.intrinsicTriangulation(),*geom).w;


    init_ps(cpm);
    cpm.visualize(cf.f_screenshots);
    cpm.runUntil(2); cpm.visualize(cf.f_screenshots);
    cpm.runUntil(4); cpm.visualize(cf.f_screenshots);

    cpm.write(cf.fev);
    copyFolder(cf.fev,cf.flatest);
    cpm.visualize();
    polyscope::show();

}

TEST(EvaluatorTest,evaluateBadTriangulation) {
    CaseFolder cf ("tc7");
    auto [mesh,geom] = readManifoldSurfaceMesh(cf.fmodels /"cheese_min.stl");
    std::tie(mesh,geom) = uniform_refine(*mesh,*geom,1,MARKING_STRATEGY::RANDOM);
    auto [meshO,geomO] = uniform_refine(*mesh,*geom,2,MARKING_STRATEGY::RANDOM);
    AdaptiveFluidSolverData staticD(DOPRI5PresetConf::LOW,DoerflerPresetConf::UNIFORM_REFINE,0.001,false,false,MARKING_STRATEGY::RANDOM,false);
    AdaptiveFluidSolverData adaptDC = staticD;
    adaptDC.dt = 0.02; adaptDC.strategy = MARKING_STRATEGY::LONGEST_EDGE; adaptDC.adaptive_space =true; adaptDC.adaptive_time=true; adaptDC.doerflerConf = DoerflerPreset(DoerflerPresetConf::LOW);
    AdaptiveFluidSolverData adaptDI = adaptDC;
    adaptDI.interpolate_harmonic_basis =true; adaptDI.use_interpolated_harmonic_basis = true;

    Comparator cpm;
    cpm.testcases = {
        makeTaylorCase("OR","Original (Yin et al., 2023)", *mesh, *geom, staticD),
        makeTaylorCase("ASAT", "Adaptive, recomputed h", *mesh, *geom, adaptDC),
        makeTaylorCase("ASATIH", "Adaptive, interpolated h", *mesh, *geom, adaptDI),

    };

    auto* tc1 = dynamic_cast<TaylorVorticesCase*>(cpm.testcases[0].get());
    int n = 9000;
    while (tc1->solver->tri.mesh().nFaces()< n) {
        tc1->solver->adapt();
    }

    init_ps(cpm);
    cpm.visualize(cf.f_screenshots);
    cpm.runUntil(1./4* 3); cpm.visualize(cf.f_screenshots);
    cpm.runUntil(2./4 * 3); cpm.visualize(cf.f_screenshots);
    cpm.runUntil(3./4 * 3); cpm.visualize(cf.f_screenshots);
    cpm.runUntil(4./4 * 3); cpm.visualize(cf.f_screenshots);

    cpm.write(cf.fev);
    copyFolder(cf.fev,cf.flatest);
    cpm.visualize();
    polyscope::show();

}


class ClosedLambCase : public TestCase{
  public:
    std::unique_ptr<ManifoldSurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geom;
    std::unique_ptr<AdaptiveFluidSolver> solver;
    Evaluator eval_t;

    DOPRI5_sample dp5s;

    DoeflerConf doerflerConf = DoerflerPreset(DoerflerPresetConf::LOW);
    DOPRI5_conf dopri5Conf = DOPRI5Preset(DOPRI5PresetConf::LOW);;
    void registerAll(Evaluator& ev,int h_size){
        // ev.reg("dt (s)", [](EvData d) { return d.dp5s.t_past; });
        // ev.reg("attempts", [](EvData d) { return d.dp5s.attempts; });
        // ev.reg(R"($\|u\|_2$)", [](EvData d) { return L2Norm(d.vel.u,d.geom); });
        // ev.reg(R"($\int_M \psi$)", [](EvData d) { return integral(d.vel.stream_function,d.geom); });
        // ev.reg(R"($\int_M w$)", [](EvData d) { return integral(d.wc.w,d.geom); });
        // ev.reg(R"($\int_M \frac{d}{dt} w$)", [](EvData d) { return integral(d.rhs.w,d.geom); });
        // ev.reg("wall time per simulation time (s)", [](EvData d) { return d.time_per_sim_sec; });
        // ev.reg(R"($\eta_R$)", [](EvData d) { return d.poison_residual_error; });
        for (int i = 0; i < h_size; ++i) {
            ev.reg("$c_IDX("+std::to_string(i) +")$", [i](EvData d) {
                return d.wc.c[i];
            });
            // ev.reg(R"($\frac{d}{dt} c_IDX()"+std::to_string(i) +")$", [i](EvData d) { return d.rhs.c[i]; });
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

    ClosedLambCase(std::string name, std::string label, MeshP&& pmesh, GeomP&& pgeom, AdaptiveFluidSolverData data)
        : TestCase(name,label), mesh(std::move(pmesh)), geom(std::move(pgeom)), solver(std::make_unique<AdaptiveFluidSolver>(*mesh,*geom,data)),
          dopri5Conf(data.dopri5Conf), doerflerConf(data.doerflerConf)
    {
        solver->wc.w = VertexData<double>(solver->tri.mesh(), 1);
        solver->wc.c.resize(2,0);
        solver->wc.c[0] = 0.5;
        registerAll(eval_t, solver->h.size());
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
            auto fullh = solver->hom.fullHarmonicBasis();
            for (auto& h:solver->h) {
                std::string name = "h_{"+std::to_string(i) + "}";
                auto* ftq = pm->addFaceTangentVectorQuantity("h"+std::to_string(i),h, e1,e2);
                ftq->setEnabled(true);
                ftq->setVectorLengthRange(2);
                ftq->setVectorLengthScale(0.04);
                polyscope::screenshot(f_screenshots/(name+".png"),true);
                ftq->setEnabled(false);

                // name = "h-df_"+std::to_string(i);
                // auto* esq = pm->addEdgeScalarQuantity(name,fullh.df[i]);
                // esq->setEnabled(true);
                // polyscope::screenshot(f_screenshots/(name+".png"),true);
                // esq->setEnabled(false);

                // name = "h-proj_df_"+std::to_string(i);
                // esq = pm->addEdgeScalarQuantity(name,fullh.proj_df[i]);
                // esq->setEnabled(true);
                // polyscope::screenshot(f_screenshots/(name+".png"),true);
                // esq->setEnabled(false);

                i++;

            }
            pm->remove();
        }
        return pm;
    }

    ManifoldSurfaceMesh& getMesh() override { return solver->tri.mesh(); }

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

std::shared_ptr<ClosedLambCase> makeClosedLambCase(
    std::string name, std::string label,
    ManifoldSurfaceMesh& pmesh, VertexPositionGeometry& pgeom,
    AdaptiveFluidSolverData d
    )
{
    // deep copy mesh and geom
    std::unique_ptr<ManifoldSurfaceMesh> mesh = pmesh.copy();
    GeomP  geom= std::make_unique<VertexPositionGeometry>(*mesh,pgeom.vertexPositions.reinterpretTo(*mesh));
    auto aCase = std::make_shared<ClosedLambCase>( name, label, std::move(mesh), std::move(geom), d );
    return aCase;
}

TEST(EvaluatorTest,evaluateBoundedTorus) {
    CaseFolder cf("tc5");
    // auto [mesh,geom] = readManifoldSurfaceMesh(cf.fmodels /"torus_bounded.obj");
    auto [mesh,geom] = readManifoldSurfaceMesh(cf.fmodels /"torus_bounded.stl");
    Comparator cpm;

    cpm.testcases = {
        makeClosedLambCase("Original","original",*mesh,*geom,AdaptiveFluidSolverData (DOPRI5PresetConf::MEDIUM,DoerflerPresetConf::LOW,0.01,false,false,MARKING_STRATEGY::PATTERN,false)),
        // makeClosedLambCase("Adaptive","adaptive",*mesh2,*geom2,AdaptiveFluidSolverData (DOPRI5PresetConf::LOW,DoerflerPresetConf::LOW,0.01,true,true,MARKING_STRATEGY::PATTERN,false)),
    };

    polyscope::init();
    cpm.visualize();
    polyscope::view::ensureViewValid();
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    polyscope::view::setWindowSize(1080,1080);
    polyscope::options::ssaaFactor = 3;
    polyscope::view::lookAt(glm::vec3{-3., 3., 2.}, glm::vec3{0., 0.5, 0.5});
    cpm.visualize(cf.f_screenshots);
    cpm.runUntil(1);
    cpm.runUntil(10); cpm.visualize(cf.f_screenshots);
    cpm.runUntil(20); cpm.visualize(cf.f_screenshots);
    cpm.runUntil(30); cpm.visualize(cf.f_screenshots);
    cpm.runUntil(75);

    cpm.write(cf.fev);
    copyFolder(cf.fev,cf.flatest);
    cpm.visualize();
    polyscope::show();
}


TEST(CleanEvaluatorTest, Clean)
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
