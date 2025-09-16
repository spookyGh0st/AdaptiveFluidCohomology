#include "Evaluator.h"
#include "Stopwatch.h"
#include "TaylorInitializer.h"
#include "eider/homotopy.h"
#include "eider/util.h"
#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include "geometrycentral/surface/meshio.h"
#include <algorithm>
#include <eider/AdaptiveFluidSolver.h>
#include <filesystem>
#include <gtest/gtest.h>
#include <implot.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

namespace fs = std::filesystem;

fs::path results_folder() {
    fs::path file_path = __FILE__;
    fs::path project_root = file_path.parent_path().parent_path().parent_path();
    fs::path results_f = project_root/ "tex" / "thesis"/ "results";
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
        if (!fs::exists(destination)) {
            fs::create_directories(destination);
        }

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
    std::string name;
    TestCase(std::string name) :name(name) {}
    virtual void runUntil(double t) = 0;
    virtual void write(fs::path) = 0;
    virtual polyscope::SurfaceMesh* visualize() = 0;
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
        ev.reg("dt", [](EvData d) { return d.dp5s.t_past; });
        ev.reg("attempts", [](EvData d) { return d.dp5s.attempts; });
        ev.reg("velocity", [](EvData d) { return L2Norm(d.vel.u,d.geom); });
        ev.reg("streamfunction", [](EvData d) { return L2Norm(d.vel.stream_function,d.geom); });
        ev.reg("w", [](EvData d) { return L2Norm(d.wc.w,d.geom); });
        ev.reg("dw", [](EvData d) { return L2Norm(d.rhs.w,d.geom); });
        for (int i = 0; i < h_size; ++i) {
            ev.reg("c"+std::to_string(i), [&i](EvData d) { return d.wc.c[i]; });
            ev.reg("dc"+std::to_string(i), [&i](EvData d) { return d.rhs.c[i]; });
        }
    }

    EvData evData(){
        ManifoldSurfaceMesh& mesh = solver->tri.mesh();
        IntrinsicGeometryInterface& geom = solver->tri.geom();
        std::vector<FaceData<double>> face_dc (solver->h.size(),FaceData<double>(mesh));
        auto rhs = evalRHS(mesh,geom,solver->wc,solver->h, solver->S,&face_dc);
        EvData data(mesh,geom,solver->velocity(), rhs,solver->wc,solver->h,dp5s,face_dc);
        return data;
    }

    TaylorVorticesCase(std::string name, MeshP&& pmesh, GeomP&& pgeom, AdaptiveFluidSolverData data)
        : TestCase(name), mesh(std::move(pmesh)), geom(std::move(pgeom)), solver(std::make_unique<AdaptiveFluidSolver>(*mesh,*geom,data))
    {
        solver->wc.w = TaylorInitializer().wc(solver->tri.intrinsicTriangulation(),*geom).w;
        registerAll(eval_t, solver->h.size());
        registerAll(eval_a, solver->h.size());
        dp5s = DOPRI5_sample(solver->wc,solver->dt,solver->dt,0);
    }

    void runUntil(double t) override{
        while (solver->elapsed_time < t){
            dp5s = solver->step();
            eval_t.onStep(evData(),dp5s.t_past);
        }
    }

    polyscope::SurfaceMesh* visualize() override{
        AdaptiveTriangulation& atri = solver->tri;
        auto vp = intrinsic_geom(atri.intrinsicTriangulation(),*geom);
        polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh(name, vp,atri.mesh().getFaceVertexList(), polyscopePermutations(atri.mesh()));
        pm->addVertexScalarQuantity("vorticity",solver->wc.w)->setMapRange({-30,30})->setEnabled(true);

        return pm;
    }

    void write(fs::path p) override {
        to_csv(dopri5Conf,p/ (name+"-dopri5.csv"));
        to_csv(doerflerConf,p/ (name+"-doerfler.csv"));
        eval_t.saveCSV_T(p/ (name + "-time.csv"));
    }
};

void refineAndReset(TaylorVorticesCase& aCase, GeomP& geom) {
    aCase.eval_a = Evaluator();
    aCase.registerAll(aCase.eval_a, aCase.solver->h.size());
    for (int i = 0; i < 16; ++i) {
        aCase.solver->adapt();
        aCase.eval_a.onStep(aCase.evData(), 1);
    }
    aCase.solver->wc.w = TaylorInitializer().wc(aCase.solver->tri.intrinsicTriangulation(), *geom).w;
}
std::shared_ptr<TaylorVorticesCase> makeTaylorCase(
    std::string name,
    MeshP&& mesh, GeomP&& geom,
    bool adapt_time, bool adapt_space,
    bool uniformRefine = false)
{
    AdaptiveFluidSolverData data {
        DOPRI5Preset(DOPRI5PresetConf::HIGH),
        DoerflerPreset(DoerflerPresetConf::LOW),
        0.01,
        adapt_time,
        adapt_space
    };

    if (uniformRefine) {
        std::tie(mesh, geom) = uniform_refine(*mesh, *geom, 2);
    }

    auto aCase = std::make_shared<TaylorVorticesCase>(
        name, std::move(mesh), std::move(geom), data
    );

    if (adapt_space) {
        refineAndReset(*aCase, aCase->geom);
    }

    return aCase;
}

std::shared_ptr<TaylorVorticesCase> taylorVortices_ASAT(MeshP&& mesh, GeomP&& geom) {
    return makeTaylorCase("Taylor Vortices AS AT", std::move(mesh), std::move(geom), true, true);
}

std::shared_ptr<TaylorVorticesCase> taylorVortices_SSST(MeshP&& mesh, GeomP&& geom) {
    return makeTaylorCase("Taylor Vortices SS ST", std::move(mesh), std::move(geom), false, false);
}

std::shared_ptr<TaylorVorticesCase> taylorVortices_ASST(MeshP&& mesh, GeomP&& geom) {
    return makeTaylorCase("Taylor Vortices AS ST", std::move(mesh), std::move(geom), false, true);
}

std::shared_ptr<TaylorVorticesCase> taylorVortices_SSAT(MeshP&& mesh, GeomP&& geom) {
    return makeTaylorCase("Taylor Vortices SS AT", std::move(mesh), std::move(geom), true, false);
}

std::shared_ptr<TaylorVorticesCase> taylorVortices_SSST_REFINED(MeshP&& mesh, GeomP&& geom) {
    return makeTaylorCase("Taylor Vortices SS ST REFINED", std::move(mesh), std::move(geom), false, false, true);
}

std::shared_ptr<TaylorVorticesCase> taylorVortices_OR(MeshP&& mesh, GeomP&& geom) {
    return makeTaylorCase("Taylor Vortices Original", std::move(mesh), std::move(geom), false, false);
}

class Comparator {
  public:
    std::vector<std::shared_ptr<TestCase>> testcases;

    void runUntil(double t) {
        for(auto& tc: testcases){
            tc->runUntil(t);
        }
    }

    void visualize() {
        double height = 0;
        for(auto& tc: testcases){
            polyscope::SurfaceMesh* pm = tc->visualize();
            pm->setEdgeColor(glm::vec3(0,0,0));
            pm->setEdgeWidth(0.3);
            pm->setPosition(glm::vec3( height,0,0));
            height += 2.2;
        }
    };

    void write(fs::path p) {
        for(auto& tc: testcases){
            tc->write(p);
        }
    }
};

TEST(EvaluatorTest, Evaluate)
{
    results_folder();
    fs::path fsrc = fs::path(__FILE__).parent_path();
    fs::path fds = fsrc / "models" / "cheese_min.stl";
    fs::path fev = run_folder(results_folder());
    fs::path flatest = fev.parent_path() / "run_latest";

    auto [mesh,geom] = readManifoldSurfaceMesh(fds.string());
    auto [meshO,geomO] = readManifoldSurfaceMesh((fsrc/"models"/"cheese_oriented.stl"));

    Comparator cpm;
    cpm.testcases = {
        taylorVortices_OR(meshO->copy(), geomO->copy()),
        taylorVortices_SSST(mesh->copy(), geom->copy()),
        taylorVortices_SSST_REFINED(mesh->copy(), geom->copy()),
        taylorVortices_ASAT(mesh->copy(), geom->copy()),
    };

    polyscope::init();
    cpm.visualize();
    polyscope::view::ensureViewValid();
    polyscope::view::fov = 30;
    polyscope::view::projectionMode = polyscope::ProjectionMode::Orthographic;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    polyscope::view::setWindowSize(1080*cpm.testcases.size()/2,1080/2);

    polyscope::options::ssaaFactor = 3;

    cpm.visualize();
    polyscope::screenshot(fev/"adaptive_a1.png",true);

    cpm.runUntil(1.5); cpm.visualize();
    polyscope::screenshot(fev/"adaptive_a2.png",true);

    cpm.runUntil(3); cpm.visualize();
    polyscope::screenshot(fev/"adaptive_a3.png",true);

    cpm.runUntil(4.5); cpm.visualize();
    polyscope::screenshot(fev/"adaptive_a4.png",true);

    cpm.runUntil(6); cpm.visualize();
    polyscope::screenshot(fev/"adaptive_a5.png",true);

    cpm.write(fev);
    copyFolder(fev,flatest);
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
