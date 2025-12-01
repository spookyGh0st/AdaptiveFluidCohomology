#include "TaylorInitializer.h"
#include "eider/util.h"
#include "eider/Stopwatch.h"
#include <algorithm>
#include <eider/AdaptiveFluidSolver.h>
#include <filesystem>
#include <gtest/gtest.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>



namespace geometrycentral::surface{
namespace fs = std::filesystem;
fs::path tmp_dir(){
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&t), "%Y-%m-%dT%H-%M-%S");
    fs::path test_case = std::filesystem::temp_directory_path() /  ss.str();
    fs::create_directory(test_case);
    return test_case;
}

VertexData<Vector3> intrinsic_vertex_positions(IntrinsicTriangulation& Tri, VertexPositionGeometry& inputG){
    auto& mesh = *Tri.intrinsicMesh;
    VertexData<Vector3> int_positions(mesh) ;
    for (Vertex v : mesh.vertices()) { int_positions[v] = Tri.vertexLocations[v].interpolate(inputG.vertexPositions); }
    return int_positions;
}

void visualize_w(AdaptiveFluidSolver& solver, VertexPositionGeometry& geom){
    AdaptiveTriangulation& atri = solver.tri;
    auto vp = intrinsic_vertex_positions(atri.intrinsicTriangulation(),geom);
    auto& m = atri.mesh();
    polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh("mesh", vp,atri.mesh().getFaceVertexList(), polyscopePermutations(atri.mesh()));
    pm->setEdgeColor(glm::vec3(0,0,0));
    pm->setEdgeWidth(0.4);
    auto* vsq = pm->addVertexScalarQuantity("vorticity",solver.wc.w);
    vsq->setColorMap("coolwarm");
    vsq->setMapRange({-30,30});
    vsq->setEnabled(true);
}


void createVideo(AdaptiveFluidSolver& solver, VertexPositionGeometry& geom, const std::string& name, double target_time,std::function<void(AdaptiveFluidSolver&)> step){
    visualize_w(solver,geom);
    polyscope::view::ensureViewValid();
    polyscope::view::fov = 30;
    polyscope::view::projectionMode = polyscope::ProjectionMode::Orthographic;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    polyscope::view::setWindowSize(1080,1080);
    polyscope::options::ssaaFactor = 3;

    auto tmp = tmp_dir();
    int i = 0;
    std::cout  << "Running Simulation in " << tmp.string() << std::endl;
    while(solver.elapsed_time < target_time){
        visualize_w(solver,geom);
        std::stringstream  ss;
        ss << "vid" << std::setw(4) << std::setfill('0') << i << ".png";
        polyscope::screenshot(tmp/ss.str());

        step(solver);
        i++;

    }


    fs::path video_dir = (std::filesystem::path(__FILE__).parent_path().parent_path().parent_path()) / "tex" / "thesis"/"figures"/"videos";;

    std::cout  << "Creating video " << tmp.string() << std::endl;
    std::stringstream  ss;
    ss << "cd \"" << tmp.string() << "\" &&";
    ss << "ffmpeg -framerate " << 1./solver.dt << " -i vid%04d.png -c:v libx264 -pix_fmt yuv420p -y "<< video_dir / (name +".mp4");
    std::cout <<"Executing: '"<< ss.str() <<"'" <<  std::endl;
    std::system(ss.str().c_str());

    std::cout << "Removing " << tmp << std::endl;
    std::filesystem::remove_all(tmp);
}

TEST(VideoTest,Grid){
    auto [mesh,geom] = readManifoldSurfaceMesh(std::filesystem::path(__FILE__).parent_path() / "models" / "grid_fine.stl");
    polyscope::init();

    AdaptiveFluidSolverData d;
    d.adaptive_space = false; d.adaptive_time = false; d.dt = 0.01;
    AdaptiveFluidSolver solver(*mesh,*geom,d);
    solver.wc.w = TaylorInitializer().wc(solver.tri.intrinsicTriangulation(),*geom).w;

    createVideo(solver,*geom,"grid",8,[](AdaptiveFluidSolver& s) {s.step(); std::fill(s.wc.c.begin(), s.wc.c.end(), 0);});
}

TEST(VideoTest,fixedHarmonicCoefficients){
    auto [mesh,geom] = readManifoldSurfaceMesh(std::filesystem::path(__FILE__).parent_path() / "models" / "cheese_oriented.stl");
    polyscope::init();

    AdaptiveFluidSolverData d;
    d.adaptive_space = false; d.adaptive_time = false; d.dt = 0.01;
    AdaptiveFluidSolver solver(*mesh,*geom,d);
    solver.wc.w = TaylorInitializer().wc(solver.tri.intrinsicTriangulation(),*geom).w;

    createVideo(solver,*geom,"constant_harm_coefficients",8,[](AdaptiveFluidSolver& s) {s.step(); std::fill(s.wc.c.begin(), s.wc.c.end(), 0);});
}

TEST(VideoTest,staticFC){
    auto [mesh,geom] = readManifoldSurfaceMesh(std::filesystem::path(__FILE__).parent_path() / "models" / "cheese_oriented.stl");
    polyscope::init();

    AdaptiveFluidSolverData d;
    d.adaptive_space = false; d.adaptive_time = false; d.dt = 0.01;
    AdaptiveFluidSolver solver(*mesh,*geom,d);
    solver.wc.w = TaylorInitializer().wc(solver.tri.intrinsicTriangulation(),*geom).w;

    createVideo(solver,*geom,"static",8,[](AdaptiveFluidSolver& s) {s.step();});
}
TEST(VideoTest,adaptiveFC){
    auto [mesh,geom] = readManifoldSurfaceMesh(std::filesystem::path(__FILE__).parent_path() / "models" / "cheese_min.stl");
    polyscope::init();

    AdaptiveFluidSolverData d;
    d.adaptive_space = true; d.adaptive_time = false; d.dt = 0.01;
    AdaptiveFluidSolver solver(*mesh,*geom,d);
    solver.wc.w = TaylorInitializer().wc(solver.tri.intrinsicTriangulation(),*geom).w;

    createVideo(solver,*geom,"adaptive",8,[](AdaptiveFluidSolver& s) {s.step();});
}
}