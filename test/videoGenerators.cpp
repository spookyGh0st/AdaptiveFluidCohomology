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
void display_progress(double v, int width = 50) {
    v = std::max(0.0, std::min(1.0, v));
    int filled = int(v * width);

    // Move cursor to start of line and clear the line
    std::printf("\33[2K\r");

    std::printf("[");
    for (int i = 0; i < filled; ++i) std::printf("=");
    for (int i = filled; i < width; ++i) std::printf(" ");
    std::printf("] %3d%%", int(std::round(v * 100)));

    std::fflush(stdout);

    if (v == 1.0) std::printf("\n");
}



void createVideo(AdaptiveFluidSolver& solver, VertexPositionGeometry& geom, const std::string& name, double target_time,std::function<void(AdaptiveFluidSolver&)> step){
    visualize_w(solver,geom);
    polyscope::view::ensureViewValid();
    polyscope::view::fov = 30;
    polyscope::view::projectionMode = polyscope::ProjectionMode::Perspective;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    polyscope::view::setWindowSize(1080,1080);
    polyscope::options::ssaaFactor = 3;
    polyscope::view::lookAt(glm::vec3{-5., 5., 2.}, glm::vec3{0., 0, 0});

    auto tmp = tmp_dir();
    int i = 0;
    std::cout  << "Running Simulation in " << tmp.string() << std::endl;
    fs::path thumbnail;
    while(solver.elapsed_time < target_time){
        display_progress(solver.elapsed_time/target_time);
        std::stringstream  ss;
        ss << "vid" << std::setw(16) << std::setfill('0') << i << ".png";
        polyscope::screenshot(tmp/ss.str());
        if (i == 0) thumbnail = tmp/ss.str();
        step(solver);
        i++;
        visualize_w(solver,geom);

    }


    fs::path video_dir = (std::filesystem::path(__FILE__).parent_path().parent_path().parent_path()) / "tex" / "thesis"/"figures"/"videos";;

    std::cout  << "Creating video " << tmp.string() << std::endl;
    std::stringstream  ss;
    ss << "cd \"" << tmp.string() << "\" &&";
    ss << "ffmpeg -framerate " << 1./solver.dt << " -i vid%16d.png -c:v libx264 -pix_fmt yuv420p -r 60 -y "<< video_dir / (name +".mp4");
    std::cout <<"Executing: '"<< ss.str() <<"'" <<  std::endl;
    std::system(ss.str().c_str());
    fs::copy_file(thumbnail, video_dir/ (name+".png"),fs::copy_options::overwrite_existing);

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
    for (int i = 0; i < 8; ++i) {
        solver.adapt();
        solver.wc.w = TaylorInitializer().wc(solver.tri.intrinsicTriangulation(),*geom).w;
    }

    createVideo(solver,*geom,"static",8,[](AdaptiveFluidSolver& s) {s.step();});
}
TEST(VideoTest,adaptiveFC){
    auto [mesh,geom] = readManifoldSurfaceMesh(std::filesystem::path(__FILE__).parent_path() / "models" / "cheese_min.stl");
    polyscope::init();

    AdaptiveFluidSolverData d;
    d.adaptive_space = true; d.adaptive_time = false; d.dt = 0.01;
    AdaptiveFluidSolver solver(*mesh,*geom,d);
    solver.wc.w = TaylorInitializer().wc(solver.tri.intrinsicTriangulation(),*geom).w;
    for (int i = 0; i < 8; ++i) {
        solver.adapt();
        solver.wc.w = TaylorInitializer().wc(solver.tri.intrinsicTriangulation(),*geom).w;
    }

    createVideo(solver,*geom,"adaptive",8,[](AdaptiveFluidSolver& s) {s.step();});
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

TEST(VideoTest,staticShearOnTorus) {

    auto [mesh,geom,param] = readParameterizedManifoldSurfaceMesh(std::filesystem::path(__FILE__).parent_path() / "models" /"torus.obj");
    VertexData<Vector2> uv(*mesh); for (Vertex v: mesh->vertices()) { uv[v] = (*param)[v.corner()]; }

    AdaptiveFluidSolverData data_comp_h(DOPRI5PresetConf::HIGH,DoerflerPresetConf::LOW,0.001,false, false,MARKING_STRATEGY::LONGEST_EDGE,false,false);
    data_comp_h.doerflerConf.threshold_refine = 2;
    data_comp_h.doerflerConf.threshold_coarse = 0.5;
    AdaptiveFluidSolver solver(*mesh,*geom, data_comp_h);

    for (int i = 0; i < 16; ++i) {
        for (Vertex v: solver.tri.mesh().vertices()){
            auto vec = solver.tri.intrinsicTriangulation().vertexLocations[v].interpolate(uv);
            solver.wc.w[v] = double_shear_layer(vec.x,vec.y);
        }
         solver.adapt();
    }
    for (Vertex v : solver.tri.mesh().vertices()) {
        auto vec = solver.tri.intrinsicTriangulation().vertexLocations[v].interpolate(uv);
        solver.wc.w[v] = double_shear_layer(vec.x, vec.y);
    }
    polyscope::init();
    createVideo(solver,*geom,"shear_static",3,[](AdaptiveFluidSolver& s) {s.step();});
}
TEST(VideoTest,adaptiveShearOnTorus) {

    auto [mesh,geom,param] = readParameterizedManifoldSurfaceMesh(std::filesystem::path(__FILE__).parent_path() / "models" /"torus.obj");
    VertexData<Vector2> uv(*mesh); for (Vertex v: mesh->vertices()) { uv[v] = (*param)[v.corner()]; }

    AdaptiveFluidSolverData data_comp_h(DOPRI5PresetConf::HIGH,DoerflerPresetConf::LOW,0.001,false,true,MARKING_STRATEGY::LONGEST_EDGE,false,false);
    data_comp_h.doerflerConf.threshold_refine = 2;
    data_comp_h.doerflerConf.threshold_coarse = 0.5;
    AdaptiveFluidSolver solver(*mesh,*geom, data_comp_h);

    for (int i = 0; i < 16; ++i) {
        for (Vertex v: solver.tri.mesh().vertices()){
            auto vec = solver.tri.intrinsicTriangulation().vertexLocations[v].interpolate(uv);
            solver.wc.w[v] = double_shear_layer(vec.x,vec.y);
        }
         solver.adapt();
    }
    for (Vertex v : solver.tri.mesh().vertices()) {
        auto vec = solver.tri.intrinsicTriangulation().vertexLocations[v].interpolate(uv);
        solver.wc.w[v] = double_shear_layer(vec.x, vec.y);
    }
    polyscope::init();
    createVideo(solver,*geom,"shear_adaptive",3,[](AdaptiveFluidSolver& s) {s.step();});
}


}