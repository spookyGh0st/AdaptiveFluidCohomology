#include "TaylorInitializer.h"
#include "eider/util.h"
#include "eider/Stopwatch.h"
#include <algorithm>
#include <eider/AdaptiveFluidSolver.h>
#include <filesystem>
#include <gtest/gtest.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include "testutil.h"
#include <thread>
#include <vector>



namespace geometrycentral::surface{
namespace fs = std::filesystem;
fs::path tmp_dir(){
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&t), "%Y-%m-%dT%H-%M-%S-");

    std::random_device dev;
    std::mt19937 prng(dev());
    std::uniform_int_distribution<uint64_t> rand(0);
    ss << std::hex << rand(prng);

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

void visualize_c_h0(AdaptiveFluidSolver& solver, VertexPositionGeometry& geom){
    AdaptiveTriangulation& atri = solver.tri;
    auto& m = atri.mesh();
    auto vp = intrinsic_vertex_positions(atri.intrinsicTriangulation(),geom);
    auto g_vis =VertexPositionGeometry(m,vp);

    polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh("mesh", vp,atri.mesh().getFaceVertexList(), polyscopePermutations(atri.mesh()));
    pm->setEdgeColor(glm::vec3(0,0,0));
    pm->setEdgeWidth(0.4);

    HalfedgeData<int> nextleft(m,0);
    for(Halfedge h: solver.hom.homologyB[0] ){
        auto nl = solver.hom.homologyB[0].nextLeft[h];
        if(nl.has_value() && nl.value()) nextleft[h] = 1;
        if(nl.has_value() && !nl.value()) nextleft[h] = -1;
    }

    auto* vsq = pm->addHalfedgeScalarQuantity("cycle",nextleft);
    vsq->setColorMap("coolwarm");
    vsq->setEnabled(true);

    g_vis.requireFaceTangentBasis();
    FaceData<Vector3> e1(m),e2(m);
    for (Face f: m.faces()) { e1[f] = g_vis.faceTangentBasis[f][0], e2[f] = g_vis.faceTangentBasis[f][1];  }

    std::string name = "h0";
    auto* ftq = pm->addFaceTangentVectorQuantity(name,solver.h[0], e1,e2);
    ftq->setVectorLengthRange(2);
    ftq->setVectorLengthScale(0.02);
    ftq->setEnabled(true);
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

void init_polyscope_2d(AdaptiveFluidSolver& solver, VertexPositionGeometry& geom) {
    if(!polyscope::isInitialized()) polyscope::init();
    visualize_w(solver,geom);
    polyscope::view::ensureViewValid();
    polyscope::view::fov = 30;
    polyscope::view::projectionMode = polyscope::ProjectionMode::Orthographic;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    polyscope::view::setWindowSize(1080,1080);
    polyscope::options::ssaaFactor = 3;
    // polyscope::view::lookAt(glm::vec3{-5., 5., 2.}, glm::vec3{0., 0, 0});
}
void init_polyscope_2d_zoom(AdaptiveFluidSolver& solver, VertexPositionGeometry& geom) {
    if(!polyscope::isInitialized()) polyscope::init();
    visualize_w(solver,geom);
    polyscope::view::ensureViewValid();
    polyscope::view::fov = 10;
    polyscope::view::projectionMode = polyscope::ProjectionMode::Orthographic;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    polyscope::view::setWindowSize(1080,1080);
    polyscope::options::ssaaFactor = 3;
    // polyscope::view::lookAt(glm::vec3{-5., 5., 2.}, glm::vec3{0., 0, 0});
}

void init_polyscope_3d(AdaptiveFluidSolver& solver, VertexPositionGeometry& geom) {
    if(!polyscope::isInitialized()) polyscope::init();
    visualize_w(solver,geom);
    polyscope::view::ensureViewValid();
    polyscope::view::fov = 30;
    polyscope::view::projectionMode = polyscope::ProjectionMode::Perspective;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    polyscope::view::setWindowSize(1080,1080);
    polyscope::options::ssaaFactor = 3;
    polyscope::view::lookAt(glm::vec3{-5., 5., 2.}, glm::vec3{0., 0, 0});
}



void createVideo(
    AdaptiveFluidSolver& solver, VertexPositionGeometry& geom,
    const std::string& name, double target_time,
    std::function<void(AdaptiveFluidSolver&)> step,
    std::function<void(AdaptiveFluidSolver&, VertexPositionGeometry&)> visualize,
    double speed = 1
){
    polyscope::removeAllStructures();

    auto tmp = tmp_dir();
    int i = 0;
    fs::path manifest_path= tmp / "manifest.txt";
    std::ofstream manifest(manifest_path);

    std::cout  << "Running Simulation in " << tmp.string() << std::endl;
    fs::path thumbnail, last_frame;
    double elapsed_time = 0;
    while(solver.elapsed_time < target_time){
        display_progress(solver.elapsed_time/target_time);
        std::stringstream  ss; ss << "vid" << std::setw(16) << std::setfill('0') << i << ".png";

        visualize(solver,geom);
        last_frame = tmp/ss.str();
        polyscope::screenshot(last_frame);
        manifest << "file '" << last_frame.string() << "'\n"  << "duration " << solver.elapsed_time/speed - elapsed_time/speed << "\n";

        if (i == 0) thumbnail = last_frame;

        elapsed_time = solver.elapsed_time; step(solver); i++;
    }
    manifest << "file '" << last_frame.string() << "'\n";
    manifest.close();


    fs::path video_dir = (std::filesystem::path(__FILE__).parent_path().parent_path().parent_path()) / "tex" / "thesis"/"figures"/"videos";;

    std::cout  << "Creating video " << tmp.string() << std::endl;
    std::stringstream  ss;
    ss << "cd \"" << tmp.string() << "\" &&";
    ss << "ffmpeg"
       << " -f concat -safe 0"              // Use the concat demuxer
       << " -i " << manifest_path.string()  // Input is the manifest file
       << " -vsync vfr"                     // Ensure variable frame rate is respected
       << " -c:v libx264"                   // Standard video codec
       << " -pix_fmt yuv420p"               // Pixel format for wide compatibility
       << " -r 60"                           // output with 60 fps
       << " -y "                            // Overwrite output file if it exists
       << video_dir / (name + ".mp4");
    std::cout <<"Executing: '"<< ss.str() <<"'" <<  std::endl;

    if (int return_code = std::system(ss.str().c_str()); return_code != 0) {
        throw std::runtime_error("External command failed with return code " + std::to_string(return_code) + ". Command was: '" + ss.str() + "'");
    }
    fs::copy_file(thumbnail, video_dir/ (name+".png"),fs::copy_options::overwrite_existing);

    std::cout << "Removing " << tmp << std::endl;
    std::filesystem::remove_all(tmp);
}

TEST(VideoTest,Grid){
    auto [mesh,geom] = readManifoldSurfaceMesh(std::filesystem::path(__FILE__).parent_path() / "models" / "grid_fine.stl");

    AdaptiveFluidSolverData d;
    d.adaptive_space = false; d.adaptive_time = false; d.dt = 0.01;
    AdaptiveFluidSolver solver(*mesh,*geom,d);
    solver.wc.w = TaylorInitializer().wc(solver.tri.intrinsicTriangulation(),*geom).w;

    init_polyscope_2d(solver,*geom);
    createVideo(solver,*geom,"grid",8,[](AdaptiveFluidSolver& s) {s.step(); std::fill(s.wc.c.begin(), s.wc.c.end(), 0);}, visualize_w);
}

TEST(VideoTest,fixedHarmonicCoefficients){
    auto [mesh,geom] = readManifoldSurfaceMesh(std::filesystem::path(__FILE__).parent_path() / "models" / "cheese_oriented.stl");

    AdaptiveFluidSolverData d;
    d.adaptive_space = false; d.adaptive_time = false; d.dt = 0.01;
    AdaptiveFluidSolver solver(*mesh,*geom,d);
    solver.wc.w = TaylorInitializer().wc(solver.tri.intrinsicTriangulation(),*geom).w;

    init_polyscope_2d(solver,*geom);
    createVideo(solver,*geom,"constant_harm_coefficients",8,[](AdaptiveFluidSolver& s) {s.step(); std::fill(s.wc.c.begin(), s.wc.c.end(), 0);}, visualize_w);
}

TEST(VideoTest,staticFC){
    auto [mesh,geom] = readManifoldSurfaceMesh(std::filesystem::path(__FILE__).parent_path() / "models" / "cheese_min.stl");

    AdaptiveFluidSolverData d;
    d.adaptive_space = false; d.adaptive_time = false; d.dt = 0.01;
    AdaptiveFluidSolver solver(*mesh,*geom,d);
    solver.wc.w = TaylorInitializer().wc(solver.tri.intrinsicTriangulation(),*geom).w;
    for (int i = 0; i < 8; ++i) {
        solver.adapt();
        solver.wc.w = TaylorInitializer().wc(solver.tri.intrinsicTriangulation(),*geom).w;
    }

    init_polyscope_2d(solver,*geom);
    createVideo(solver,*geom,"static",8,[](AdaptiveFluidSolver& s) {s.step();}, visualize_w);
}
TEST(VideoTest,adaptiveFC){
    auto [mesh,geom] = readManifoldSurfaceMesh(std::filesystem::path(__FILE__).parent_path() / "models" / "cheese_min.stl");

    AdaptiveFluidSolverData d;
    d.adaptive_space = true; d.adaptive_time = true; d.dt = 0.01;
    AdaptiveFluidSolver solver(*mesh,*geom,d);
    solver.wc.w = TaylorInitializer().wc(solver.tri.intrinsicTriangulation(),*geom).w;
    for (int i = 0; i < 8; ++i) {
        solver.adapt();
        solver.wc.w = TaylorInitializer().wc(solver.tri.intrinsicTriangulation(),*geom).w;
    }

    init_polyscope_2d(solver,*geom);
    createVideo(solver,*geom,"adaptive",8,[](AdaptiveFluidSolver& s) {s.step();}, visualize_w);
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

TEST(VideoTest,ShearOnTorusStatic) {

    auto [mesh,geom,param] = readParameterizedManifoldSurfaceMesh(std::filesystem::path(__FILE__).parent_path() / "models" /"torus.obj");
    VertexData<Vector2> uv(*mesh); for (Vertex v: mesh->vertices()) { uv[v] = (*param)[v.corner()]; }

    AdaptiveFluidSolverData data_comp_h(DOPRI5PresetConf::HIGH,DoerflerPresetConf::LOW,0.0001,false, false,MARKING_STRATEGY::LONGEST_EDGE,false,false);
    data_comp_h.doerflerConf.threshold_refine = 1;
    data_comp_h.doerflerConf.threshold_coarse = 0.2;
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
    init_polyscope_3d(solver,*geom);
    createVideo(solver,*geom,"shear_static",3,[](AdaptiveFluidSolver& s) {s.step();},visualize_w,0.2);
}
TEST(VideoTest,ShearOnTorusAdaptive) {

    auto [mesh,geom,param] = readParameterizedManifoldSurfaceMesh(std::filesystem::path(__FILE__).parent_path() / "models" /"torus.obj");
    VertexData<Vector2> uv(*mesh); for (Vertex v: mesh->vertices()) { uv[v] = (*param)[v.corner()]; }

    AdaptiveFluidSolverData data_comp_h(DOPRI5PresetConf::LOW,DoerflerPresetConf::LOW,0.001,false,true,MARKING_STRATEGY::LONGEST_EDGE,false,false);
    data_comp_h.doerflerConf.threshold_refine = 1;
    data_comp_h.doerflerConf.threshold_coarse = 0.2;
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
    init_polyscope_3d(solver,*geom);
    createVideo(solver,*geom,"shear_adaptive",3,[](AdaptiveFluidSolver& s) {s.step();},visualize_w, 0.2);
}

void performAdaptiveRefinement(AdaptiveFluidSolver& solver, ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, int iterations) {
    solver.wc.w = TaylorInitializer().wc(solver.tri.intrinsicTriangulation(), geom).w;
    for (int i = 0; i < iterations; ++i) {
        solver.adapt();
        solver.wc.w = TaylorInitializer().wc(solver.tri.intrinsicTriangulation(), geom).w;
    }
}

TEST(VideoTest, EvaluationS) {
    auto [mesh, geom] = readManifoldSurfaceMesh(std::filesystem::path(__FILE__).parent_path() / "models" / "cheese_min.stl");
    auto [meshO, geomO] = uniform_refine(*mesh, *geom, 4);

    AdaptiveFluidSolver solver_static(*meshO, *geomO, static_solver_data);
    solver_static.wc.w = TaylorInitializer().wc(solver_static.tri.intrinsicTriangulation(), *geomO).w;

    init_polyscope_2d(solver_static, *geomO);

    createVideo(solver_static, *geomO, "evaluation_s", 6, [](AdaptiveFluidSolver& s) {s.step();}, visualize_w);
}
TEST(VideoTest, EvaluationC) {
    auto [mesh, geom] = readManifoldSurfaceMesh(std::filesystem::path(__FILE__).parent_path() / "models" / "cheese_min.stl");

    AdaptiveFluidSolverData data_comp_h(DOPRI5PresetConf::LOW, DoerflerPresetConf::LOW, 0.01, true, true, MARKING_STRATEGY::PATTERN, false, false);
    AdaptiveFluidSolver solver_adapt_c(*mesh, *geom, data_comp_h);

    performAdaptiveRefinement(solver_adapt_c, *mesh, *geom, 8);

    init_polyscope_2d(solver_adapt_c, *geom);
    createVideo(solver_adapt_c, *geom, "evaluation_c", 6, [](AdaptiveFluidSolver& s) {s.step();}, visualize_w);
}
TEST(VideoTest, EvaluationI) {
    auto [mesh, geom] = readManifoldSurfaceMesh(std::filesystem::path(__FILE__).parent_path() / "models" / "cheese_min.stl");

    AdaptiveFluidSolverData data_comp_h(DOPRI5PresetConf::LOW, DoerflerPresetConf::LOW, 0.01, true, true, MARKING_STRATEGY::PATTERN, false, false);
    AdaptiveFluidSolverData data_interp_ha = data_comp_h;
    data_interp_ha.interpolate_harmonic_basis = true;
    data_interp_ha.use_interpolated_harmonic_basis = true;

    AdaptiveFluidSolver solver_adapt_i(*mesh, *geom, data_interp_ha);

    performAdaptiveRefinement(solver_adapt_i, *mesh, *geom, 8);

    init_polyscope_2d(solver_adapt_i, *geom);
    createVideo(solver_adapt_i, *geom, "evaluation_i", 6, [](AdaptiveFluidSolver& s) {s.step();}, visualize_w);
}

TEST(VideoTest, EvaluationCH) {
    auto [mesh, geom] = readManifoldSurfaceMesh(std::filesystem::path(__FILE__).parent_path() / "models" / "cheese_min.stl");

    AdaptiveFluidSolverData data_comp_h(DOPRI5PresetConf::LOW, DoerflerPresetConf::LOW, 0.01, true, true, MARKING_STRATEGY::PATTERN, false, false);
    AdaptiveFluidSolver solver_adapt_c(*mesh, *geom, data_comp_h);

    performAdaptiveRefinement(solver_adapt_c, *mesh, *geom, 8);

    init_polyscope_2d_zoom(solver_adapt_c, *geom);
    createVideo(solver_adapt_c, *geom, "evaluation_c-h", 6, [](AdaptiveFluidSolver& s) {s.step();}, visualize_c_h0);
}

std::pair<MeshP, GeomP> irregular_geometry(){
    auto [mesh, geom] = readManifoldSurfaceMesh(std::filesystem::path(__FILE__).parent_path() / "models" / "cheese_min.stl");
    return uniform_refine(*mesh,*geom,2,MARKING_STRATEGY::RANDOM);
}

TEST(VideoTest, Irregular_S){
    auto[mesh,geom] = irregular_geometry();
    AdaptiveFluidSolverData staticD(DOPRI5PresetConf::LOW,DoerflerPresetConf::UNIFORM_REFINE,0.001,false,false,MARKING_STRATEGY::RANDOM,false);
    AdaptiveFluidSolver solver(*mesh,*geom,staticD);
    performAdaptiveRefinement(solver, *mesh, *geom, 0);
    init_polyscope_2d(solver, *geom);
    createVideo(solver, *geom, "cheese_irregular_s", 3, [](AdaptiveFluidSolver& s) {s.step();}, visualize_w);
}

TEST(VideoTest, Irregular_C){
    auto[mesh,geom] = irregular_geometry();
    AdaptiveFluidSolverData data(DOPRI5PresetConf::LOW,DoerflerPresetConf::LOW,0.01,true, true,MARKING_STRATEGY::LONGEST_EDGE,false);
    AdaptiveFluidSolver solver(*mesh,*geom,data);
    performAdaptiveRefinement(solver, *mesh, *geom, 8);
    init_polyscope_2d(solver, *geom);
    createVideo(solver, *geom, "cheese_irregular-c", 3, [](AdaptiveFluidSolver& s) {s.step();}, visualize_w);
}

TEST(VideoTest, Irregular_I){
    auto[mesh,geom] = irregular_geometry();
    AdaptiveFluidSolverData data(DOPRI5PresetConf::LOW,DoerflerPresetConf::LOW,0.01,true, true,MARKING_STRATEGY::LONGEST_EDGE,true, true);
    AdaptiveFluidSolver solver(*mesh,*geom,data);
    performAdaptiveRefinement(solver, *mesh, *geom, 8);
    init_polyscope_2d(solver, *geom);
    createVideo(solver, *geom, "cheese_irregular-i", 3, [](AdaptiveFluidSolver& s) {s.step();}, visualize_w);
}

TEST(VideoTest, Irregular_CH) {
    auto[mesh,geom] = irregular_geometry();
    AdaptiveFluidSolverData data(DOPRI5PresetConf::LOW,DoerflerPresetConf::LOW,0.01,true, true,MARKING_STRATEGY::LONGEST_EDGE,false);
    AdaptiveFluidSolver solver(*mesh,*geom,data);
    performAdaptiveRefinement(solver, *mesh, *geom, 8);
    init_polyscope_2d_zoom(solver, *geom);
    createVideo(solver, *geom, "cheese_irregular-c-h", 3, [](AdaptiveFluidSolver& s) {s.step();}, visualize_c_h0);
}


}