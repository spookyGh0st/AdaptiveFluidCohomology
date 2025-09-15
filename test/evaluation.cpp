#include "Stopwatch.h"
#include <eider/AdaptiveFluidSolver.h>
#include "Evaluator.h"
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

void run(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom){
    Evaluator ev;
    IntegerCoordinatesIntrinsicTriangulation ic(mesh,geom);
    AdaptiveTriangulation atri(ic);
    DOPRI5_conf dopri5Conf;
    DoeflerConf doeflerConf;
    wc_wrapper wc;
    AdaptiveFluidSolver solver(atri,wc,dopri5Conf,doeflerConf);
}

TEST(EvaluatorTest, Evaluate)
{
    results_folder();
    fs::path fsrc = fs::path(__FILE__).parent_path();
    fs::path fds = fsrc / "models" / "cheese_min.stl";
    fs::path fev = run_folder(results_folder());
    fs::path flatest = fev.parent_path() / "run_latest";

    auto [mesh,geom] = readManifoldSurfaceMesh(fds.string());
    polyscope::init();

    // Register your mesh with Polyscope
    auto* pm1 = polyscope::registerSurfaceMesh("cheese_min", geom->vertexPositions,mesh->getFaceVertexList(), polyscopePermutations(*mesh));
    auto* pm2 = polyscope::registerSurfaceMesh("cheese_max", geom->vertexPositions,mesh->getFaceVertexList(), polyscopePermutations(*mesh));
    pm1->setPosition(glm::vec3(-1.1,0,0));
    pm2->setPosition(glm::vec3( 1.1,0,0));
    // Set some options and save a named screenshot

    polyscope::view::ensureViewValid();
    polyscope::view::fov = 30;
    polyscope::view::projectionMode = polyscope::ProjectionMode::Orthographic;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    polyscope::view::setWindowSize(2160/2,1080/2);

    polyscope::screenshot(fev/"im.png",true);
    to_csv(DoeflerConf(),fev/"doerfler.csv");
    to_csv(DOPRI5_conf(),fev/"dopri.csv");
    copyFolder(fev,flatest);
    ImPlot::CreateContext();
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
