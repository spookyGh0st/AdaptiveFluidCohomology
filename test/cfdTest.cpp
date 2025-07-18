#include <eider/homology.h>
#include <eider/cfd.h>

#include <filesystem>
#include "geometrycentral/surface/meshio.h"
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <implot.h>

#include <gtest/gtest.h>

#include "eider/util.h"
#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include "geometrycentral/surface/surface_point.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;


wc_wrapper init_wc(SurfaceMesh& mesh, VertexPositionGeometry& geo, std::vector<FaceData<Vector2>> h) {
    wc_wrapper wc;
    wc.w = VertexData<double>(mesh,1);
    wc.c = std::vector<double>(h.size(), 0);
    wc.c[0] = 0;
    wc.c[1] = 0.5;
    return wc;
}

wc_wrapper init_taylor(SurfaceMesh& mesh, VertexPositionGeometry& geo, std::vector<FaceData<Vector2>> h, double vorticity_distance, Vector2 offset, Vector2 cuttof, Vector2 cuttof_offset) {
    VertexData<double> w(mesh,0);
    double k = 2*PI / vorticity_distance;
    double A = 1;
    for (Vertex v : mesh.vertices())
    {
        Vector3 p = geo.vertexPositions[v];
        double x = p.x  + offset.x, y = p.y + offset.y;
        if (abs(p.x+cuttof_offset.x) > cuttof.x || abs(p.y+cuttof_offset.y) > cuttof.y) { continue;}
        w[v] = 2 * A *k * cos(k*x) * cos(k*y);
    }
    wc_wrapper wc;
    wc.w = w;
    wc.c = std::vector<double>(h.size(), 0);
    return wc;
}

TEST(cfdTest, testSec15)
{
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" /"torus_bounded.obj";
    auto [m,g] = readManifoldSurfaceMesh(fds.string());
    std::vector<FaceData<Vector2>> h= orthonormal_hom_basis(*m,*g);

    StreamFunctionSolver S;
    S.compute(*m,*g);
    wc_wrapper wc = init_wc(*m, *g, h);
    g->requireHalfedgeVectorsInFace();
    velocity_wrapper vel = velocity(*m,*g,wc,h, S);

    g->requireFaceTangentBasis();
    FaceData<Vector3> e1(*m),e2(*m);
    for (Face f: m->faces()) { e1[f] = g->faceTangentBasis[f][0], e2[f] = g->faceTangentBasis[f][1]; }

    polyscope::init();
    ImPlot::CreateContext();
    polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh("M", g->vertexPositions,m->getFaceVertexList());
    pm->addVertexScalarQuantity("vorticity",wc.w)->setEnabled(true);
    pm->addFaceTangentVectorQuantity("velocity",vel.u,e1,e2)->setEnabled(true);
    std::size_t i = 0;
    for (const auto& b: h) {
        pm->addFaceTangentVectorQuantity("Hom basis " + std::to_string(i),b,e1,e2);
        i++;
    }

    float dt = 0.001;
    bool running = false, fix_c = false;
    std::vector<float> c1({float(wc.c[0])}), c2({float(wc.c[1])}), t({0});
    polyscope::state::userCallback = [&]() {
        if (ImGui::Button("reset")) {
            wc = init_wc(*m,*g,h);
            vel = velocity(*m,*g,wc,h, S);
            pm->addVertexScalarQuantity("vorticity",wc.w);
            pm->addFaceTangentVectorQuantity("velocity",vel.u,e1,e2);
        };
        ImGui::InputFloat("delta time",&dt,0.001,0.01);
        ImGui::Checkbox("Run", &running);
        if (running || ImGui::Button("Advance")) {
            wc = RK4Step(*m,*g,h,wc, dt, S);
            vel = velocity(*m,*g,wc,h, S);
            pm->addVertexScalarQuantity("vorticity",wc.w);
            pm->addFaceTangentVectorQuantity("velocity",vel.u,e1,e2);
            c1.push_back(wc.c[0]); c2.push_back(wc.c[1]);
            t.push_back(t.back() + dt);
        }
      if (ImPlot::BeginPlot("Homology coefficient")) {
        ImPlot::SetupAxis(ImAxis_X1, "Time",ImPlotAxisFlags_AutoFit);
        ImPlot::SetupAxis(ImAxis_Y1, "coefficient", ImPlotAxisFlags_AutoFit);
        ImPlot::PlotLine("c1", t.data(), c1.data(),int(t.size()));
        ImPlot::PlotLine("c2", t.data(), c2.data(),int(t.size()));
        ImPlot::EndPlot();
      }
    };// specify the callback
    polyscope::show();

}

inline double wedge(Vector2 a, Vector2 b){
  return a[0] * b[1] - a[1] * b[0];
}


TEST(cfdTest, IntrinsicHole)
{
    using namespace geometrycentral::surface;
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" /"grid_hole.stl";
    auto [parent_m,parent_g] = readManifoldSurfaceMesh(fds.string());
    IntegerCoordinatesIntrinsicTriangulation int_T (*parent_m,*parent_g);
    int_T.flipToDelaunay();
    int_T.delaunayRefine(20);
    int_T.refreshQuantities();

    ManifoldSurfaceMesh& m = *int_T.intrinsicMesh;
    m.compress();
    VertexData<Vector3> int_positions(m) ;
    for (Vertex v : m.vertices()) {
        int_positions[v] = int_T.vertexLocations[v].interpolate(parent_g->vertexPositions);
    }
    VertexPositionGeometry g(m, int_positions);

    g.requireFaceTangentBasis();
    g.requireHalfedgeVectorsInFace();
    auto h = orthonormal_hom_basis(m,g);

    FaceData<Vector3> e1(m),e2(m);
    for (Face f: m.faces()) { e1[f] = g.faceTangentBasis[f][0], e2[f] = g.faceTangentBasis[f][1]; }

    StreamFunctionSolver S;
    S.compute(m,g);

    // wc_wrapper wc = init_wc(*m, *g, h);
    float dt = 0.001, v_dist  = 0.5, x_add = v_dist/4, y_add = 0, x_cuttof = v_dist/2, y_cuttof = v_dist/4, x_cutt_offset = 0, y_cutt_offset = 0.5;
    wc_wrapper wc = init_taylor(m, g, h, v_dist, Vector2(x_add, y_add), Vector2(x_cuttof, y_cuttof), Vector2(x_cutt_offset, y_cutt_offset));
    velocity_wrapper vel = velocity(m,g,wc,h, S);

    std::vector<float> t_data ({0}); std::vector<std::vector<float>> c_data(h.size()), c_should(h.size());
    for (int i = 0; i < h.size(); ++i) {
      c_data[i].push_back(float(wc.c[i]));
      c_should[i].push_back(float(wc.c[i]));
    }

    polyscope::init();
    ImPlot::CreateContext();

    polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh("M", g.vertexPositions,m.getFaceVertexList(), polyscopePermutations(m));
    auto vispm = [&](){
      pm->addVertexScalarQuantity("vorticity",wc.w);
      pm->addFaceTangentVectorQuantity("velocity",vel.u,e1,e2);
      pm->addFaceScalarQuantity("error",vel.residual);
      pm->addVertexScalarQuantity("stream_function",vel.stream_function);
      VertexData<double> curlU(m,0);
      for (Vertex v: m.vertices()) { curlU[v] = curl(g,v,vel.u) ;}
      pm->addVertexScalarQuantity("vorticity new",curlU);
    };
    vispm();
    {
      std::size_t i = 0;
      for (const auto& b: h) {
        pm->addFaceTangentVectorQuantity("Hom basis " + std::to_string(i),b,e1,e2);
        i++;
      }
    }

    DOPRI5_conf conf=  DOPRI5_conf();

    bool running = false, fix_c = false;
    polyscope::state::userCallback = [&]() {
        ImGui::InputFloat("vorticity distance",&v_dist,0.125,0.5);
        ImGui::InputFloat("x_add",&x_add,0.125,0.5); ImGui::InputFloat("y_add",&y_add,0.125,0.5);
        ImGui::InputFloat("x_cuttof",&x_cuttof,0.125,0.5); ImGui::InputFloat("y_cuttoff",&y_cuttof,0.125,0.5);
        ImGui::InputFloat("x_cut offset",&x_cutt_offset,0.125,0.5); ImGui::InputFloat("y_cut offset",&y_cutt_offset,0.125,0.5);
        if (ImGui::Button("reset")) {
            wc = init_taylor(m, g, h, v_dist, Vector2(x_add,y_add), Vector2(x_cuttof, y_cuttof), Vector2(x_cutt_offset, y_cutt_offset));
            vel = velocity(m,g,wc,h, S);
            t_data = std::vector<float>({0});
            for (int j = 0; j < c_data.size(); ++j) {
              c_data[j].clear(); c_data[j].push_back(float(wc.c[j]));
              c_should[j].clear(); c_should[j].push_back(float(wc.c[j]));
            }
            vispm();
        };
        ImGui::InputDouble("Absolut Error",&conf.Atol_i,0,0,"%.10f");
        ImGui::InputDouble("relative Error",&conf.Rtol_i,0,0,"%.10f");
        ImGui::InputDouble("fac Min",&conf.facmin);
        ImGui::InputDouble("fac max",&conf.faxmax);
        ImGui::InputFloat("delta time",&dt,0.001,0.01);
        ImGui::Checkbox("Run", &running);
        ImGui::Checkbox("Fix c", &fix_c);
        if (running || ImGui::Button("Advance")) {
            auto tmpc = wc.c; double tmpdt;
            auto dp5_sample = adaptive_step(m,g,h,wc, dt, S,conf);
            wc = dp5_sample.wc;
            if (fix_c) wc.c = tmpc;
            dt = float(dp5_sample.t_future);
            vel = velocity(m,g,wc,h, S);
            vispm();

            // update plots arrays
            t_data.push_back(t_data.back() + float(dp5_sample.t_past));
            for (int j = 0; j < wc.c.size(); ++j) {
              c_data[j].push_back(float(wc.c[j]));

              double c_sum = 0;
              for (Face f: m.faces())
                c_sum += dot(h[j][f],vel.u[f]) * g.faceAreas[f];
              c_should[j].push_back(float(c_sum));
            }

        }
      if (ImPlot::BeginPlot("Homology coefficient")) {
        ImPlot::SetupAxis(ImAxis_X1, "Time",ImPlotAxisFlags_AutoFit);
        ImPlot::SetupAxis(ImAxis_Y1, "coefficient", ImPlotAxisFlags_AutoFit);
        for (int j = 0; j < c_data.size() ; ++j) {
          ImPlot::PlotLine(("c" + std::to_string(j)).c_str(), t_data.data(), c_data[j].data(),int(t_data.size()));
          ImPlot::PlotLine(("c_precise" + std::to_string(j)).c_str(), t_data.data(), c_should[j].data(),int(t_data.size()));
        }
        ImPlot::EndPlot();
      }
        for (int i = 0; i< wc.c.size(); i++) {
            ImGui::Text("c%d: %f",i,wc.c[i]);
        }

    };// specify the callback
    polyscope::show();
}

wc_wrapper init_torus_taylor(SurfaceMesh& m, IntrinsicGeometryInterface& g, CornerData<Vector2>& p, std::vector<FaceData<Vector2>> h, double vorticity_distance, double A) {
    wc_wrapper wc; wc.w = VertexData<double>(m);
    double k = 2 * PI / vorticity_distance;
    for (Vertex v : m.vertices())
    {
        double x = (p)[v.corner()].x, y = (p)[v.corner()].y;
        wc.w[v] = 2 * A * k * cos(k * x) * cos(k * y);
    }
    wc.c.resize(h.size(), 0);
    for (int i = 0; i < h.size(); i++)
    {
        for (Face f : m.faces())
        {
            double x = 0, y= 0;
            for (Corner c: f.adjacentCorners()) { x += (p)[c].x * (1./3), y = (p)[c].y * (1./3); }
            Vector2 u = h[0][f], v = Vector2(A * cos(k * x) * sin(k * y), A * sin(k * x) * cos(k * y));
            wc.c[i] += (u.x * v.y - u.y * v.x) * g.faceAreas[f];
        }
    }
    wc.c[0] = 1;
    return wc;
}

TEST(cfdTest, testTorus)
{
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" /"torus.obj";
    std::unique_ptr<ManifoldSurfaceMesh> m; std::unique_ptr<VertexPositionGeometry> g;
    std::unique_ptr<CornerData<Vector2>> p;
    std::tie(m,g,p) = readParameterizedManifoldSurfaceMesh(fds.string());
    std::vector<FaceData<Vector2>> h= orthonormal_hom_basis(*m,*g);
    StreamFunctionSolver S;
    S.compute(*m,*g);

    float v_dist = 1, A = 1;
    wc_wrapper wc = init_torus_taylor(*m, *g, *p,h, v_dist,A);

    g->requireHalfedgeVectorsInFace();
    velocity_wrapper vel = velocity(*m,*g,wc,h, S);

    g->requireFaceTangentBasis();
    FaceData<Vector3> e1(*m),e2(*m);
    for (Face f: m->faces()) { e1[f] = g->faceTangentBasis[f][0], e2[f] = g->faceTangentBasis[f][1]; }

    polyscope::init();
    polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh("M", g->vertexPositions,m->getFaceVertexList());
    pm->addVertexScalarQuantity("vorticity",wc.w)->setEnabled(true);
    pm->addFaceTangentVectorQuantity("velocity",vel.u,e1,e2)->setEnabled(true);
    std::size_t i = 0;
    for (const auto& b: h) {
        pm->addFaceTangentVectorQuantity("Hom basis " + std::to_string(i),b,e1,e2);
        i++;
    }

    bool running = false;
    float dt = 0.001;
    polyscope::state::userCallback = [&]() {
        ImGui::InputFloat("vorticity distance",&v_dist,0.125,0.5);
        ImGui::SliderFloat("A",&A,0,10);
        if (ImGui::Button("reset")) {
            wc = init_torus_taylor(*m, *g, *p,h, v_dist,A);
            vel = velocity(*m,*g,wc,h, S);
            pm->addVertexScalarQuantity("vorticity",wc.w);
            pm->addFaceTangentVectorQuantity("velocity",vel.u,e1,e2);
        };
        ImGui::InputFloat("delta time",&dt,0.01,0.1);
        ImGui::Checkbox("Run", &running);
        if (running || ImGui::Button("Advance")) {
            wc = RK4Step(*m,*g,h,wc, dt, S);
            vel = velocity(*m,*g,wc,h, S);
            pm->addVertexScalarQuantity("vorticity",wc.w);
            pm->addFaceTangentVectorQuantity("velocity",vel.u,e1,e2);
        }
        for (int i = 0; i< wc.c.size(); i++) {
            ImGui::Text("c%d: %f",i,wc.c[i]);
            if (i < wc.c.size() - 1) {}
        }

    };// specify the callback
    polyscope::show();

}

