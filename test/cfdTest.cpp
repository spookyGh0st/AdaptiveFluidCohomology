#include <eider/homology.h>
#include <eider/cfd.h>

#include <filesystem>
#include "geometrycentral/surface/meshio.h"
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include <gtest/gtest.h>

#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"


using namespace geometrycentral;
using namespace geometrycentral::surface;

double curl(Vertex v, IntrinsicGeometryInterface& geom, FaceData<Vector2> u){
    double curl_sum = 0.0f;
    double area_sum = 0.0f;
    for (Face f : v.adjacentFaces()) {
        double Af = geom.faceAreas[f];
        Vector2 vf = u[f];

        double face_curl = 0;
        for (Halfedge he : f.adjacentHalfedges())
        {
            Vector2 efp = geom.halfedgeVectorsInFace[he.next()];
            face_curl += dot(vf,efp.rotate90());
        }
        face_curl = face_curl / (2.0f * Af);
        curl_sum += face_curl * Af;
        area_sum += Af;
    }
    return curl_sum / area_sum;
}
VertexData<double> curl(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom, FaceData<Vector2> u){
    VertexData<double> c(mesh);
    geom.requireHalfedgeVectorsInFace();
    for (Vertex v : mesh.vertices()) {
        c[v] = curl(v, geom, u);
    }
    geom.unrequireHalfedgeVectorsInFace();
    return c;
}

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
    fds = fds.parent_path()/ "models" /"torus_bounded_max.stl";
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
        }
        for (int i = 0; i< wc.c.size(); i++) {
            ImGui::Text("c%d: %f",i,wc.c[i]);
        }

    };// specify the callback
    polyscope::show();

}


TEST(cfdTest, testVorticesThroughHole)
{
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" /"punctured_disk.stl";
    auto [m,g] = readManifoldSurfaceMesh(fds.string());
    std::vector<FaceData<Vector2>> h= orthonormal_hom_basis(*m,*g);

    StreamFunctionSolver S;
    S.compute(*m,*g);

    // wc_wrapper wc = init_wc(*m, *g, h);
    float dt = 0.001, v_dist  = 0.5, x_add = v_dist/4, y_add = 0, x_cuttof = v_dist/2, y_cuttof = v_dist/4, x_cutt_offset = 0, y_cutt_offset = 0.5;
    wc_wrapper wc = init_taylor(*m, *g, h, v_dist, Vector2(x_add, y_add), Vector2(x_cuttof, y_cuttof), Vector2(x_cutt_offset, y_cutt_offset));
    g->requireHalfedgeVectorsInFace();
    velocity_wrapper vel = velocity(*m,*g,wc,h, S);

    g->requireFaceTangentBasis();
    FaceData<Vector3> e1(*m),e2(*m);
    for (Face f: m->faces()) { e1[f] = g->faceTangentBasis[f][0], e2[f] = g->faceTangentBasis[f][1]; }

    polyscope::init();
    polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh("M", g->vertexPositions,m->getFaceVertexList());
    pm->addVertexScalarQuantity("vorticity",wc.w)->setEnabled(true);
    pm->addFaceTangentVectorQuantity("velocity",vel.u,e1,e2)->setEnabled(true);
    pm->addFaceScalarQuantity("residual",vel.residual);
    std::size_t i = 0;
    for (const auto& b: h) {
        pm->addFaceTangentVectorQuantity("Hom basis " + std::to_string(i),b,e1,e2);
        i++;
    }

    bool running = false, fix_c = false;
    polyscope::state::userCallback = [&]() {
        ImGui::InputFloat("vorticity distance",&v_dist,0.125,0.5);
        ImGui::InputFloat("x_add",&x_add,0.125,0.5); ImGui::InputFloat("y_add",&y_add,0.125,0.5);
        ImGui::InputFloat("x_cuttof",&x_cuttof,0.125,0.5); ImGui::InputFloat("y_cuttoff",&y_cuttof,0.125,0.5);
        ImGui::InputFloat("x_cut offset",&x_cutt_offset,0.125,0.5); ImGui::InputFloat("y_cut offset",&y_cutt_offset,0.125,0.5);
        if (ImGui::Button("reset")) {
            wc = init_taylor(*m, *g, h, v_dist, Vector2(x_add,y_add), Vector2(x_cuttof, y_cuttof), Vector2(x_cutt_offset, y_cutt_offset));
            vel = velocity(*m,*g,wc,h, S);
            pm->addVertexScalarQuantity("vorticity",wc.w);
            pm->addFaceTangentVectorQuantity("velocity",vel.u,e1,e2);
            pm->addFaceScalarQuantity("residual",vel.residual);
        };
        ImGui::InputFloat("delta time",&dt,0.001,0.01);
        ImGui::Checkbox("Run", &running);
        ImGui::Checkbox("Fix c", &fix_c);
        if (running || ImGui::Button("Advance")) {
            auto tmpc = wc.c;
            wc = RK4Step(*m,*g,h,wc, dt, S);
            if (fix_c) wc.c = tmpc;
            vel = velocity(*m,*g,wc,h, S);
            pm->addVertexScalarQuantity("vorticity",wc.w);
            pm->addFaceTangentVectorQuantity("velocity",vel.u,e1,e2);
            pm->addFaceScalarQuantity("residual",vel.residual);
        }
        for (int i = 0; i< wc.c.size(); i++) {
            ImGui::Text("c%d: %f",i,wc.c[i]);
        }

    };// specify the callback
    polyscope::show();

}


TEST(cfdTest, IntrinsicHole)
{
    using namespace geometrycentral::surface;
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" /"grid_hole.stl";
    auto [parent_m,parent_g] = readManifoldSurfaceMesh(fds.string());
    IntegerCoordinatesIntrinsicTriangulation int_T (*parent_m,*parent_g);
    int_T.flipToDelaunay();
    int_T.delaunayRefine(10);
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

    polyscope::init();
    polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh("M", g.vertexPositions,m.getFaceVertexList(), polyscopePermutations(m));
    pm->addVertexScalarQuantity("vorticity",wc.w)->setEnabled(true);
    pm->addFaceTangentVectorQuantity("velocity",vel.u,e1,e2)->setEnabled(true);
    std::size_t i = 0;
    for (const auto& b: h) {
        pm->addFaceTangentVectorQuantity("Hom basis " + std::to_string(i),b,e1,e2);
        i++;
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
            pm->addVertexScalarQuantity("vorticity",wc.w);
            pm->addFaceTangentVectorQuantity("velocity",vel.u,e1,e2);
        };
        ImGui::InputDouble("Absolut Error",&conf.Atol_i,0,0,"%.10f");
        ImGui::InputDouble("relative Error",&conf.Rtol_i,0,0,"%.10f");
        ImGui::InputDouble("fac Min",&conf.facmin);
        ImGui::InputDouble("fac max",&conf.faxmax);
        ImGui::InputFloat("delta time",&dt,0.001,0.01);
        ImGui::Checkbox("Run", &running);
        ImGui::Checkbox("Fix c", &fix_c);
        if (running || ImGui::Button("Advance")) {
            auto tmpc = wc.c;
            std::tie(wc,dt) = adaptive_step(m,g,h,wc, dt, S,conf);
            if (fix_c) wc.c = tmpc;
            vel = velocity(m,g,wc,h, S);
            pm->addVertexScalarQuantity("vorticity",wc.w);
            pm->addFaceTangentVectorQuantity("velocity",vel.u,e1,e2);
            pm->addFaceScalarQuantity("error",vel.residual);
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

