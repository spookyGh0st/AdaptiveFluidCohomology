#include <eider/homology.h>
#include <eider/cfd.h>

#include <filesystem>
#include "geometrycentral/surface/meshio.h"
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include <gtest/gtest.h>

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


TEST(cfdTest, testFluidSim)
{
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" /"torus_bounded_max.stl";
    auto [m,g] = readManifoldSurfaceMesh(fds.string());
    std::vector<FaceData<Vector2>> h= orthonormal_hom_basis(*m,*g);

    StreamFunctionSolver S;
    S.compute(*m,*g);

    wc_wrapper wc = init_wc(*m, *g, h);
    FaceData<Vector2> u = velocity(*m,*g,wc,h, S);

    g->requireFaceTangentBasis();
    FaceData<Vector3> e1(*m),e2(*m);
    for (Face f: m->faces()) { e1[f] = g->faceTangentBasis[f][0], e2[f] = g->faceTangentBasis[f][1]; }

    polyscope::init();
    polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh("M", g->vertexPositions,m->getFaceVertexList());
    pm->addVertexScalarQuantity("vorticity",wc.w)->setEnabled(true);
    pm->addFaceTangentVectorQuantity("velocity",u,e1,e2)->setEnabled(true);
    std::size_t i = 0;
    for (const auto& b: h) {
        pm->addFaceTangentVectorQuantity("Hom basis " + std::to_string(i),b,e1,e2);
        i++;
    }

    float dt = 0.01;
    polyscope::state::userCallback = [&]() {
        if (ImGui::Button("reset")) {
            wc = init_wc(*m, *g, h);
            u = velocity(*m,*g,wc,h, S);
            pm->addVertexScalarQuantity("vorticity",wc.w);
            pm->addFaceTangentVectorQuantity("velocity",u,e1,e2);
        };
        ImGui::SliderFloat("delta time",&dt,0,1); ImGui::SameLine();
        if (!ImGui::Button("Advance")) {
            wc = RK4Step(*m,*g,h,wc, dt, S);
            u = velocity(*m,*g,wc,h, S);
            pm->addVertexScalarQuantity("vorticity",wc.w);
            pm->addFaceTangentVectorQuantity("velocity",u,e1,e2);
        }
        for (int i = 0; i< wc.c.size(); i++) {
            ImGui::Text("c%d: %f",i,wc.c[i]);
            if (i < wc.c.size() - 1) {ImGui::SameLine();}
        }

    };// specify the callback
    polyscope::show();

}

