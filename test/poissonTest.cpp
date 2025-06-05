#include <eider/homology.h>
#include <eider/poisson.h>

#include <filesystem>
#include "geometrycentral/surface/meshio.h"
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include <gtest/gtest.h>

using namespace geometrycentral;
using namespace geometrycentral::surface;


TEST(poissonTest, testLaplace)
{
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" /"grid.stl";
    auto [m,g] = readManifoldSurfaceMesh(fds.string());

    VertexData<double> f(*m, 0);
    for (Vertex v: m->vertices()) f[v] = g->vertexPositions[v].x;
    VertexData<double> g2(*m, 1);
    for (Vertex v: m->vertices()) g2[v] = std::sin(g->vertexPositions[v].x);
    auto S = StreamFunctionSolver();
    S.compute(*m, *g);
    S.solve(*m,*g,f,g2);

    polyscope::init();
    polyscope::SurfaceMesh* pm = polyscope::registerSurfaceMesh("M", g->vertexPositions,m->getFaceVertexList());
    pm->addVertexScalarQuantity("f",f)->setEnabled(true);
    pm->addVertexScalarQuantity("g",g2);
    std::size_t i = 0;
    polyscope::show();

}
