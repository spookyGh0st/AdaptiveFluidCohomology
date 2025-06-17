#include <eider/homology.h>
#include <eider/poisson.h>

#include <filesystem>
#include "geometrycentral/surface/meshio.h"
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include <gtest/gtest.h>

#include "eider/util.h"

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

TEST(poissonTest, testLaplaceInverse)
{
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" /"grid_fine.stl";
    auto [m,g] = readManifoldSurfaceMesh(fds.string());

    VertexData<double> f(*m, 0);
    // for (Vertex v: m->vertices()) f[v] = g->vertexPositions[v].x;
    VertexData<double> g2(*m, 1);
    for (Vertex v: m->vertices()) if (v.isBoundary()) g2[v] = 0;
    auto S = StreamFunctionSolver();
    S.compute(*m, *g);
    S.solve(*m,*g,f,g2);

    g->requireHalfedgeCotanWeights();
    g->requireEdgeCotanWeights();
    g->requireVertexDualAreas();
    g->requireCotanLaplacian();

    VertexData<double> gMatrix(*m, g->cotanLaplacian * f.toVector()), gVertex(*m);
    for (Vertex v: m->vertices()) gVertex[v] = laplacian(*g,v,f);

    for (Vertex v: m->vertices()) {
        ASSERT_NEAR(gMatrix[v], gVertex[v], 0.01);
    }
}
