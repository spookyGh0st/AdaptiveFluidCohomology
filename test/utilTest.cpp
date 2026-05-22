#include <gtest/gtest.h>

#include "eider/util.h"
#include <filesystem>
#include <geometrycentral/surface/meshio.h>

using namespace geometrycentral;
using namespace geometrycentral::surface;

TEST(utilTest, testGradient) {
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path() / "models" / "grid.stl";
    auto [m, g] = readManifoldSurfaceMesh(fds.string());
    g->requireFaceTangentBasis();
    g->requireHalfedgeVectorsInFace();
    VertexData<double> f(*m);
    for (Vertex v : m->vertices()) {
        f[v] = g->vertexPositions[v].x;
    }
    for (Face face : m->faces()) {
        Vector2 gr = grad(*g, face, f);
        Vector3 gr3d = g->faceTangentBasis[face][0] * gr[0] + g->faceTangentBasis[face][1] * gr[1];
        ASSERT_NEAR(gr3d.x, 1, 1e-12);
        ASSERT_NEAR(gr3d.y, 0, 1e-12);
        ASSERT_NEAR(gr3d.z, 0, 1e-12);
    }
}
