#include <eider/afem.h>

#include <filesystem>
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include <gtest/gtest.h>

#include "eider/util.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

TEST(afemTest, testDefaultBehauvior)
{
  std::filesystem::path fds(__FILE__);
  fds = fds.parent_path()/ "models" /"L.stl";
  auto [pm,pg] = readManifoldSurfaceMesh(fds.string());
  auto intrTri = IntegerCoordinatesIntrinsicTriangulation(*pm,*pg);
  ManifoldSurfaceMesh& m = *intrTri.intrinsicMesh;
  VertexData<double> f(m,0);
  Face face = m.face(0);
  f[*face.adjacentVertices().begin()] = 1;
  intrTri.faceInsertionCallbackList.push_back([&](Face face, Vertex v) {f[v] = 2; });
  Vertex v = intrTri.splitFace(face,Vector3(0.33,0.333,0.33));

  ASSERT_EQ(f[v], 2);
}

TEST(afemTest, testL)
{
  std::filesystem::path fds(__FILE__);
  fds = fds.parent_path()/ "models" /"L.stl";
  auto [pm,pg] = readManifoldSurfaceMesh(fds.string());
  auto intrTri = IntegerCoordinatesIntrinsicTriangulation(*pm,*pg);

  ManifoldSurfaceMesh& m = *intrTri.intrinsicMesh;
  m.compress();

  VertexData<double> f(m, 1);
  adaptMesh(intrTri,f,0.1);

  VertexData<Vector3> int_positions(m) ;
  for (Vertex v : m.vertices()) {
    int_positions[v] = intrTri.vertexLocations[v].interpolate(pg->vertexPositions);
  }
  VertexPositionGeometry g(m, int_positions);

  polyscope::init();
  polyscope::SurfaceMesh* polym = polyscope::registerSurfaceMesh("M", g.vertexPositions,m.getFaceVertexList(), polyscopePermutations(m));
  polym->addVertexScalarQuantity("f",f);
  polyscope::show();
}
