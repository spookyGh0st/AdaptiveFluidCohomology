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

TEST(afemTest, testSplitEdgePath)
{
  std::filesystem::path fds(__FILE__);
  fds = fds.parent_path()/ "models" /"quad.stl";
  auto [pm,pg] = readManifoldSurfaceMesh(fds.string());
  ASSERT_EQ(pm->nFaces(),2);
  ASSERT_EQ(pm->nEdges(),5);

  auto intrTri = IntegerCoordinatesIntrinsicTriangulation(*pm,*pg);
  ManifoldSurfaceMesh& mesh = *intrTri.intrinsicMesh;

  std::array<Edge,4> bEdges {}; int b_i = 0;
  Edge cEdge;
  for (Edge e: mesh.edges()) {
    if (e.isBoundary()) bEdges[b_i++] = e; else cEdge = e;
  }
  Vertex start_v = cEdge.halfedge().vertex();
  EdgeData<Halfedge> next(mesh,Halfedge());
  next[bEdges[0]]  = cEdge.halfedge();
  next[cEdge] = bEdges[1].halfedge();
  next[bEdges[1]] = bEdges[0].halfedge(); //points back

  intrTri.edgeSplitCallbackList.push_back([&](Edge e, Halfedge he1, Halfedge he2)
  {
    ASSERT_EQ(he1.vertex(), he2.vertex()) << "new he originate from the same vertex";
    // he1 and he2 lay along the old edge
    ASSERT_TRUE(he1.tipVertex() == start_v || he2.tipVertex() == start_v);
    ASSERT_EQ(next[e],  bEdges[1].halfedge()) << "original edge still points to the same he";

    std::array<Halfedge, 4> cyc_edges {
      he1.next(),
      he1.next().next().twin().next(),
      he2.next(),
      he2.next().next().twin().next(),
    };
    for (Halfedge be: cyc_edges) ASSERT_TRUE(be.edge().isBoundary());

    Halfedge in, out;
    for (Halfedge be: cyc_edges) {
      if (next[be.edge()] == Halfedge()) continue;
      if (next[be.edge()].edge() == e) in = be;
      else out = be;
    }
    ASSERT_NE(in,Halfedge());
    ASSERT_NE(out,Halfedge());
    ASSERT_NE(in,out);

    // 4 cases:
    if (in.tailVertex() == out.tipVertex())
      std::cout << "case 1" << std::endl;
    else if (in.tipVertex() == out.tailVertex())
      std::cout << "case 2" << std::endl;
    else
      std::cout << "case 3" << std::endl;
    });
  intrTri.splitEdge(cEdge.halfedge(),0.5);
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
