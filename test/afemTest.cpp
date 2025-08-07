#include <eider/afem.h>

#include <filesystem>
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include <gtest/gtest.h>

#include "eider/homology.h"
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
  fds = fds.parent_path()/ "models" /"grid_hole.stl";
  auto [pm,pg] = readManifoldSurfaceMesh(fds.string());

  auto intrTri = IntegerCoordinatesIntrinsicTriangulation(*pm,*pg);
  intrTri.flipToDelaunay();
  intrTri.delaunayRefine(20);
  intrTri.refreshQuantities();
  ManifoldSurfaceMesh& mesh = *intrTri.intrinsicMesh;
  mesh.compress();
  auto homotopy_b= homotopy_basis(mesh,intrTri,mesh.face(0));
  auto homology_b = singular_homology_basis(mesh,homotopy_b);
  for (int h_idx = 0; h_idx< homology_b.size(); h_idx++) {
    intrTri.edgeSplitCallbackList.push_back([&,h_idx](Edge e, Halfedge he1, Halfedge he2) {
            onSplit(e, he1, he2, homology_b[h_idx].next, &homology_b[h_idx].start_e); });
  }

  auto refine_mesh = [&]()
  {
    VertexData<double> f(mesh, 1);
    VertexData<double> u(mesh, 0);
    StreamFunctionSolver S;
    S.compute(mesh, intrTri);
    S.solve(mesh, intrTri, u, f);
    FaceData<double> res = poisson_residual_error_sqr(mesh, intrTri, u, f);
    auto faces = select_doerfler(mesh, res, 0.5);
    refine(intrTri, faces);
  };


  // TODO edge split
  bool normalize_vectors = false;


  auto vis_mesh = [&]()
  {
    mesh.compress();
    VertexData<Vector3> int_positions(mesh) ;
    for (Vertex v : mesh.vertices()) { int_positions[v] = intrTri.vertexLocations[v].interpolate(pg->vertexPositions); }
    VertexPositionGeometry g(mesh, int_positions);

    polyscope::SurfaceMesh* polym = polyscope::registerSurfaceMesh("M", g.vertexPositions,mesh.getFaceVertexList(), polyscopePermutations(mesh));
    for (int h_idx = 0; h_idx< homology_b.size(); h_idx++)
    {
      auto& h = homology_b[h_idx];
      auto& n = h.next;
      HalfedgeData<int> nextInt(mesh,0);
      for (Edge e: mesh.edges())
        if (n[e] != Halfedge()) nextInt[n[e]] =1;
      polym->addHalfedgeScalarQuantity("h_b " + std::to_string(h_idx),nextInt);
    }

    auto orth_h_basis = orthonormal_hom_basis(mesh,intrTri,homology_b);
    g.requireFaceTangentBasis();
    FaceData<Vector3> e1(mesh),e2(mesh);
    for (Face f: mesh.faces()) { e1[f] = g.faceTangentBasis[f][0], e2[f] = g.faceTangentBasis[f][1];  }
    std::size_t i = 0;
    for (auto& b: orth_h_basis) {
      if (normalize_vectors) for (Face f: mesh.faces()) { b[f] *= 1./g.faceArea(f); }
      polym->addFaceTangentVectorQuantity("Hom basis" + std::to_string(i),b,e1,e2);
      i++;
    }
  };

  polyscope::init();
  vis_mesh();

  polyscope::state::userCallback = [&]()
  {
    if (ImGui::Button("Refine"))
    {
      refine_mesh();
      vis_mesh();
    }
    if (ImGui::Checkbox("Normalize",&normalize_vectors)) { vis_mesh(); }
  };

  polyscope::show();
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
