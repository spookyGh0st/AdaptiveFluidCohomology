#include "Stopwatch.h"

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

wc_wrapper init_taylor(SurfaceMesh& mesh, const VertexPositionGeometry& geo, const std::vector<FaceData<Vector2>>& h, double vorticity_distance, Vector2 offset, Vector2 cuttof, Vector2 cuttof_offset) {
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
VertexPositionGeometry intrinsic_geom(IntrinsicTriangulation& Tri, VertexPositionGeometry& inputG){
    auto& mesh = *Tri.intrinsicMesh;
    VertexData<Vector3> int_positions(mesh) ;
    for (Vertex v : mesh.vertices()) { int_positions[v] = Tri.vertexLocations[v].interpolate(inputG.vertexPositions); }
    return VertexPositionGeometry(mesh, int_positions);
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
    AdaptiveTriangulation atri(intrTri);
    ManifoldSurfaceMesh& mesh = atri.mesh();
    mesh.compress();
    auto homotopy_b= greedy_homotopy_basis(mesh,intrTri,arbitrary_base_face(mesh));
    auto homology_b = singular_homology_basis(mesh,homotopy_b);
    auto harmonic_b = orthonormal_hom_basis(mesh,intrTri,homology_b);

    float dt = 0.001, v_dist  = 0.5, x_add = v_dist/4, y_add = 0, x_cuttof = v_dist/2, y_cuttof = v_dist/4, x_cutt_offset = 0, y_cutt_offset = 0.5;
    wc_wrapper wc = init_taylor(mesh, intrinsic_geom(intrTri,*pg), harmonic_b, v_dist, Vector2(x_add, y_add), Vector2(x_cuttof, y_cuttof), Vector2(x_cutt_offset, y_cutt_offset));

    double error_threshold = 0.0001;
    double theta = 0.1;
    auto refine_mesh = [&]()
    {
      // TODO: adaptMesh(intrTri,wc,homology_b,theta, error_threshold);
      harmonic_b = orthonormal_hom_basis(mesh,intrTri,homology_b);
    };
    auto coarse_mesh = [&]()
    {
      atri.coarse([](Vertex v) {return true;});
      mesh.compress();
      harmonic_b = orthonormal_hom_basis(mesh,intrTri,homology_b);
    };


    // TODO edge split
    bool normalize_vectors = false;

    DOPRI5_conf conf=  DOPRI5_conf();
    bool running = false, fix_c = false, force_c = false, adaptive_run = false;

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
          auto& n = h.nextLeft;
          HalfedgeData<int> nextInt(mesh,0);
          for (Halfedge e: mesh.halfedges()){
              if (n[e].has_value() && *n[e]) nextInt[e] =1;
              if (n[e].has_value() && !(*n[e])) nextInt[e] =2;
          }

          polym->addHalfedgeScalarQuantity("h_b " + std::to_string(h_idx),nextInt);
      }


      g.requireFaceTangentBasis();
      FaceData<Vector3> e1(mesh),e2(mesh);
      for (Face f: mesh.faces()) { e1[f] = g.faceTangentBasis[f][0], e2[f] = g.faceTangentBasis[f][1];  }
      std::size_t i = 0;
      for (auto& b: harmonic_b) {
          if (normalize_vectors) for (Face f: mesh.faces()) { b[f] *= 1./g.faceArea(f); }
          polym->addFaceTangentVectorQuantity("Hom basis" + std::to_string(i),b,e1,e2);
          i++;
      }

      auto S =StreamFunctionSolver();
      S.compute(mesh,intrTri);
      auto v = velocity(mesh,intrTri,wc,harmonic_b,S);
      polym->addFaceTangentVectorQuantity("velocity",v.u,e1,e2);
      polym->addVertexScalarQuantity("vorticity",wc.w);
      polym->addFaceScalarQuantity("residual",v.residual);
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
      if (ImGui::Button("Coarse")) {
          coarse_mesh(); vis_mesh();
      }
      if (ImGui::Checkbox("Normalize",&normalize_vectors)) { vis_mesh(); }
      ImGui::Checkbox("Run", &running);
      ImGui::Checkbox("Refine while running", &adaptive_run);
      if (running || ImGui::Button("Advance")) {
          auto tmpc = wc.c; double tmpdt;
          StreamFunctionSolver S;
          S.compute(mesh,intrTri);
          auto dp5_sample = adaptive_step(mesh,intrTri,harmonic_b,wc, dt, S,conf);
          wc = dp5_sample.wc;
          dt = float(dp5_sample.t_future);

          if (adaptive_run) refine_mesh();
          vis_mesh();
      }
      ImGui::Text("dt: %.6f", dt);
      for (auto& c: wc.c) {
          ImGui::Text("c: %.6f", c);
      }
    };

    polyscope::show();
}

TEST(afemTest, testPathConsistency){
    std::filesystem::path fds(__FILE__);
    fds = fds.parent_path()/ "models" / "quad.stl";
    auto [parent_m,parent_g] = readManifoldSurfaceMesh(fds.string());

    IntegerCoordinatesIntrinsicTriangulation icit(*parent_m,*parent_g);
    AdaptiveTriangulation atri(icit);
    ManifoldSurfaceMesh& m = atri.mesh();
    IntrinsicGeometryInterface& g = atri.geom();

    auto homotopy_b= greedy_homotopy_basis(m,g, arbitrary_base_face(m));
    auto homologyBasis = singular_homology_basis(m,homotopy_b);

    auto& h = homologyBasis;
    for (int h_idx = 0; h_idx < h.size(); ++h_idx) {
        icit.edgeSplitCallbackList.emplace_back([&,h_idx](Edge e, Halfedge he1, Halfedge he2) {
          onSplit(e, he1, he2, h[h_idx].nextLeft); });
    }
    for (int h_idx = 0; h_idx < h.size(); ++h_idx) {
        icit.edgeCollapseCallbackList.push_back([&,h_idx](Halfedge he) {
          onCollapse(he, h[h_idx].nextLeft); });
    }

    bool first = true;
    auto vis = [& ](){
      m.compress();
      VertexData<Vector3> int_positions(m) ;
      for (Vertex v : m.vertices()) {
          int_positions[v] = icit.vertexLocations[v].interpolate(parent_g->vertexPositions);
      }



      polyscope::SurfaceMesh* pm_A = polyscope::registerSurfaceMesh("original", parent_g->vertexPositions,parent_m->getFaceVertexList());
      polyscope::SurfaceMesh* pm_B = polyscope::registerSurfaceMesh("Intrinsic", int_positions,m.getFaceVertexList());
      pm_B->setAllPermutations(polyscopePermutations(m));
      pm_A->addFaceScalarQuantity("index",parent_m->getFaceIndices());
      pm_B->addFaceScalarQuantity("index",m.getFaceIndices());
      pm_B->addEdgeScalarQuantity("edgeCoords",icit.normalCoordinates.edgeCoords);
      pm_B->addHalfedgeScalarQuantity("roundabouts",icit.normalCoordinates.roundabouts);
      pm_B->addVertexScalarQuantity("roundabouts degree",icit.normalCoordinates.roundaboutDegrees);
      pm_B->addEdgeScalarQuantity("lenghts",icit.edgeLengths);
      pm_A->setEdgeWidth(1);
      pm_B->setEdgeWidth(1);
      if(first){
          pm_A->translate(glm::vec3(-1,0,0));
          pm_B->translate(glm::vec3( 1,0,0)); first = false;
      }
      HalfedgeData<int> d(m,0);
      for (int i = 0; i < homologyBasis.size(); ++i) {
          for (Halfedge e: m.halfedges()){
              if(homologyBasis[i].nextLeft[e].has_value()) {
                  if (*homologyBasis[i].nextLeft[e]) d[e] = i*2+1; else d[e] = i*2+2;
                  if (!e.isInterior()) d[e.twin()] = i*2+3;
              }
          }
      }
      pm_B->addHalfedgeScalarQuantity("Hom", d);
    };


    polyscope::init();
    std::vector<Face> faces {  };
    for (Face f: m.faces()) faces.push_back(f);
    vis();

    polyscope::state::userCallback = [&]()
    {
      if (ImGui::Button("Refine"))
      {
          std::vector<Face> faces {  };
          for (Face f: m.faces()) faces.push_back(f);
          atri.refine(faces);
          vis();
      }
      if (ImGui::Button("Coarse")) {
          atri.coarse([](Vertex v)->bool { return true;});
          vis();
      }
    };

    polyscope::show();
}

