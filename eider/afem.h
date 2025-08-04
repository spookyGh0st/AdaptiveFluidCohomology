#pragma  once

#include "refine.h"
#include "poisson.h"
#include "geometrycentral/surface/transfer_functions.h"
#include "geometrycentral/surface/common_subdivision.h"

using namespace geometrycentral::surface;


void onSplit(EdgeData<Halfedge>& next, Edge e, Halfedge he1, Halfedge he2);


inline void adaptMesh(IntrinsicTriangulation& Tri,VertexData<double>& f, double theta) {
  auto A = AttributeTransfer(Tri);

  ManifoldSurfaceMesh &mesh = *Tri.intrinsicMesh;
  IntrinsicGeometryInterface &geom = Tri;

  Tri.edgeSplitCallbackList.emplace_back([&](Edge,Halfedge he1,Halfedge he2) {
    Vertex vb = he1.tailVertex(), ve = he2.tipVertex();
    double e1 = Tri.edgeLengths[he1.edge()], e2 = Tri.edgeLengths[he2.edge()];
    double f1 = e1 / (e1+e2), f2 = 1 - f1;
    f[he2.vertex()] = f[vb] * f1 + f[ve] * f2;
  });

  for (int i = 0; i < 10; ++i) {
    VertexData<double> u(mesh, 0);
    StreamFunctionSolver S;
    S.compute(mesh, Tri);
    S.solve(mesh, geom, u, f);
    FaceData<double> res = poisson_residual_error_sqr(mesh, geom, u, f);
    auto faces = select_doerfler(mesh, res, theta);
    refine(Tri, faces);
  }
}

