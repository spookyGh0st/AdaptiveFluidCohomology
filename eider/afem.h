#pragma once

#include "geometrycentral/surface/common_subdivision.h"
#include "geometrycentral/surface/transfer_functions.h"
#include "poisson.h"
#include "refine.h"
#include "cfd.h"
#include "homology.h"

using namespace geometrycentral::surface;

void onSplit(Edge e, Halfedge he1, Halfedge he2, HalfedgeData<std::optional<bool>> &nextLeft);

void onCollapse(Halfedge he, HalfedgeData<std::optional<bool>> &nextLeft);

inline void adaptMesh(IntrinsicTriangulation &Tri, wc_wrapper& wc, Homology_basis& h, double theta, double threshold) {
    ManifoldSurfaceMesh &mesh = *Tri.intrinsicMesh;
    IntrinsicGeometryInterface &geom = Tri;

    // TODO: DEBUG only
    Tri.edgeSplitCallbackList.resize(1);

    Tri.edgeSplitCallbackList.emplace_back([&](Edge, Halfedge he1, Halfedge he2) {
      Vertex vb = he1.tipVertex(), ve = he2.tipVertex();
      double e1 = Tri.edgeLengths[he1.edge()], e2 = Tri.edgeLengths[he2.edge()];
      double f1 = e1 / (e1 + e2), f2 = 1 - f1;
      wc.w[he2.vertex()] = wc.w[vb] * f1 + wc.w[ve] * f2;
    });

    for (int h_idx = 0; h_idx < h.size(); ++h_idx) {
        Tri.edgeSplitCallbackList.push_back([&,h_idx](Edge e, Halfedge he1, Halfedge he2) {
          onSplit(e, he1, he2, h[h_idx].nextLeft); });
    }
    for (int h_idx = 0; h_idx < h.size(); ++h_idx) {
        Tri.edgeCollapseCallbackList.push_back([&,h_idx](Halfedge he) {
          onCollapse(he, h[h_idx].nextLeft); });
    }

    VertexData<double> u(mesh, 0);
    StreamFunctionSolver S;
    S.compute(mesh, Tri);
    S.solve(mesh, geom, u, wc.w);
    FaceData<double> res = poisson_residual_error_sqr(mesh, geom, u, wc.w);
    auto faces = select_doerfler(mesh, res, theta, threshold);
    refine(Tri, faces);
}
