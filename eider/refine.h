#pragma once
#include <geometrycentral/surface/intrinsic_triangulation.h>
#include <vector>

namespace geometrycentral::surface {

void refine(IntrinsicTriangulation &tri, std::vector<Face> faces);

// get halfedge that was used in the edgesplit resulting in vertex v
Halfedge coarse_halfedge(Vertex v);

void coarse(IntrinsicTriangulation& m, const std::function<bool(Vertex)>& f);

/// Gives the residual error for Δu = f in Ω, u = 0 in ∂Ω
FaceData<double> poisson_residual_error_sqr(
    ManifoldSurfaceMesh &mesh,
    IntrinsicGeometryInterface &geom,
    const VertexData<double> &u,
    const VertexData<double> &f);

std::vector<Face> select_doerfler(ManifoldSurfaceMesh &mesh, FaceData<double> residual, double theta, double threshold);
}
