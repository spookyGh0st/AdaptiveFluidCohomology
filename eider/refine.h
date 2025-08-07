#pragma once
#include <geometrycentral/surface/intrinsic_triangulation.h>
#include <vector>

namespace gcs = geometrycentral::surface;

void refine(gcs::IntrinsicTriangulation &tri, std::vector<gcs::Face> faces);
/// Gives the residual error for Δu = f in Ω, u = 0 in ∂Ω
gcs::FaceData<double> poisson_residual_error_sqr(
    gcs::ManifoldSurfaceMesh &mesh,
    gcs::IntrinsicGeometryInterface &geom,
    const gcs::VertexData<double> &u,
    const gcs::VertexData<double> &f);
std::vector<gcs::Face> select_doerfler(gcs::ManifoldSurfaceMesh &mesh, gcs::FaceData<double> residual, double theta);
