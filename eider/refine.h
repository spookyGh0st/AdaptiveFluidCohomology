#pragma once
#include <geometrycentral/surface/intrinsic_triangulation.h>
#include <vector>

namespace gcs=geometrycentral::surface;

void refine(gcs::IntrinsicTriangulation& tri, std::vector<gcs::Face> faces);
gcs::FaceData<double> poisson_residual_error(gcs::ManifoldSurfaceMesh& mesh, gcs::IntrinsicGeometryInterface& geom, const gcs::VertexData<double>& f, const gcs::VertexData<double>& u);
