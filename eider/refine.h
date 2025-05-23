#pragma once
#include <geometrycentral/surface/intrinsic_triangulation.h>
#include <vector>

namespace gcs=geometrycentral::surface;

void refine(gcs::IntrinsicTriangulation& tri, const std::vector<gcs::Face>& faces);
