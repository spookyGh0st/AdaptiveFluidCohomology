#pragma once

#include <geometrycentral/surface/intrinsic_triangulation.h>
#include "homology.h"

using namespace geometrycentral::surface;

void onSplit(Edge e, Halfedge he1, Halfedge he2, HalfedgeData<std::optional<bool>> &nextLeft);

void onCollapse(Halfedge he, HalfedgeData<std::optional<bool>> &nextLeft);

class AdaptiveHomologyBasis{
    Homology_basis homologyB;
    ManifoldSurfaceMesh& mesh;
    IntrinsicGeometryInterface& geom;

  public:
    explicit AdaptiveHomologyBasis(IntrinsicTriangulation& icit);
    [[nodiscard]] Harmonic_basis harmonicBasis() const;
};
