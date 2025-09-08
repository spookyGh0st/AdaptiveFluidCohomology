#pragma once

#include <geometrycentral/surface/intrinsic_triangulation.h>
#include "homology.h"

namespace geometrycentral::surface {

void onSplit(Edge e, Halfedge he1, Halfedge he2, HalfedgeData<std::optional<bool>> &nextLeft);

void onCollapse(Halfedge he, HalfedgeData<std::optional<bool>> &nextLeft);

class AdaptiveHomologyBasis{
    ManifoldSurfaceMesh& mesh;
    IntrinsicGeometryInterface& geom;

  public:
    PressureProjectionSolver pp_solver;
    std::vector<EdgeData<double>> pf_guess;
    std::vector<EdgeData<double>> pf_guess_L2;
    Homology_basis homologyB;
    explicit AdaptiveHomologyBasis(IntrinsicTriangulation& icit);
    [[nodiscard]] Harmonic_basis harmonicBasis();
};

}
