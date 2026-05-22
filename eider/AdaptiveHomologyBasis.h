#pragma once

#include "homology.h"
#include <geometrycentral/surface/intrinsic_triangulation.h>

namespace geometrycentral::surface {

void onSplit(Edge e, Halfedge he1, Halfedge he2, HalfedgeData<std::optional<bool>> &nextLeft);

void onCollapse(Halfedge he, HalfedgeData<std::optional<bool>> &nextLeft);

struct Harmonic_Data {
    Homology_basis hom;
    std::vector<EdgeData<double>> df;
    std::vector<EdgeData<double>> proj_df;
    Harmonic_basis h_unorth;
    Harmonic_basis h_orth;
    Harmonic_Data(Homology_basis hom) : hom(hom), df(hom.size()), proj_df(hom.size()), h_unorth(hom.size()), h_orth(hom.size()) {}
};

class AdaptiveHomologyBasis {
    ManifoldSurfaceMesh &mesh;
    IntrinsicGeometryInterface &geom;

  public:
    PressureProjectionSolver pp_solver;
    std::vector<EdgeData<double>> pf_guess;
    std::vector<EdgeData<double>> pf_guess_L2;
    Homology_basis homologyB;
    explicit AdaptiveHomologyBasis(IntrinsicTriangulation &icit);
    [[nodiscard]] Harmonic_basis harmonicBasis();
    [[nodiscard]] Harmonic_Data fullHarmonicBasis();
};

} // namespace geometrycentral::surface
