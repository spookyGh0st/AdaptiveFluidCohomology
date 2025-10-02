#pragma once

#include <utility>

#include "AdaptiveHomologyBasis.h"
#include "cfd.h"
#include "poisson.h"
#include "refine.h"

namespace geometrycentral::surface{



struct AdaptiveFluidSolverData{
    DOPRI5_conf dopri5Conf = DOPRI5Preset(DOPRI5PresetConf::HIGH);
    DoeflerConf doerflerConf = DoerflerPreset(DoerflerPresetConf::LOW);
    double dt = 0.0001;
    bool adaptive_time = true;
    bool adaptive_space = true;
    MARKING_STRATEGY strategy = MARKING_STRATEGY::PATTERN;
};

class AdaptiveFluidSolver {
  public:
    AdaptiveTriangulation tri;
    DOPRI5_conf conf;
    DoeflerConf doerflerConf;

    AdaptiveHomologyBasis hom;
    Harmonic_basis h;

    wc_wrapper wc;
    double dt, elapsed_time = 0;
    StreamFunctionSolver S;

    bool adapt_time, adapte_space;

    velocity_wrapper velocity();

    AdaptiveFluidSolver(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, const AdaptiveFluidSolverData& data);

    AdaptiveFluidSolverData data() const;

    void adapt();

    DOPRI5_sample step();
};

}

