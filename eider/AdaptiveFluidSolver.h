#pragma once

#include <utility>

#include "AdaptiveHomologyBasis.h"
#include "cfd.h"
#include "poisson.h"
#include "refine.h"

namespace geometrycentral::surface{

class AdaptiveFluidSolver {
  public:
    AdaptiveTriangulation& tri;
    wc_wrapper wc;
    DOPRI5_conf conf;
    DoeflerConf doerflerConf;

    AdaptiveHomologyBasis hom;
    Harmonic_basis h;

    double dt, elapsed_time = 0;
    StreamFunctionSolver S;

    velocity_wrapper velocity() const;

    AdaptiveFluidSolver(AdaptiveTriangulation &tri, wc_wrapper wc, const DOPRI5_conf &conf, const DoeflerConf &doerflerConf);

    void adapt();

    void step();
};

}

