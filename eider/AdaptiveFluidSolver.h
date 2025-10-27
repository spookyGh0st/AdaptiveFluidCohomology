#pragma once

#include <utility>

#include "AdaptiveHomologyBasis.h"
#include "cfd.h"
#include "poisson.h"
#include "refine.h"

namespace geometrycentral::surface{

struct AdaptiveFluidSolverData {

    AdaptiveFluidSolverData( DOPRI5_conf d5c, DoeflerConf doc, double dt_, bool adaptive_time_, bool adaptive_space_, MARKING_STRATEGY strategy_, bool interpolate_harmonic_basis_)
        : dopri5Conf(d5c), doerflerConf(doc), dt(dt_), adaptive_time(adaptive_time_), adaptive_space(adaptive_space_), strategy(strategy_), interpolate_harmonic_basis(interpolate_harmonic_basis_) {}

    AdaptiveFluidSolverData(
        DOPRI5PresetConf d5c = DOPRI5PresetConf::HIGH,
        DoerflerPresetConf doc = DoerflerPresetConf::LOW,
        double dt_ = 0.0001,
        bool adaptive_time_ = true,
        bool adaptive_space_ = true,
        MARKING_STRATEGY strategy_ = MARKING_STRATEGY::PATTERN,
        bool interpolate_harmonic_basis_ = true
    )
        : dopri5Conf(DOPRI5Preset(d5c)),
          doerflerConf(DoerflerPreset(doc)),
          dt(dt_),
          adaptive_time(adaptive_time_),
          adaptive_space(adaptive_space_),
          strategy(strategy_),
          interpolate_harmonic_basis(interpolate_harmonic_basis_) {}

    DOPRI5_conf dopri5Conf = DOPRI5Preset(DOPRI5PresetConf::HIGH);
    DoeflerConf doerflerConf = DoerflerPreset(DoerflerPresetConf::LOW);
    double dt = 0.0001;
    bool adaptive_time = true;
    bool adaptive_space = true;
    MARKING_STRATEGY strategy = MARKING_STRATEGY::PATTERN;
    bool interpolate_harmonic_basis = true;
};

class AdaptiveFluidSolver {
  public:
    AdaptiveTriangulation tri;
    DOPRI5_conf conf;
    DoeflerConf doerflerConf;

    AdaptiveHomologyBasis hom;
    Harmonic_basis h; Harmonic_basis h_interpolated;

    wc_wrapper wc;
    double dt, elapsed_time = 0;
    StreamFunctionSolver S;

    bool adapt_time, adapte_space, interpolate_h;

    velocity_wrapper velocity();

    AdaptiveFluidSolver(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, const AdaptiveFluidSolverData& data);

    AdaptiveFluidSolverData data() const;

    void adapt();

    DOPRI5_sample step();
};

}

