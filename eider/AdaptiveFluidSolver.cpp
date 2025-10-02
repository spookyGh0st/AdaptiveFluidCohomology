#include "AdaptiveFluidSolver.h"

namespace geometrycentral::surface {

velocity_wrapper AdaptiveFluidSolver::velocity() {
    return ::geometrycentral::surface::velocity(tri.mesh(), tri.geom(), wc, h, S);
}

AdaptiveFluidSolver::AdaptiveFluidSolver(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, const AdaptiveFluidSolverData& data)
    : tri(mesh,geom,data.strategy), hom(tri.intrinsicTriangulation()), h(hom.harmonicBasis()), S(tri.mesh(),tri.geom()),
      conf(data.dopri5Conf), doerflerConf(data.doerflerConf), adapt_time(data.adaptive_time), adapte_space(data.adaptive_space), dt(data.dt),
      wc(VertexData<double>(tri.mesh()),std::vector<double>(hom.homologyB.size(),0)) // empty wc
{
    // On split, lerp w s.t. guess is closer to actual solution.
    tri.intrinsicTriangulation().edgeSplitCallbackList.push_back([&](Edge e, Halfedge he1, Halfedge he2) {
        this->wc.w[he1.vertex()] = 0.5 * this->wc.w[he1.tipVertex()] + 0.5 * this->wc.w[he2.tipVertex()];
    });
}

void AdaptiveFluidSolver::adapt() {
    tri.mesh().compress();

    AdaptiveTransfer transfer(tri.intrinsicTriangulation(), wc.w);

    VertexData<double> u(tri.mesh(), 0);
    // TODO: duplicate, but im unsure about indices
    S.compute(tri.mesh(), tri.geom());
    S.solve(tri.mesh(), tri.geom(), u, wc.w);
    FaceData<double> eta = poisson_residual_error_sqr(tri.mesh(), tri.geom(), u, wc.w);
    std::array<std::vector<Face>, 2> ref_coarse = select_doerfler(tri.mesh(), eta, doerflerConf);

    tri.refine(ref_coarse[0], &transfer);
    tri.coarse(ref_coarse[1], &transfer);

    // First, transfer the values to the new mesh, this requires noncompressed vertices!
    wc.w = transfer.transfer();

    // next, compress the mesh and update required quantities
    tri.mesh().compress();
    tri.intrinsicTriangulation().refreshQuantities();

    // on the new mesh, now recompute the rest of the mesh
    h = hom.harmonicBasis();
    S.compute(tri.mesh(), tri.geom());
}

DOPRI5_sample AdaptiveFluidSolver::step() {
    if(adapte_space) adapt();

    DOPRI5_sample dps;
    if (adapt_time){
        dps = adaptive_step(tri.mesh(), tri.geom(), h, wc, dt, S, conf);
        wc = dps.wc;
        dt = dps.t_future;
    }else {
        wc = RK4Step(tri.mesh(),tri.geom(),h,wc,dt,S);
        dps.wc = wc; dps.t_past = dt; dps.t_future = dt;
    }
    elapsed_time += dps.t_past;
    return dps;
}
AdaptiveFluidSolverData AdaptiveFluidSolver::data() const {
    return AdaptiveFluidSolverData(conf,doerflerConf,dt,adapt_time,adapte_space);
}
}
