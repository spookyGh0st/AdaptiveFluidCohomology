#pragma once

#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/surface_mesh.h"
#include <geometrycentral/utilities/vector2.h>

#include "poisson.h"

namespace geometrycentral::surface {
struct wc_wrapper {
    VertexData<double> w;
    std::vector<double> c;
};

struct velocity_wrapper {
    VertexData<double> stream_function;
    FaceData<Vector2> u;
};

velocity_wrapper velocity(
    ManifoldSurfaceMesh &mesh,
    IntrinsicGeometryInterface &geom,
    const wc_wrapper &wc,
    const std::vector<FaceData<Vector2>> &h,
    const StreamFunctionSolver &S);

wc_wrapper evalRHS(
    ManifoldSurfaceMesh &mesh,
    IntrinsicGeometryInterface &geom,
    const wc_wrapper &wc,
    const std::vector<FaceData<Vector2>> &h,
    const StreamFunctionSolver &S,
    std::vector<FaceData<double>>* face_dc = nullptr
);

wc_wrapper RK4Step(
    ManifoldSurfaceMesh &mesh,
    IntrinsicGeometryInterface &geom,
    const std::vector<FaceData<Vector2>> &h,
    const wc_wrapper &x,
    double dt,
    const StreamFunctionSolver &S);

/**
 * @brief Configuration parameters for DOPRI5 adaptive solver.
 */
struct DOPRI5_conf {
    /**
     * @brief Absolute tolerance.
     */
    double Atol_i = 1e-12;

    /** * @brief Relative tolerance. */
    double Rtol_i = 1e-6;

    /** * @brief Max factor to increase step size (1.5–5.0 typical). */
    double faxmax = 3.0;
    /** @brief Min factor to reduce step size (e.g. 0.1–0.2). */
    double facmin = 0.1;
};

struct DOPRI5_sample {
    wc_wrapper wc{};
    double t_past{};   // past step
    double t_future{}; // new stepsize
    int attempts = 0;
};

DOPRI5_sample adaptive_step(
    ManifoldSurfaceMesh &mesh,
    IntrinsicGeometryInterface &geom,
    const std::vector<FaceData<Vector2>> &h,
    const wc_wrapper &x,
    double dt,
    const StreamFunctionSolver &S,
    const DOPRI5_conf &conf);
} // namespace geometrycentral::surface
