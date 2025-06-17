#pragma once

#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include <geometrycentral/utilities/vector2.h>

#include "poisson.h"

namespace geometrycentral::surface
{
    struct wc_wrapper
    {
        VertexData<double> w;
        std::vector<double> c;
    };

    struct velocity_wrapper {
        FaceData<Vector2> u;
        FaceData<double> residual;
    };

    velocity_wrapper velocity(
        ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
        const wc_wrapper& wc, const std::vector<FaceData<Vector2>>& h, const StreamFunctionSolver& S
    );

    wc_wrapper RK4Step(
        ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
        const std::vector<FaceData<Vector2>>& h,
        const wc_wrapper& x, double dt, const StreamFunctionSolver& S
    );
}
