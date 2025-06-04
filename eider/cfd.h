#pragma once

#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include <geometrycentral/utilities/vector2.h>

namespace geometrycentral::surface
{
    struct wc_wrapper
    {
        VertexData<double> w;
        std::vector<double> c;
    };

    FaceData<Vector2> velocity(
        SurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
        const wc_wrapper& wc, const std::vector<FaceData<Vector2>>& h
    );

    wc_wrapper RK4Step(
        SurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
        const std::vector<FaceData<Vector2>>& h,
        const wc_wrapper& x, double dt
    );
}
