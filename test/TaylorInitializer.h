#pragma  once
#include "geometrycentral/surface/intrinsic_triangulation.h"
#include <eider/cfd.h>

namespace geometrycentral::surface {

using EigenVec2 = Eigen::Vector2f;
// Convert GC Vector2 -> Eigen::Vector2d
inline EigenVec2 toEigen(const Vector2& v) {
    return EigenVec2{v.x, v.y};
}

// Convert Eigen::Vector2d -> GC Vector2
inline Vector2 toGC(const EigenVec2 & v) {
    return Vector2{v.x(), v.y()};
}

struct TaylorInitializer {
    EigenVec2 center = { 0.0,-0.5 };
    double v_dist = 0.5;
    EigenVec2 offset = {v_dist / 4.0, 0.0};
    EigenVec2 box_min = {0,0};
    EigenVec2 box_max = {1,1};
    TaylorInitializer() { set_vortexPair(v_dist, toGC(center)); }

    void set_vortexPair(double new_v_dist, const Vector2& center);

    wc_wrapper wc(IntrinsicTriangulation &intTri, VertexPositionGeometry &pg);
    wc_wrapper wc(ManifoldSurfaceMesh &mesh, VertexData<Vector2> &uv);

    void callback();
};
}

