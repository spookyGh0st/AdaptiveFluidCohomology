#include "TaylorInitializer.h"
#include "imgui.h"

namespace geometrycentral::surface {

// Intrinsic geometry -> 2D UV coordinates scaled into [0,1)x[0,1)
VertexData<Vector2> xy_plane(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& inputG) {
    VertexData<Vector2> uv(mesh);
    Vector2 min = Vector2::infinity();
    Vector2 max = -Vector2::infinity();
    for (Vertex v : mesh.vertices()) {
        Vector3 p = inputG.vertexPositions[v];
        uv[v] = Vector2(p.x, p.y);
        min = componentwiseMin(uv[v],min);
        max = componentwiseMax(uv[v],max);
    }

    Vector2 scale = max - min;
    if (scale.x == 0)
        scale.x = 1;
    if (scale.y == 0)
        scale.y = 1;

    // Normalize into [-1,1)x[-1,1)
    double mscale = std::max(scale.x,scale.y)/2;
    for (Vertex v : mesh.vertices()) {
        uv[v].x = (uv[v].x) / mscale;
        uv[v].y = (uv[v].y) / mscale;
    }

    return uv;
}

// Intrinsic geometry -> 2D UV coordinates scaled into [0,1)x[0,1)
VertexData<Vector2> xy_plane(IntrinsicTriangulation& Tri, VertexPositionGeometry& inputG) {
    auto &mesh = *Tri.intrinsicMesh;
    VertexData<Vector2> uv(mesh);


    Vector2 min = Vector2::infinity();
    Vector2 max = -Vector2::infinity();
    for (Vertex v : mesh.vertices()) {
        Vector3 p = Tri.vertexLocations[v].interpolate(inputG.vertexPositions);
        uv[v] = Vector2(p.x, p.y);
        min = componentwiseMin(uv[v],min);
        max = componentwiseMax(uv[v],max);
    }

    Vector2 scale = max - min;
    if (scale.x == 0)
        scale.x = 1;
    if (scale.y == 0)
        scale.y = 1;

    // Normalize into [-1,1)x[-1,1)
    double mscale = std::max(scale.x,scale.y)/2;
    for (Vertex v : mesh.vertices()) {
        uv[v].x = (uv[v].x) / mscale;
        uv[v].y = (uv[v].y) / mscale;
    }

    return uv;
}

wc_wrapper init_taylor(SurfaceMesh &mesh, VertexData<Vector2> geo, double vorticity_distance, const EigenVec2 &offset, const EigenVec2 &box_min, const EigenVec2 &box_max) {
    VertexData<double> w(mesh, 0.0);
    double k = 2 * PI / vorticity_distance;
    double A = 1.0;

    for (Vertex v : mesh.vertices()) {
        EigenVec2 pos(geo[v].x,geo[v].y);
        // Bounding box check
        if (pos.x() < box_min.x() || pos.x() > box_max.x() ||
            pos.y() < box_min.y() || pos.y() > box_max.y()) {
            continue;
        }

        // Apply offset
        EigenVec2 shifted = pos + offset;

        // Taylor vortex field
        w[v] = 2 * A * k * std::cos(k * shifted.x()) * std::cos(k * shifted.y());
    }

    wc_wrapper wc;
    wc.w = w;
    return wc;
}

void TaylorInitializer::callback() {
    if (ImGui::TreeNode("Initial Taylor Vortices")) {
        ImGui::InputDouble("vorticity distance", &v_dist, 0.125, 0.5);
        ImGui::InputFloat2("center", center.data());
        set_vortexPair(v_dist, toGC(center));
        ImGui::TreePop();
    }
}
wc_wrapper TaylorInitializer::wc(IntrinsicTriangulation &intTri, VertexPositionGeometry &pg) {
    return init_taylor(*intTri.intrinsicMesh, xy_plane(intTri,pg), v_dist, offset, box_min, box_max);
}
wc_wrapper TaylorInitializer::wc(ManifoldSurfaceMesh &mesh, VertexData<Vector2> &uv) {
    return init_taylor(mesh, uv, v_dist, offset, box_min, box_max);
}

void TaylorInitializer::set_vortexPair(double new_v_dist, const Vector2 &center) {
    v_dist = new_v_dist;
    // Region: 1 wavelength wide in x, 2 wavelengths tall in y
    box_min = {center.x - v_dist * 0.5, center.y - v_dist*0.25};
    box_max = {center.x + v_dist * 0.5, center.y + v_dist*0.25};

    offset = EigenVec2(v_dist/4.f,0) ;
}
wc_wrapper TaylorInitializer::wc(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &pg) {
    return init_taylor(mesh, xy_plane(mesh,pg), v_dist, offset, box_min, box_max);
}

}
