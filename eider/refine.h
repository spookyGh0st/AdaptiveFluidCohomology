#pragma once
#include <geometrycentral/surface/intrinsic_triangulation.h>
#include <vector>

namespace geometrycentral::surface {

struct IncrementingIndex {
    static constexpr std::size_t invalidIdx = std::numeric_limits<std::size_t>::max();
    std::size_t& operator[](Face f) {
        if (idx[f] == invalidIdx) idx[f] = current++;
        return idx[f];
    }
    IncrementingIndex(ManifoldSurfaceMesh& mesh): idx(mesh,invalidIdx) {}
private:
    FaceData<std::size_t> idx;
    std::size_t current = 0;;
};

class AdaptiveTriangulation {
public:
    AdaptiveTriangulation(IntrinsicTriangulation& tri);
    Halfedge vertex_bisection(Halfedge he);
    void refine(std::vector<Face> faces);


    // get halfedge that was used in the edgesplit resulting in vertex v
    Halfedge coarse_halfedge(Vertex v);
    Halfedge vertex_biunion(Halfedge he);
    void coarse(const std::vector<Face> &f);
    ManifoldSurfaceMesh& mesh() const { return *tri.intrinsicMesh; }
    IntrinsicGeometryInterface& geom() { return tri; }
    IncrementingIndex idx;
    CornerData<bool> marked_corner;
    Halfedge getRefinementEdge(Face f);
private:
    void setRefinementEdge(Halfedge he);
    IntrinsicTriangulation& tri;
};

/// Gives the residual error for Δu = f in Ω, u = 0 in ∂Ω
FaceData<double> poisson_residual_error_sqr(
    ManifoldSurfaceMesh &mesh,
    IntrinsicGeometryInterface &geom,
    const VertexData<double> &u,
    const VertexData<double> &f);

std::vector<Face> select_doerfler(ManifoldSurfaceMesh &mesh, FaceData<double> residual, double theta, double threshold);
}
