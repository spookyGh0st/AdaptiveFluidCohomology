#pragma once
#include "transfer.h"

#include <geometrycentral/surface/intrinsic_triangulation.h>
#include <vector>

namespace geometrycentral::surface {

struct IncrementingIndex {
    static constexpr std::size_t invalidIdx = std::numeric_limits<std::size_t>::max();
    std::size_t &operator[](Face f);
    explicit IncrementingIndex(ManifoldSurfaceMesh &mesh);

    // Compress the Index List, such that it goes from 0 to nFaces, but still retains same order
    void compress();

  private:
    FaceData<std::size_t> idx;
    std::size_t current = 0;
    ;
};

class AdaptiveTriangulation {
  public:
    AdaptiveTriangulation(IntrinsicTriangulation &tri);
    Halfedge vertex_bisection(Halfedge he, AdaptiveTransfer *transfer = nullptr);
    void refine(std::vector<Face> faces, AdaptiveTransfer *transfer = nullptr);

    // get halfedge that was used in the edgesplit resulting in vertex v
    Halfedge coarse_halfedge(Vertex v);
    Halfedge vertex_biunion(Halfedge he, AdaptiveTransfer *transfer = nullptr);
    void coarse(const std::vector<Face> &f, AdaptiveTransfer *transfer = nullptr);

    ManifoldSurfaceMesh &mesh() const { return *tri.intrinsicMesh; }
    IntrinsicGeometryInterface &geom() const { return tri; }
    IntrinsicTriangulation &intrinsicTriangulation() const { return tri; }

    IncrementingIndex idx;
    CornerData<bool> marked_corner;
    Halfedge getRefinementEdge(Face f);

    Eigen::SparseMatrix<double> M_restriction;

  private:
    void setRefinementEdge(Halfedge he);
    IntrinsicTriangulation &tri;
};

/// Gives the residual error for Δu = f in Ω, u = 0 in ∂Ω
FaceData<double> poisson_residual_error_sqr(
    ManifoldSurfaceMesh &mesh,
    IntrinsicGeometryInterface &geom,
    const VertexData<double> &u,
    const VertexData<double> &f);

/**
 * @brief Configuration parameters for the Dörfler marking criterion.
 *
 * The Dörfler criterion selects a subset of elements (faces) such that the
 * cumulative residual contribution exceeds a prescribed fraction θ of the
 * total residual. This is commonly used in adaptive mesh refinement/coarsening.
 *
 * - θ_refine ∈ (0,1]  : fraction of residual to mark for refinement
 * - θ_coarse ∈ (0,1]  : fraction of residual to mark for coarsening
 * - threshold_refine  : absolute lower bound for refinement marking
 * - threshold_coarse  : absolute upper bound for coarsening marking
 */
struct DoeflerConf {
    double theta_coarse = 0.8;                                              ///< θ_coarse: fraction of residual for coarsening
    double theta_refine = 0.2;                                              ///< θ_refine: fraction of residual for refinement
    double threshold_refine = std::numeric_limits<double>::epsilon(); ///< refinement cutoff
    double threshold_coarse = std::numeric_limits<double>::max();     ///< coarsening cutoff
};

/**
 * @brief Selects faces for refinement and coarsening using the Dörfler criterion.
 *
 * Given per-face residual values, this function partitions the set of faces into
 * two groups:
 * - Faces marked for refinement, whose cumulative residual ≥ θ_refine ⋅ (∑ residuals),
 *   subject to threshold_refine.
 * - Faces marked for coarsening, whose cumulative residual ≥ θ_coarse ⋅ (∑ residuals),
 *   subject to threshold_coarse.
 *
 * @param mesh      Input manifold surface mesh
 * @param residual  Per-face residual indicators
 * @param conf      Dörfler criterion configuration
 *
 * @return { refine_faces, coarse_faces } as an array of two vectors of faces
 */
std::array<std::vector<Face>, 2> select_doerfler(ManifoldSurfaceMesh &mesh, FaceData<double> residual, const DoeflerConf& conf);

}

