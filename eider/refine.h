#pragma once
#include "transfer.h"

#include <geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h>
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
    AdaptiveTriangulation(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface& geom);
    Halfedge vertex_bisection(Halfedge he, AdaptiveTransfer *transfer = nullptr);
    void refine(std::vector<Face> faces, AdaptiveTransfer *transfer = nullptr);

    // get halfedge that was used in the edgesplit resulting in vertex v
    Halfedge coarse_halfedge(Vertex v);
    Halfedge vertex_biunion(Halfedge he, AdaptiveTransfer *transfer = nullptr);
    void coarse(const std::vector<Face> &f, AdaptiveTransfer *transfer = nullptr);

    [[nodiscard]] const ManifoldSurfaceMesh &mesh() const { return *tri.intrinsicMesh; }
    [[nodiscard]] const IntrinsicGeometryInterface &geom() const { return tri; }
    [[nodiscard]] const IntrinsicTriangulation &intrinsicTriangulation() const { return tri; }
    ManifoldSurfaceMesh &mesh() { return *tri.intrinsicMesh; }
    IntrinsicGeometryInterface &geom() { return tri; }
    IntegerCoordinatesIntrinsicTriangulation &intrinsicTriangulation() { return tri; }

    Halfedge getRefinementEdge(Face f);
  private:
    IntegerCoordinatesIntrinsicTriangulation tri;
    void setRefinementEdge(Halfedge he);
  public:
    IncrementingIndex idx;
    CornerData<bool> marked_corner;
    Eigen::SparseMatrix<double> M_restriction;
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
 * The Dörfler criterion selects a subset of faces such that the
 * cumulative residual contribution exceeds a prescribed fraction θ of the
 * total residual, i.e. the minimal set of elements M_r ⊆ T_h such that
 * ∑_{T ∈ M_r} η_T²  ≥  θ ⋅ ∑_{T ∈ T_h} η_T², where
 * η_T : local error indicator for element T
 * θ   : bulk parameter, with 0 < θ < 1
 * T_h : current mesh
 *
 * Similarly, coarsed elements are the subsets of faces, such that the
 * cumulative residual contributions is smaller then a fraction of the
 * total residual.
 *
 * Both Inversing and coarsing additionally respect a absolut threshhold.
 *
 * - θ_coarse ∈ (0,1]  : fraction of residual to mark for coarsing
 * - θ_refine ∈ (0,1]  : fraction of residual to mark for refinement
 * - threshold_refine  : absolute lower bound for refinement marking,
 *                  i.e. only elements with higher residual are refined
 * - threshold_coarse  : absolute upper bound for coarsening marking
 *                  i.e. only elements with lower residual are coarsed
 */
struct DoeflerConf {
    double theta_coarse = 0.1;                                              ///< θ_coarse: fraction of residual for coarsening
    double theta_refine = 0.1;                                              ///< θ_refine: fraction of residual for refinement
    double threshold_refine = 0.000001;                                     ///< refinement cutoff
    double threshold_coarse = 0.000001;           ///< coarsening cutoff
};

enum class DoerflerPresetConf{ LOW, MEDIUM, HIGH,VERY_HIGH, UNIFORM_REFINE, UNIFORM_COARSE };
DoeflerConf DoerflerPreset(DoerflerPresetConf preset);

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

