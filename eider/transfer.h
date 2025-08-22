#pragma once

#include <geometrycentral/surface/intrinsic_triangulation.h>
#include <vector>

namespace geometrycentral::surface {

// We make use of the fact that even if vertices are dead, they still store the vertexData.
// This now allows us to map different triangulations to their indices
struct AdaptiveTriplet{
    Vertex i;
    Vertex j;
    double value;
    AdaptiveTriplet(const Vertex &i, const Vertex &j, double value) : i(i), j(j), value(value) {}
    Eigen::Triplet<double> toEigen(const VertexData<std::size_t>& idx_i, const VertexData<std::size_t>& idx_j) const{
        return Eigen::Triplet<double>(int(idx_i[i]),int(idx_j[j]),value);
    }
};

/// Class to facilitate Transfering Vertexdata between adaptive triangulations
/// The Pipeline is -- base T -> Refined T -> coarse T --
/// We make use of the fact that the Refined T is the Common Subdivison of the Base T and coarse T.
/// We then use the L2 Projection to find the nearest function in the Garlekin view.
/// To do this, this class builds the required matrices and indices while the Mesh is adapted
class AdaptiveTransfer {
  public:
    /// Initialize Adaptive Transfer with Vertexdata at base Triangulation
    AdaptiveTransfer(IntrinsicTriangulation& tri, VertexData<double> f_B);

    void startRefine();

    void refineEdge(Vertex vi, Vertex vj, Vertex vp);

    void endRefine();

    void startCoarse();

    void coarseEdge(Vertex vi, Vertex vj, Vertex vp);
    void endCoarse();

    VertexData<double> transfer();
  private:
    ManifoldSurfaceMesh& mesh; IntrinsicGeometryInterface& geom;
    VertexData<std::size_t> base_Idx, refined_idx, coarse_idx;
    // Interpolation Matrixes from base/coarse T to Refined T
    std::vector<AdaptiveTriplet> triplets_B, triplets_C;
    // Garlekin Mass Matrix
    SparseMatrix<double> GLMM, P_B, P_C;

    // Amount of base/refined/coarse Vertices
    std::size_t nB = 0, nR = 0, nC = 0;

    Vector<double> vecfB;
};

}