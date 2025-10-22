#pragma once

#include <geometrycentral/surface/intrinsic_triangulation.h>
#include <vector>

namespace geometrycentral::surface {

// We make use of the fact that even if vertices are dead, they still store the vertexData.
// This now allows us to map different triangulations to their indices
template <typename E>
struct AdaptiveTriplet{
    E i;
    E j;
    double value;
    AdaptiveTriplet(const E &i, const E &j, double value) : i(i), j(j), value(value) {}
    Eigen::Triplet<double,Eigen::Index> toEigen(const MeshData<E,std::size_t>& idx_i, const MeshData<E,std::size_t>& idx_j) const{
        return Eigen::Triplet(Eigen::Index(idx_i[i]),Eigen::Index(idx_j[j]),value);
    }
};

struct TransferSolver {
    TransferSolver(SparseMatrix<double> P_B, SparseMatrix<double> P_A, SparseMatrix<double> M)
    : solver(P_B.transpose() * M * P_B), PBT_M_PA(P_B.transpose()*M*P_A) {}
    Eigen::ConjugateGradient<SparseMatrix<double>> solver;
    SparseMatrix<double>PBT_M_PA;
    Vector<double> solveWithGuess(const Vector<double> &fA, const Vector<double> &guess) const {
        return solver.solveWithGuess(PBT_M_PA*fA,guess);
    }
};


/// Class to facilitate Transfering Vertexdata between adaptive triangulations
/// The Pipeline is -- base T -> Refined T -> coarse T --
/// We make use of the fact that the Refined T is the Common Subdivison of the Base T and coarse T.
/// We then use the L2 Projection to find the nearest function in the Garlekin view.
/// To do this, this class builds the required matrices and indices while the Mesh is adapted
template <typename E>
class AdaptiveTransfer {
    using T = double;
public:
    AdaptiveTransfer(IntrinsicTriangulation& Tri,MeshData<E,T>& f_B);
    virtual ~AdaptiveTransfer() = default;
    virtual void startRefine() = 0;
    virtual void refineEdge(Vertex vi, Vertex vj, Vertex vp) = 0;
    virtual void endRefine();
    virtual void startCoarse()=0;
    virtual void coarseEdge(Vertex vi, Vertex vj, Vertex vp) =0;
    virtual void endCoarse();

    MeshData<E,double> transfer() const;
protected:
    // function at base

    ManifoldSurfaceMesh& mesh; IntrinsicGeometryInterface& geom;
    MeshData<E,std::size_t> base_Idx, refined_idx, coarse_idx;
    // Interpolation Matrixes from base/coarse T to Refined T
    std::vector<AdaptiveTriplet<E>> triplets_B, triplets_C;
    // Garlekin Mass Matrix, permutation matrices
    SparseMatrix<double> GLMM, P_B, P_C;

    // Amount of base/refined/coarse elements
    std::size_t nB = 0, nR = 0, nC = 0;

    MeshData<E,T>& f_B;
    Vector<double> vecfB;
};


template <typename E> MeshData<E, double> scalar_transfer(
    const TransferSolver &solver,
    MeshData<E, std::size_t> coarse_indices,
    MeshData<E, double> &f_base,
    const Vector<double> &f_base_vec)
{
    Vector<double> f_coarse_vec = solver.solveWithGuess(f_base_vec, f_base.toVector(coarse_indices));
    assert(f_coarse_vec.allFinite());
    return MeshData<E,double>(coarse_indices.mesh, f_coarse_vec, coarse_indices);
}

// elementwise scalar transfer
// Note that since tangent spaces are not equivalent,
// this is not a transfer of vectors
template <typename E> MeshData<E, Vector2> array_transfer(
    ManifoldSurfaceMesh &mesh,
    const TransferSolver &solver,
    MeshData<E, std::size_t> coarse_indices,
    MeshData<E, Vector2> &f_base,
    const Vector<Vector2> &f_base_vec) {
    Vector<Vector2> f_coarse_vec(nElements<E>(mesh), Vector2::zero());

    // Element wise L2 Projection
    for (int i = 0; i < 2; ++i) {
        Eigen::VectorXd f_base_vec_i(f_base_vec.size());
        for (size_t j = 0; i < f_base_vec.size(); ++i) { f_base_vec_i[j] = f_base_vec[j][i]; }
        Vector<double> f_coarse_vec_i = solver.solveWithGuess(f_base_vec, f_base.toVector(coarse_indices));
        assert(f_coarse_vec_i.allFinite());
        for (size_t j = 0; i < f_base_vec.size(); ++i) { f_coarse_vec[j][i] = f_coarse_vec_i[j]; }
    }

    return MeshData<E, Vector2>(mesh, f_coarse_vec, coarse_indices);
}


template <typename E>
MeshData<E,std::size_t> getElementIndices(ManifoldSurfaceMesh& mesh);

template <> inline MeshData<Vertex,std::size_t> getElementIndices<Vertex>(ManifoldSurfaceMesh& mesh) {
    return  mesh.getVertexIndices();
}
template <> inline MeshData<Face,std::size_t> getElementIndices<Face>(ManifoldSurfaceMesh& mesh) {
    return  mesh.getFaceIndices();
}

template <typename E> std::size_t nElements(ManifoldSurfaceMesh& mesh);
template <> inline std::size_t nElements<Vertex>(ManifoldSurfaceMesh& mesh) {
    return  mesh.nVertices();
}
template <> inline std::size_t nElements<Face>(ManifoldSurfaceMesh& mesh) {
    return  mesh.nFaces();
}

template <typename E>
AdaptiveTransfer<E>::AdaptiveTransfer(IntrinsicTriangulation &tri, MeshData<E, double> &f_B)
    : mesh(*tri.intrinsicMesh), geom(tri),
    base_Idx(getElementIndices<E>(mesh)), nB(nElements<E>(mesh)), vecfB(f_B.toVector(base_Idx)), f_B(f_B)
{
    triplets_B.reserve(nB*2);
    triplets_C.reserve(nB*2);
}

template <typename E>
void AdaptiveTransfer<E>::endRefine() {
    nR = nElements<E>(mesh);

    geom.refreshQuantities();
    refined_idx = getElementIndices<E>(mesh);

    std::vector<Eigen::Triplet<double,Eigen::Index>> triplets;
    triplets.reserve(triplets_B.size());
    for (const AdaptiveTriplet<E>& at: triplets_B) {
        triplets.emplace_back(at.toEigen(refined_idx,base_Idx));
    }
    triplets_B.clear(); // clear elements
    P_B = SparseMatrix<double>(nR,nB);
    P_B.setFromTriplets(triplets.begin(),triplets.end());
}

template <typename E>
void AdaptiveTransfer<E>::endCoarse() {
    coarse_idx = getElementIndices<E>(mesh);
    nC = nElements<E>(mesh);

    std::vector<Eigen::Triplet<double,Eigen::Index>> triplets;
    triplets.reserve(triplets_C.size());
    for (const AdaptiveTriplet<E>& at: triplets_C) {
        triplets.emplace_back(at.toEigen(refined_idx,coarse_idx));
    }
    triplets_C.clear(); // clear elements
    P_C = SparseMatrix<double>(nR,nC);
    P_C.setFromTriplets(triplets.begin(),triplets.end());

}

template <typename E>
MeshData<E,double> AdaptiveTransfer<E>::transfer() const {
    assert(P_B.rows() == nR && P_B.cols() == nB);
    assert(P_B.cols() == vecfB.size());
    assert(P_C.rows() == nR && P_C.cols() == nC);
    assert(GLMM.cols() == nR);
    assert(GLMM.rows() == nR);
    auto solver = TransferSolver(P_C,P_B,GLMM);
    Vector<double> vec_fC = solver.solveWithGuess(vecfB,f_B.toVector(coarse_idx));
    return MeshData<E,double>(mesh,vec_fC,coarse_idx);
}


class AdaptiveVertexTransfer: public AdaptiveTransfer<Vertex>{
  public:
    /// Initialize Adaptive Transfer with Vertexdata at base Triangulation
    AdaptiveVertexTransfer(IntrinsicTriangulation& tri, VertexData<double>& f_B);

    void startRefine() override;

    void refineEdge(Vertex vi, Vertex vj, Vertex vp) override;

    void endRefine() override;

    void startCoarse() override;

    void coarseEdge(Vertex vi, Vertex vj, Vertex vp) override;
    void endCoarse() override;
};

class AdaptiveFaceTransfer: public AdaptiveTransfer<Face>{
  public:
    /// Initialize Adaptive Transfer with Vertexdata at base Triangulation
    AdaptiveFaceTransfer(IntrinsicTriangulation& tri);

    void startRefine() override;

    void refineEdge(Vertex vi, Vertex vj, Vertex vp) override;

    void endRefine() override;

    void startCoarse() override;

    void coarseEdge(Vertex vi, Vertex vj, Vertex vp) override;
    void endCoarse() override;
private:
    void setSplitHeVec(Halfedge he, HalfedgeData<Vector2>& hev);
    void refineEdge(Halfedge hpjd);
    CornerData<bool>* markedCorner;
    FaceData<Vector2> r;
    FaceData<Vector2> splitHeVec;
};
}