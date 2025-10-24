#pragma once

#include <geometrycentral/surface/intrinsic_triangulation.h>
#include <vector>
#include <optional>

namespace geometrycentral::surface {

// We make use of the fact that even if vertices are dead, they still store the vertexData.
// This now allows us to map different triangulations to their indices
template <typename E, typename T>
struct AdaptiveTriplet{
    E i;
    E j;
    T value;
    AdaptiveTriplet(const E &i, const E &j, T value) : i(i), j(j), value(value) {}
    Eigen::Triplet<T,Eigen::Index> toEigen(const MeshData<E,std::size_t>& idx_i, const MeshData<E,std::size_t>& idx_j) const{
        return Eigen::Triplet(Eigen::Index(idx_i[i]),Eigen::Index(idx_j[j]),value);
    }
};

template <typename T>
struct TransferSolver {
    TransferSolver(SparseMatrix<T> P_B, SparseMatrix<T> P_A, SparseMatrix<T> M)
        : solver(P_B.adjoint() * M * P_B), PBT_M_PA(P_B.adjoint()*M*P_A) { }
    Eigen::ConjugateGradient<SparseMatrix<T>> solver;
    SparseMatrix<T>PBT_M_PA;
    Vector<T> solveWithGuess(const Vector<T> &fA, const Vector<T> &guess) const {
        return solver.solveWithGuess(PBT_M_PA*fA,guess);
    }
};

class AdaptiveTransfer {
   public:
     virtual void startRefine() = 0;
     virtual void refineEdge(Vertex vi, Vertex vj, Vertex vp) = 0;
     virtual void endRefine() = 0;
     virtual void startCoarse()=0;
     virtual void coarseEdge(Vertex vi, Vertex vj, Vertex vp) =0;
     virtual void endCoarse() = 0;
};

struct Tri {
    // Ordered such that 0 is edge to be split/joined
    Face f;
    std::array<Halfedge,3> halfedges;
    std::array<Vector2,3> vectors;
    int8_t idx_f_he, idx_split_he;

    Halfedge face_he() const { return halfedges[idx_f_he]; }
    Halfedge edge_he() const { return halfedges[idx_split_he]; }
    Vector2 e1() const { return vectors[idx_f_he].normalize(); }
    bool isDead() const { return f.isDead(); }

    Tri(Halfedge he, const HalfedgeData<Vector2>& he_tang)
        : f(he.face()),
          halfedges({he, he.next(), he.next().next()}),
          idx_split_he(0)
    {
        int8_t i = 0;
        for (Halfedge he : f.adjacentHalfedges()) {
            vectors[i] = he_tang[he];
            if (he == f.halfedge()) idx_f_he = i;
            i++;
        }
    }

    Tri() = default;
};

struct Side {
    std::array<Tri, 2> tris;
    Side() = default;

    Side(Halfedge pj, const HalfedgeData<Vector2> &he_tang)
        : tris({
              Tri(pj.prevOrbitFace().twin().prevOrbitFace(), he_tang),
              Tri(pj, he_tang),
          })
    {}
};

struct Diamond {
    Vertex vi, vj, vp;
    std::array<std::optional<Side>, 2> sides;

    Diamond(Halfedge pj, const HalfedgeData<Vector2> &he_tang) {
        vi = pj.prevOrbitFace().twin().prevOrbitFace().vertex();
        vj = pj.tipVertex();
        vp = pj.tailVertex();

        sides[0] = Side(pj, he_tang);
        if (pj.twin().isInterior())
            sides[1] = Side(pj.twin().next().twin().next(), he_tang);
    }
};

struct Quad {
    Vertex vi, vj;
    std::array<std::optional<Tri>,2> tris;
    std::array<Vector2,2> v_kp;

    Quad(Halfedge ij, const HalfedgeData<Vector2> &he_tang) {
        v_kp[0] = (he_tang[ij] * 0.5 - (he_tang[ij] + he_tang[ij.next()])).normalize();
        vi = ij.tailVertex();
        vj = ij.tipVertex();
        tris[0] = Tri(ij, he_tang);

        if (ij.twin().isInterior()) {
            Halfedge ji = ij.twin();
            tris[1] = Tri(ij.twin(), he_tang);
            v_kp[1] = (he_tang[ji] * 0.5 - (he_tang[ji] + he_tang[ji.next()])).normalize();
        }
    }
};


/// Class to facilitate Transfering Vertexdata between adaptive triangulations
/// The Pipeline is -- base T -> Refined T -> coarse T --
/// We make use of the fact that the Refined T is the Common Subdivison of the Base T and coarse T.
/// We then use the L2 Projection to find the nearest function in the Garlekin view.
/// To do this, this class builds the required matrices and indices while the Mesh is adapted
template <typename E,typename T>
class AdaptiveTransferL2: public AdaptiveTransfer{
public:
  AdaptiveTransferL2(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom);
    virtual ~AdaptiveTransferL2() = default;
protected:
    // function at base

    ManifoldSurfaceMesh& mesh; IntrinsicGeometryInterface& geom;
    MeshData<E,std::size_t> base_Idx, refined_idx, coarse_idx;
    // Interpolation Matrixes from base/coarse T to Refined T
    std::vector<AdaptiveTriplet<E,T>> triplets_B, triplets_C;
    // Garlekin Mass Matrix, permutation matrices
    SparseMatrix<T> GLMM, P_B, P_C;

    // Amount of base/refined/coarse elements
    std::size_t nB = 0, nR = 0, nC = 0;

};


class AdaptiveVertexTransfer: public AdaptiveTransferL2<Vertex,double>{
  public:
    /// Initialize Adaptive Transfer with Vertexdata at base Triangulation
    AdaptiveVertexTransfer(IntrinsicTriangulation& tri, VertexData<double>& f_B);

    void startRefine() override;

    void refineEdge(Vertex vi, Vertex vj, Vertex vp) override;

    void endRefine() override;

    void startCoarse() override;

    void coarseEdge(Vertex vi, Vertex vj, Vertex vp) override;
    void endCoarse() override;

    VertexData<double> transfer() const;
  private:
    VertexData<double>& f_B;
    Vector<double> vecfB;
};

class AdaptiveFaceTransfer: public AdaptiveTransferL2<Face,std::complex<double>>{
  public:
    /// Initialize Adaptive Transfer with Vertexdata at base Triangulation
    AdaptiveFaceTransfer(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom, FaceData<Vector2> &f_B, CornerData<bool> &marked_corners);

    void startRefine() override;
    void refineEdge(Vertex vi, Vertex vj, Vertex vp) override;
    void refineEdge(const Quad &T0, const Diamond &T1);
    void endRefine() override;
    void startCoarse() override;
    void coarseEdge(Vertex vi, Vertex vj, Vertex vp) override;
    void coarseEdge(const Quad &T0, const Diamond &T1);
    void endCoarse() override;

    FaceData<Vector2> transfer() const;
private:
    using complex_t = std::complex<double>;
    SparseMatrix<complex_t> M_CS_Lumped();
    void setSplitHeVec(Halfedge he, HalfedgeData<Vector2>& hev);
    void refineEdge(Halfedge hpjd);
    FaceData<Vector2> r;
    FaceData<Vector2> splitHeVec;
    FaceData<complex_t> fB_complex;

    FaceData<Vector2>& f_B;
    FaceData<complex_t> f_B_complex;
    Vector<complex_t> vecfB;
    void refineSide(const Tri &t, const Side &s, Vector2 v_kp);
    void coarseSide(const Tri &t, const Side &s, Vector2 v_kp);
};
}