#include "transfer.h"

#include "truthy_iterator.h"

#include <geometrycentral/numerical/linear_solvers.h>
#include <cassert>

namespace geometrycentral::surface {

/************************************************************
 *                    L2 Transfer
 ************************************************************/

template <typename E>
MeshData<E,std::size_t> getElementIndices(ManifoldSurfaceMesh& mesh);

template <> inline MeshData<Vertex,std::size_t> getElementIndices<Vertex>(ManifoldSurfaceMesh& mesh) {
    return  mesh.getVertexIndices();
}
template <> inline MeshData<Face,std::size_t> getElementIndices<Face>(ManifoldSurfaceMesh& mesh) {
    return  mesh.getFaceIndices();
}

template <typename E> std::size_t nElements(SurfaceMesh& mesh);
template <> inline std::size_t nElements<Vertex>(SurfaceMesh& mesh) {
    return  mesh.nVertices();
}
template <> inline std::size_t nElements<Face>(SurfaceMesh& mesh) {
    return  mesh.nFaces();
}

template <typename E, typename T>
AdaptiveTransferL2<E,T>::AdaptiveTransferL2(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom)
    : mesh(mesh), geom(geom),
      base_Idx(getElementIndices<E>(mesh)), nB(nElements<E>(mesh))
{
    triplets_B.reserve(nB*2);
    triplets_C.reserve(nB*2);
}


template <typename E, typename T>
void AdaptiveTransferL2<E,T>::endCoarse() {
    coarse_idx = getElementIndices<E>(mesh);
    nC = nElements<E>(mesh);

    std::vector<Eigen::Triplet<T,Eigen::Index>> triplets;
    triplets.reserve(triplets_C.size());
    for (const AdaptiveTriplet<E,T>& at: triplets_C) {
        triplets.emplace_back(at.toEigen(refined_idx,coarse_idx));
    }
    triplets_C.clear(); // clear elements
    P_C = SparseMatrix<T>(nR,nC);
    P_C.setFromTriplets(triplets.begin(),triplets.end());

}

VertexData<double> AdaptiveVertexTransfer::transfer() const {
    assert(P_B.rows() == nR && P_B.cols() == nB);
    assert(P_B.cols() == vecfB.size());
    assert(P_C.rows() == nR && P_C.cols() == nC);
    assert(GLMM.cols() == nR);
    assert(GLMM.rows() == nR);
    auto solver = TransferSolver(P_C,P_B,GLMM);
    Vector<double> vec_fC = solver.solveWithGuess(vecfB,f_B.toVector(coarse_idx));
    return VertexData<double>(mesh,vec_fC,coarse_idx);
}

template class AdaptiveTransferL2<Vertex,double>;
template class AdaptiveTransferL2<Face,Vector2>;



/************************************************************
 *                    Vertex L2 Transfer
 ************************************************************/

// TODO: change to mesh, geom
AdaptiveVertexTransfer::AdaptiveVertexTransfer(IntrinsicTriangulation &tri, VertexData<double> &f_B)
    : AdaptiveTransferL2(*tri.intrinsicMesh, tri), f_B(f_B), vecfB(f_B.toVector(base_Idx)) { }

void AdaptiveVertexTransfer::startRefine() {
    for (Vertex v: mesh.vertices()) {
        triplets_B.emplace_back(v,v,1);
    }
}
void AdaptiveVertexTransfer::refineEdge(Vertex vi, Vertex vj, Vertex vp) {
    triplets_B.emplace_back(vp,vi,0.5);
    triplets_B.emplace_back(vp,vj,0.5);
}
void AdaptiveVertexTransfer::endRefine() {
    nR = mesh.nVertices();
    geom.refreshQuantities();
    refined_idx = getElementIndices<Vertex>(mesh);

    std::vector<Eigen::Triplet<double,Eigen::Index>> triplets;
    triplets.reserve(triplets_B.size());
    for (const AdaptiveTriplet<Vertex,double>& at: triplets_B) {
        triplets.emplace_back(at.toEigen(refined_idx,base_Idx));
    }
    triplets_B.clear(); // clear elements
    P_B = SparseMatrix<double>(nR,nB);
    P_B.setFromTriplets(triplets.begin(),triplets.end());

    geom.requireVertexGalerkinMassMatrix();
    GLMM = geom.vertexGalerkinMassMatrix;
    geom.unrequireVertexGalerkinMassMatrix();
    assert(GLMM.rows() == nR);
    assert(GLMM.cols() == nR);
}

void AdaptiveVertexTransfer::startCoarse() {
    if(nR != 0){
        // Assert we did not change any elements since end of refinement
        assert(refined_idx.raw() == mesh.getVertexIndices().raw());
        assert(nR == mesh.nVertices());
    } else {
        refined_idx = mesh.getVertexIndices();
        nR = mesh.nVertices();
    }
}
void AdaptiveVertexTransfer::coarseEdge(Vertex vi, Vertex vj, Vertex vp) {
    triplets_C.emplace_back(vp,vi,0.5);
    triplets_C.emplace_back(vp,vj,0.5);
}
void AdaptiveVertexTransfer::endCoarse() {
    for (Vertex v: mesh.vertices()) {
        triplets_C.emplace_back(v,v,1);
    }

    coarse_idx = getElementIndices<Vertex>(mesh);
    nC = nElements<Vertex>(mesh);

    std::vector<Eigen::Triplet<double,Eigen::Index>> triplets;
    triplets.reserve(triplets_C.size());
    for (const AdaptiveTriplet<Vertex,double>& at: triplets_C) {
        triplets.emplace_back(at.toEigen(refined_idx,coarse_idx));
    }
    triplets_C.clear(); // clear elements
    P_C = SparseMatrix<double>(nR,nC);
    P_C.setFromTriplets(triplets.begin(),triplets.end());
}

/************************************************************
 *                    Adaptive Face Transfer
 ************************************************************/
template <typename E> MeshData<E,std::complex<double>> toComplex(MeshData<E,Vector2> d){
    MeshData<E,std::complex<double>> result(*d.getMesh());
    for (E e : iterateElements<E>(d.getMesh())) { result[e] = d[e]; }
    return result;
}
template <typename E> MeshData<E,Vector2> toVector(MeshData<E,std::complex<double>> d){
    MeshData<E,Vector2> result(*d.getMesh());
    for (E e : iterateElements<E>(d.getMesh())) { result[e] = Vector2::fromComplex(d[e]); }
    return result;
}
template <typename E> Eigen::VectorXcd toVectorC(MeshData<E,Vector2> d, MeshData<E,std::size_t> indexer){
    using complex_t = std::complex<double>;
    Eigen::VectorXcd result(nElements<E>(*d.getMesh()));
    for (E e : iterateElements<E>(d.getMesh())) {
        if (indexer[e] != std::numeric_limits<size_t>::max()) {
            result(indexer[e]) = complex_t(d[e].x,d[e].y);
        }
    }
    return result;
}




AdaptiveFaceTransfer::AdaptiveFaceTransfer(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom, FaceData<Vector2> &f_B, CornerData<bool> &marked_corners)
    : AdaptiveTransferL2<Face, complex_t>(mesh, geom),
      f_B(f_B), f_B_complex(toComplex(f_B)), vecfB(f_B_complex.toVector(base_Idx)),
      r(mesh,Vector2::fromAngle(0)),
      splitHeVec(mesh,Vector2::infinity())
{
    geom.requireHalfedgeVectorsInFace();
    for (Face f: mesh.faces()) {
        for(Halfedge he: f.adjacentHalfedges()){
            if(marked_corners[he.oppositeCorner()])  setSplitHeVec(he,geom.halfedgeVectorsInFace);
        }
    }
}

void AdaptiveFaceTransfer::setSplitHeVec(Halfedge he, HalfedgeData<Vector2>& hev) {
    Vector2 v_kp= hev[he]*0.5 - (hev[he] + hev[he.next()]);
    splitHeVec[he.face()] = v_kp.normalize();
}

void AdaptiveFaceTransfer::refineEdge(Quad& T0, Diamond& T1) {
    // TODO:
}


void AdaptiveFaceTransfer::refineEdge(Halfedge h_pj) {
    if (!h_pj.isInterior()) return;
    HalfedgeData<Vector2> &hev = geom.halfedgeVectorsInFace;
    Halfedge h_ip = h_pj.prevOrbitFace().twin().prevOrbitFace();
    Halfedge h_kp = h_pj.prevOrbitFace(), h_pk = h_ip.next();
    assert(h_ip.face() == h_pk.face());
    assert(h_pj.face() == h_kp.face());
    Face f1 = h_ip.face(),f2 = h_pj.face();
    assert(f1.halfedge() == h_pk);
    assert(f2.halfedge() == h_kp);
    assert(splitHeVec[f2] == Vector2::infinity());
    Vector2 v_kp = splitHeVec[f1], v_pk = - v_kp;
    assert(v_kp != Vector2::infinity());
    Vector2 e1 = hev[f1.halfedge()].normalize();
    Vector2 e2 = hev[f2.halfedge()].normalize();
    assert(e2.y == 0);
    Vector2 r1 = e1/v_pk, r2 = e2/v_kp;
    r[f1] *= r1; r[f2] *= r2;

    triplets_B.emplace_back(f1,f1,r1);
    triplets_B.emplace_back(f2,f1,r2);

    setSplitHeVec(h_pj.next(),hev);
    setSplitHeVec(h_ip.prevOrbitFace(),hev);
}


void AdaptiveFaceTransfer::refineEdge(Vertex vi, Vertex vj, Vertex vp) {
    Halfedge hi, hj;
    for (Halfedge he: vp.outgoingHalfedges()) {
        if (he.tipVertex() == vi) hi = he;
        if (he.tipVertex() == vj) hj = he;
    }

    assert(hi != Halfedge() && hj != Halfedge());
    for (Halfedge he: {hi,hj}) {
        refineEdge(he);
    }
}


void AdaptiveFaceTransfer::startRefine() {
    for (Face f: mesh.faces()) {
        triplets_B.emplace_back(f,f,complex_t(1,0));
    }
}
void AdaptiveFaceTransfer::endRefine() {
    nR = nElements<Face>(mesh);
    geom.refreshQuantities();
    refined_idx = getElementIndices<Face>(mesh);

    std::vector<Eigen::Triplet<complex_t,Eigen::Index>> triplets;
    triplets.reserve(triplets_B.size());
    for (const AdaptiveTriplet<Face,complex_t>& at: triplets_B) {
        triplets.emplace_back(at.toEigen(refined_idx,base_Idx));
    }
    triplets_B.clear(); // clear elements
    P_B = SparseMatrix<complex_t>(nR,nB);

    P_B.setFromTriplets(triplets.begin(),triplets.end(),[]( const complex_t& a, const complex_t& b) { return a * b; });


    GLMM = M_CS_Lumped();
    assert(GLMM.rows() == nR);
    assert(GLMM.cols() == nR);

}
void AdaptiveFaceTransfer::startCoarse() {
    if(nR != 0){
        assert(nR == mesh.nFaces());
    } else {
        refined_idx = mesh.getFaceIndices();
        nR = mesh.nFaces();
    }
}
void AdaptiveFaceTransfer::coarseEdge(Vertex vi, Vertex vj, Vertex vp) {
    assert(vp.isDead());
    Halfedge hij;
    for (Halfedge he: vi.outgoingHalfedges()) { if (he.tipVertex()==vj) hij = he; }
    assert(hij!=Halfedge());

}
void AdaptiveFaceTransfer::endCoarse() {
    AdaptiveTransferL2::endCoarse();

    for (Face f: mesh.faces()) {
        triplets_C.emplace_back(f,f,complex_t(1,0));
    }

    coarse_idx = getElementIndices<Face>(mesh);
    nC = nElements<Face>(mesh);

    std::vector<Eigen::Triplet<complex_t,Eigen::Index>> triplets;
    triplets.reserve(triplets_C.size());
    for (const AdaptiveTriplet<Face,complex_t>& at: triplets_C) {
        triplets.emplace_back(at.toEigen(refined_idx,coarse_idx));
    }
    triplets_C.clear(); // clear elements
    P_C = SparseMatrix<complex_t>(nR,nC);
    P_C.setFromTriplets(triplets.begin(),triplets.end());
}

FaceData<Vector2> AdaptiveFaceTransfer::transfer() const {
    assert(P_B.rows() == nR && P_B.cols() == nB);
    assert(P_B.cols() == vecfB.size());
    assert(P_C.rows() == nR && P_C.cols() == nC);
    assert(GLMM.cols() == nR);
    assert(GLMM.rows() == nR);
    auto solver = TransferSolver(P_C,P_B,GLMM);
    Vector<complex_t> vec_fC = solver.solveWithGuess(vecfB,f_B_complex.toVector(coarse_idx));
    FaceData<complex_t> fC_complex(mesh,vec_fC,coarse_idx);
    return toVector(fC_complex);
}

SparseMatrix<AdaptiveFaceTransfer::complex_t> AdaptiveFaceTransfer::M_CS_Lumped() {
    std::vector<Eigen::Triplet<complex_t, Eigen::Index>> triplets;
    geom.requireFaceAreas();
    for (Face f : mesh.faces()) {
        int i = refined_idx[f];
        double v = geom.faceAreas[f];
        triplets.emplace_back(i, i, complex_t(v, 0));
    }
    geom.unrequireFaceAreas();
    SparseMatrix<complex_t> M(mesh.nFaces(), mesh.nFaces());
    M.setFromTriplets(triplets.begin(), triplets.end());
    return M;
}
}
