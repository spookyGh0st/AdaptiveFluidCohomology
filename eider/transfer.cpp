#include "transfer.h"
#include <geometrycentral/numerical/linear_solvers.h>
#include <cassert>

namespace geometrycentral::surface {

AdaptiveVertexTransfer::AdaptiveVertexTransfer(IntrinsicTriangulation &tri, VertexData<double> &f_B)
    : AdaptiveTransfer(tri, f_B) { }

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
    AdaptiveTransfer::endRefine();
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
    AdaptiveTransfer::endCoarse();
}

void AdaptiveFaceTransfer::setSplitHeVec(Halfedge he, HalfedgeData<Vector2>& hev) {
    Vector2 v_kp= hev[he]*0.5 - (hev[he] + hev[he.next()]);
    splitHeVec[he.face()] = v_kp.normalize();
}

void AdaptiveFaceTransfer::refineEdge(Halfedge h_pj) {
    if (!h_pj.isInterior()) return;
    HalfedgeData<Vector2> &hev = geom.halfedgeVectorsInFace;
    Halfedge h_ip = h_pj.prevOrbitFace().twin().prevOrbitFace();
    Halfedge h_kp = h_pj.prevOrbitFace(), h_pk = h_ip.next();
    assert(h_ip.face() == h_pk.face());
    assert(h_pj.face() == h_kp.face());
    Face f1 = h_ip.face(),f2 = h_pj.face();
    assert(splitHeVec[f2] == Vector2::infinity());
    Vector2 v_kp = splitHeVec[f1], v_pk = - v_pk;
    assert(v_kp != Vector2::infinity());
    Vector2 e1 = hev[f1.halfedge()].normalize();
    Vector2 e2 = hev[f2.halfedge()].normalize();
    assert(e2.y == 0);
    Vector2 r1 = e1/v_pk, r2 = e2/v_kp;
    r[f1] *= r1; r[f2] *= r2;

    setSplitHeVec(h_pj.next(),hev);
    setSplitHeVec(h_ip.prevOrbitFace(),hev);
}


void AdaptiveFaceTransfer::refineEdge(Vertex vi, Vertex vj, Vertex vp) {
    HalfedgeData<Vector2> &hev = geom.halfedgeVectorsInFace;
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

}
