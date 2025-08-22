#include "transfer.h"
#include <geometrycentral/numerical/linear_solvers.h>

namespace geometrycentral::surface {

void AdaptiveTransfer::startRefine() {
    for (Vertex v: mesh.vertices()) {
        triplets_B.emplace_back(v,v,1);
    }
}
void AdaptiveTransfer::refineEdge(Vertex vi, Vertex vj, Vertex vp) {
    triplets_B.emplace_back(vp,vi,0.5);
    triplets_B.emplace_back(vp,vj,0.5);
}
void AdaptiveTransfer::endRefine() {
    refined_idx = mesh.getVertexIndices();
    nR = mesh.nVertices();

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(triplets_B.size());
    for (const AdaptiveTriplet& at: triplets_B) {
        triplets.emplace_back(at.toEigen(refined_idx,base_Idx));
    }
    triplets_B.clear(); // clear elements
    P_B = SparseMatrix<double>(nR,nB);
    P_B.setFromTriplets(triplets.begin(),triplets.end());
    geom.requireVertexGalerkinMassMatrix();
    GLMM = geom.vertexGalerkinMassMatrix;
    geom.unrequireVertexGalerkinMassMatrix();
}
void AdaptiveTransfer::startCoarse() {
    if(nR != 0){
        // Assert we did not change any elements since end of refinement
        assert(refined_idx.raw() == mesh.getVertexIndices().raw());
        assert(nR == mesh.nVertices());
    } else {
        refined_idx = mesh.getVertexIndices();
        nR = mesh.nVertices();
    }
}
void AdaptiveTransfer::coarseEdge(Vertex vi, Vertex vj, Vertex vp) {
    triplets_C.emplace_back(vp,vi,0.5);
    triplets_C.emplace_back(vp,vj,0.5);
}
void AdaptiveTransfer::endCoarse() {
    coarse_idx = mesh.getVertexIndices();
    nC = mesh.nVertices();

    for (Vertex v: mesh.vertices()) {
        triplets_C.emplace_back(v,v,1);
    }
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(triplets_C.size());
    for (const AdaptiveTriplet& at: triplets_C) {
        triplets.emplace_back(at.toEigen(refined_idx,coarse_idx));
    }
    triplets_C.clear(); // clear elements
    P_C = SparseMatrix<double>(nR,nC);
    P_C.setFromTriplets(triplets.begin(),triplets.end());
}
AdaptiveTransfer::AdaptiveTransfer(IntrinsicTriangulation &tri, VertexData<double> f_B)
    :mesh(*tri.intrinsicMesh), geom(tri), base_Idx ( mesh.getVertexIndices()), nB(mesh.nVertices()), vecfB(f_B.toVector(base_Idx)) {

}

VertexData<double> AdaptiveTransfer::transfer() {
    assert(P_B.rows() == nR && P_B.cols() == nB);
    assert(P_B.cols() == vecfB.size());
    assert(P_C.rows() == nR && P_C.cols() == nC);
    assert(GLMM.cols() == nR);
    assert(GLMM.rows() == nR);

    SparseMatrix<double> mat = P_C.transpose() * GLMM * P_C;
    auto AtoB_L2_Solver_F = std::make_unique<SquareSolver<double>>(mat);
    Vector<double> vec = P_C.transpose() * GLMM * P_B * vecfB;
    Vector<double> vec_fC = AtoB_L2_Solver_F->solve(vec);
    return VertexData<double>(mesh,vec_fC,coarse_idx);
}

}
