#include "homology.h"
#include <algorithm>
#include <bitset>
#include <numeric>
#include <tuple>
#include <vector>

#include "geometrycentral/surface/barycentric_vector.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"

namespace geometrycentral::surface {

std::vector<Singular_Circle> singular_homology_basis(ManifoldSurfaceMesh &mesh, std::vector<std::vector<Halfedge>> hom_basis) {
    std::vector<Singular_Circle> sing_basis(hom_basis.size());
    for (std::size_t i = 0; i < hom_basis.size(); ++i) {
        auto h = reduce_co_loop(mesh, hom_basis[i]);
        sing_basis[i] = Singular_Circle(mesh, h);
    }
    return sing_basis;
}

EdgeData<double> delta_form(ManifoldSurfaceMesh &mesh, const Singular_Circle &co_loop) {
    EdgeData<double> delta(mesh, 0);
    for (Halfedge he : co_loop) {
        delta[he.edge()] = he.orientation() ? 1 : -1;
    }
    return delta;
}

EdgeData<double> delta_form(ManifoldSurfaceMesh &mesh,
                            const std::vector<Halfedge> &co_loop) {
    EdgeData<double> delta(mesh, 0);
    for (Halfedge he : co_loop) {
        delta[he.edge()] = he.orientation() ? 1 : -1;
    }
    return delta;
}

void PressureProjectionSolver::compute(IntrinsicGeometryInterface &geom) {
    geom.requireDECOperators();
    A = geom.d0;
    AT = A.transpose();
    solver.compute(AT * A);
    // TODO: QR - solver.compute(A);
    geom.unrequireDECOperators();
}

EdgeData<double>
PressureProjectionSolver::solve(ManifoldSurfaceMesh &mesh,
                                const EdgeData<double> &co_loop) const {
    Eigen::VectorXd x = co_loop.toVector();
    // TODO: QR - Eigen::VectorXd c = solver.solve(x);
    Eigen::VectorXd c = solver.solve(AT * x);
    return EdgeData<double>(mesh, x - A * c);
}

// Function to compute the Whitney interpolated vector field
FaceData<Vector2> whitney_interpolation(ManifoldSurfaceMesh &mesh,
                                        IntrinsicGeometryInterface &geom,
                                        EdgeData<double> &h) {
    geom.requireHalfedgeVectorsInFace();
    geom.requireFaceAreas();
    FaceData<Vector2> vField(mesh, Vector2::zero());

    for (Face f : mesh.faces()) {
        Vector2 vT = Vector2::zero();
        double A = geom.faceAreas[f];
        assert(A > 0);
        for (Halfedge he : f.adjacentHalfedges()) {
            Edge e = he.edge();
            double wij = (he.orientation()) ? h[e] : -h[e];
            Vector2 ni = geom.halfedgeVectorsInFace[he.next()].rotate90(),
                    nj =
                        geom.halfedgeVectorsInFace[he.next().next()].rotate90();
            vT += wij / 6 * (nj - ni);
        }
        assert(vT.isFinite());
        vField[f] = vT;
    }

    return vField;
}

using VectorX2d = Eigen::Matrix<Vector2, -1, 1>;
using MatrixX2d = Eigen::Matrix<Vector2, -1, -1>;

using InnerProductFn =
    std::function<double(const VectorX2d &, const VectorX2d &)>;

void modifiedGramSchmidt(const MatrixX2d &A, MatrixX2d &Q, Eigen::MatrixXd &R, const InnerProductFn &innerProduct) {
    const int m = A.rows();
    const int n = A.cols();

    Q = MatrixX2d::Zero(m, n);
    R = Eigen::MatrixXd::Zero(n, n);

    for (int j = 0; j < n; ++j) {
        VectorX2d v = A.col(j);

        for (int i = 0; i < j; ++i) {
            R(i, j) = innerProduct(Q.col(i), A.col(j));
            VectorX2d w = Q.col(i);
            for (int k = 0; k < v.size(); ++k) {
                v[k] -= w[k] * R(i, j);
            }
        }

        R(j, j) = std::sqrt(innerProduct(v, v));
        if (R(j, j) == 0.0) {
            throw std::runtime_error("Linearly dependent column detected");
        }
        for (int k = 0; k < v.size(); ++k) {
            Q.col(j)[k] = v[k] / R(j, j);
        }
    }
}

std::vector<FaceData<Vector2>>
orthonormalize(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom, const std::vector<FaceData<Vector2>> &X) {
    MatrixX2d matrix(mesh.nFaces(), X.size());
    // TODO: Think about indexing
    for (std::size_t m = 0; m < X.size(); m++) {
        for (std::size_t n = 0; n < mesh.nFaces(); n++) {
            matrix(n, m) =
                X[m][mesh.face(n)] * std::sqrt(geom.faceAreas[mesh.face(n)]);
        }
    }
    MatrixX2d Q{};
    Eigen::MatrixXd R{};
    auto inner_product = [](const VectorX2d &a, const VectorX2d &b) -> double {
        double s = 0;
        for (long i = 0; i < a.size(); i++)
            s += dot(a[i], b[i]);
        return s;
    };
    modifiedGramSchmidt(matrix, Q, R, inner_product);
    std::vector<FaceData<Vector2>> h{};
    for (std::size_t m = 0; m < X.size(); m++) {
        FaceData<Vector2> &hm = h.emplace_back(mesh);
        for (std::size_t n = 0; n < mesh.nFaces(); n++) {
            Face f = mesh.face(n);
            hm[f] = Q(n, m) * (1. / std::sqrt(geom.faceAreas[f]));
        }
    }
    return h;
}

std::vector<FaceData<Vector2>>
orthonormal_hom_basis(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom, const std::vector<Singular_Circle> &homology_b) {
    PressureProjectionSolver pp_solver{};
    pp_solver.compute(geom);
    std::vector<FaceData<Vector2>> h(homology_b.size());
    for (std::size_t i = 0; i < homology_b.size(); i++) {
        auto &basis = homology_b[i];
        auto df = delta_form(mesh, basis);
        assert(df.raw().allFinite());
        EdgeData<double> pf = pp_solver.solve(mesh, df);
        assert(pf.raw().allFinite());
        h[i] = whitney_interpolation(mesh, geom, pf);
        assert(h[i].raw().allFinite());
    }
    return orthonormalize(mesh, geom, h);
}

std::vector<FaceData<Vector2>>
orthonormal_hom_basis(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom) {
    assert(&mesh == &geom.mesh);
    auto homotopy_b = greedy_homotopy_basis(mesh, geom, mesh.face(0));
    auto homology_b = singular_homology_basis(mesh, homotopy_b);
    return orthonormal_hom_basis(mesh, geom, homology_b);
}
} // namespace geometrycentral::surface
