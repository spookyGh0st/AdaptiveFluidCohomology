#include "homology.h"
#include "homotopy.h"
#include <algorithm>
#include <bitset>
#include <numeric>
#include <tuple>
#include <vector>
#include <Eigen/SVD>
#include <Eigen/src/SVD/BDCSVD.h>


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

void checkMatrixCriteria(const Eigen::MatrixXd &d) {
    int m = d.rows();
    int n = d.cols();

    // Compute rank with SVD
    Eigen::BDCSVD<Eigen::MatrixXd> svd = Eigen::BDCSVD<Eigen::MatrixXd>(d, Eigen::ComputeThinU | Eigen::ComputeThinV);
    double tol = 1e-10;  // tolerance threshold
    int rank = (svd.singularValues().array() > tol).count();

    std::cout << "Matrix is " << m << "x" << n << "\n";
    std::cout << "Rank = " << rank << "\n";

    if (m >= n && rank == n) {
        std::cout << "-> Full column rank: use (d^T d)^(-1) d^T\n";
    } else if (m <= n && rank == m) {
        std::cout << "-> Full row rank: use d^T (d d^T)^(-1)\n";
    } else {
        std::cout << "-> Rank deficient: must use SVD or iterative solver for pseudoinverse\n";
    }
    std::cout << std::endl;
}

void PressureProjectionSolver::compute(IntrinsicGeometryInterface &geom) {
    geom.required0();
    checkMatrixCriteria(geom.d0);
    A = geom.d0;
    AT = A.transpose();
    solver.compute(AT * A);
    geom.unrequired0();
}

EdgeData<double>
PressureProjectionSolver::solve(ManifoldSurfaceMesh &mesh,
                                const EdgeData<double> &co_loop) const {
    Eigen::VectorXd x = co_loop.toVector();
    Eigen::VectorXd c = solver.solve(AT * x);
    return EdgeData<double>(mesh, x - A * c);
}

void AdaptivePressureProjectionSolver::compute(IntrinsicGeometryInterface &geom) {
    geom.required0();
    A = geom.d0;
    AT = A.transpose();
    solver.compute(AT * A);
    geom.unrequired0();
}

EdgeData<double> AdaptivePressureProjectionSolver::solveWithGuess(ManifoldSurfaceMesh &mesh, const EdgeData<double> &co_loop, VertexData<double> *guess) {
    if (guess == nullptr)
        guess = new VertexData<double>(mesh, 0);
    Eigen::VectorXd x = co_loop.toVector();
    assert(x.size() == mesh.nEdges());
    Eigen::VectorXd rhs = AT *x;
    Eigen::VectorXd guessV = guess->toVector();
    assert(guess->size() == rhs.size());
    Eigen::VectorXd c = solver.solveWithGuess(AT * x, guess->toVector());
    auto info = solver.info();
    std::cout << solver.iterations() << std::endl;
    if (info != Eigen::Success) {
        std::string msg = "Eigen iterative solver failed.\n";
        msg += "  Solver type: LeastSquaresConjugateGradient (or your solver)\n";
        msg += "  Iterations: " + std::to_string(solver.iterations()) + "\n";
        msg += "  Max iterations: " + std::to_string(solver.maxIterations()) + "\n";
        msg += "  Final error (rel. residual): " + std::to_string(solver.error()) + "\n";
        msg += "  Issue: ";
        if (info == Eigen::NoConvergence) {
            msg += "No convergence (did not reach tolerance).\n";
            msg += "). Retrying with regularization...\n";
            // --- Fallback 1: add Tikhonov regularization ---
            const double lambda = 1e-6; // tune as needed
            SparseMatrix<double> Areg = AT * A + lambda * SparseMatrix<double>(A.cols(), A.cols());
            Eigen::BiCGSTAB<SparseMatrix<double>> fallbackSolver;
            fallbackSolver.compute(Areg);
            c = fallbackSolver.solve(AT * x);

            if (fallbackSolver.info() != Eigen::Success) {
                throw std::runtime_error(msg);
            }
        } else if (info == Eigen::NumericalIssue) {
            msg += "Numerical issue encountered (instability or breakdown).\n";
            throw std::runtime_error(msg);
        } else if (info == Eigen::InvalidInput) {
            msg += "Invalid input or improper solver use.\n";
            throw std::runtime_error(msg);
        } else {
            msg += "Unknown error code.\n";
            throw std::runtime_error(msg);
        }
    }

    *guess = VertexData<double>(mesh,c);
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
        vField[f] = vT / A;
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

// TODO: Use Hausedolder Decomp, it is faster and has better numerical properties.
std::vector<FaceData<Vector2>>
orthonormalize(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom, const std::vector<FaceData<Vector2>> &X) {
    geom.requireFaceAreas();
    MatrixX2d matrix(mesh.nFaces(), X.size());
    auto f_idx = mesh.getFaceIndices();
    for (std::size_t m = 0; m < X.size(); m++) {
        for (Face f: mesh.faces()){
            matrix(f_idx[f], m) = X[m][f];
        }
    }
    MatrixX2d Q{};
    Eigen::MatrixXd R{};
    Eigen::VectorXd face_vec = geom.faceAreas.toVector(f_idx);
    auto inner_product = [&face_vec](const VectorX2d &a, const VectorX2d &b) -> double {
        double s = 0;
        for (long i = 0; i < a.size(); i++)
            s += dot(a[i], b[i]) * face_vec[i];
        return s;
    };
    modifiedGramSchmidt(matrix, Q, R, inner_product);
    std::vector<FaceData<Vector2>> h{};
    for (std::size_t m = 0; m < X.size(); m++) {
        FaceData<Vector2> &hm = h.emplace_back(mesh);
        for (Face f: mesh.faces()){
            hm[f] = Q(f_idx[f], m);
        }
    }
    geom.unrequireFaceAreas();
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
    std::vector<Homotopy_cycle> homotopy_b = greedy_homotopy_basis(mesh,geom,arbitrary_base_face(mesh));
    auto homology_b = singular_homology_basis(mesh, homotopy_b);
    return orthonormal_hom_basis(mesh, geom, homology_b);
}
} // namespace geometrycentral::surface
