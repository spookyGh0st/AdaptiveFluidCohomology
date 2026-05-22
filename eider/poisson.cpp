#include "poisson.h"

#include <Eigen/SparseCholesky> // For SimplicialLDLT
#include <Eigen/SparseCore>
#include <vector>

#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/surface_mesh.h>

namespace geometrycentral::surface {

std::string eigenInfoToString(const Eigen::ComputationInfo &info) {
    switch (info) {
    case Eigen::Success:
        return "Success";
    case Eigen::NumericalIssue:
        return "NumericalIssue: The provided data did not satisfy the prerequisites of the solver. For example, the matrix may be singular or not positive-definite for a Cholesky factorization.";
    case Eigen::NoConvergence:
        return "NoConvergence: The iterative solver did not converge to a solution within the given number of iterations.";
    case Eigen::InvalidInput:
        return "InvalidInput: The inputs (e.g., matrix, vector sizes) are invalid.";
    default:
        return "UnknownError: An unknown error occurred.";
    }
}

void StreamFunctionSolver::compute(SurfaceMesh &mesh, IntrinsicGeometryInterface &geom) {
    boundaryVertices.clear();
    interiorVertices.clear();

    if (mesh.hasBoundary()) {
        compute_dirichlet(mesh, geom);
    } else {
        compute_zero_mean(mesh, geom);
    }
}

void StreamFunctionSolver::solve(SurfaceMesh &mesh, IntrinsicGeometryInterface &geom, VertexData<double> &u, const VertexData<double> &f) const {
    if (mesh.hasBoundary()) {
        solve_dirichlet(u, f);
    } else {
        solve_zero_mean(mesh, u, f);
    }
}

void StreamFunctionSolver::compute_dirichlet(SurfaceMesh &mesh, IntrinsicGeometryInterface &geom) {
    assert(boundaryVertices.empty()); // must always be created new
    geom.requireCotanLaplacian();

    size_t nVertices = mesh.nVertices();
    if (nVertices == 0)
        return; // Nothing to do on an empty mesh.

    // === 1. Partition vertices and establish index mappings ===
    // We need separate numbering for interior and boundary vertices to build the sub-system.
    globalToInteriorIndex = VertexData(mesh, INVALID_IND);
    globalToBoundaryIndex = VertexData(mesh, INVALID_IND);

    size_t nInterior = 0;
    size_t nBoundary = 0;

    for (Vertex v : mesh.vertices()) {
        if (v.isBoundary()) {
            boundaryVertices.push_back(v);
            globalToBoundaryIndex[v] = nBoundary++;
        } else {
            interiorVertices.push_back(v);
            globalToInteriorIndex[v] = nInterior++;
        }
    }

    if (nInterior == 0)
        return; // Solution is just the boundary values already in f.
    if (nBoundary == 0) {
        throw std::runtime_error("solve_laplace error: Mesh has no boundary, Dirichlet problem is ill-posed.");
    }

    // === 2. Build the full Cotangent Laplacian Matrix (N x N) ===
    // This operator discretizes the Laplace-Beltrami operator (up to sign/normalization conventions).
    Eigen::SparseMatrix<double> &L_full = geom.cotanLaplacian;

    // === 3. Extract Submatrices L_II and L_IB ===
    // Required because the linear system separates known boundary values from unknown interior values.
    // We build using triplets because direct sparse submatrix slicing is not supported like dense.
    std::vector<Eigen::Triplet<double>> triplets_II, triplets_IB;
    triplets_II.reserve(nInterior * 7); // Pre-allocate assuming ~7 non-zeros per row (typical mesh sparsity)
    triplets_IB.reserve(nInterior * 3); // Guess fewer connections to boundary

    for (int k = 0; k < L_full.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(L_full, k); it; ++it) {
            Vertex rowV = mesh.vertex(it.row());
            Vertex colV = mesh.vertex(it.col());

            int interiorRowIdx = globalToInteriorIndex[rowV];
            int interiorColIdx = globalToInteriorIndex[colV];
            int boundaryColIdx = globalToBoundaryIndex[colV];

            // Only consider rows corresponding to interior vertices (equations we need to solve)
            if (interiorRowIdx != INVALID_IND) {
                if (interiorColIdx != INVALID_IND) {
                    // Interaction between two interior vertices -> goes into L_II
                    triplets_II.emplace_back(interiorRowIdx, interiorColIdx, it.value());
                } else if (boundaryColIdx != INVALID_IND) {
                    // Interaction between an interior and a boundary vertex -> goes into L_IB
                    triplets_IB.emplace_back(interiorRowIdx, boundaryColIdx, it.value());
                }
            }
        }
    }

    Eigen::SparseMatrix<double> L_II(nInterior, nInterior);
    L_IB = Eigen::SparseMatrix<double>(nInterior, nBoundary);
    L_II.setFromTriplets(triplets_II.begin(), triplets_II.end());
    L_IB.setFromTriplets(triplets_IB.begin(), triplets_IB.end());

    // === 4. Build the boundary value vector x_B ===
    // These are the known values on the right side, scaled by L_IB.

    std::vector<Eigen::Triplet<double>> triplets_MI;
    geom.requireVertexGalerkinMassMatrix();
    for (int k = 0; k < geom.vertexGalerkinMassMatrix.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(geom.vertexGalerkinMassMatrix, k); it; ++it) {
            Vertex rowV = mesh.vertex(it.row());
            Vertex colV = mesh.vertex(it.col());

            int interiorRowIdx = globalToInteriorIndex[rowV];
            int interiorColIdx = globalToInteriorIndex[colV];

            // Only consider rows corresponding to interior vertices (equations we need to solve)
            if (interiorRowIdx != -1) {
                if (interiorColIdx != -1) {
                    triplets_MI.emplace_back(interiorRowIdx, interiorColIdx, it.value());
                }
            }
        }
    }
    M_II = Eigen::SparseMatrix<double>(nInterior, nInterior);
    M_II.setFromTriplets(triplets_MI.begin(), triplets_MI.end());

    solver.compute(L_II);

    geom.unrequireCotanLaplacian();
    geom.unrequireVertexGalerkinMassMatrix();
}

void StreamFunctionSolver::solve_dirichlet(VertexData<double> &u,
                                           const VertexData<double> &f) const {
    if (interiorVertices.empty())
        return;
    Eigen::VectorXd x_B(boundaryVertices.size()), B_I(interiorVertices.size());
    for (Vertex v : boundaryVertices) {
        x_B(globalToBoundaryIndex[v]) = u[v];
    }
    for (Vertex v : interiorVertices) {
        B_I[globalToInteriorIndex[v]] = f[v];
    }

    Eigen::VectorXd rhs = (M_II * B_I) - (L_IB * x_B);
    Eigen::VectorXd x_I = solver.solve(rhs);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Solving dirichlet system failed.");
    }

    for (Vertex v : interiorVertices) {
        u[v] = x_I(globalToInteriorIndex[v]);
    }
}

void StreamFunctionSolver::compute_zero_mean(SurfaceMesh &mesh, IntrinsicGeometryInterface &geom) {
    size_t n = mesh.nVertices();
    if (n == 0)
        return;

    geom.requireCotanLaplacian();
    geom.requireVertexGalerkinMassMatrix();

    L_full = geom.cotanLaplacian;
    M_full = geom.vertexGalerkinMassMatrix;

    // Build mass vector
    massVector = Eigen::VectorXd::Zero(n);
    for (int k = 0; k < M_full.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(M_full, k); it; ++it) {
            massVector[it.row()] += it.value();
        }
    }

    // Build augmented matrix A
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(L_full.nonZeros() + 2 * n);

    for (int k = 0; k < L_full.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(L_full, k); it; ++it) {
            triplets.emplace_back(it.row(), it.col(), it.value());
        }
    }

    for (int i = 0; i < n; ++i) {
        double m_i = massVector[i];
        triplets.emplace_back(i, n, m_i);
        triplets.emplace_back(n, i, m_i);
    }

    A_zero_mean = Eigen::SparseMatrix<double>(n + 1, n + 1);
    A_zero_mean.setFromTriplets(triplets.begin(), triplets.end());

    solver.compute(A_zero_mean);
    if (solver.info() == Eigen::Success) {
        use_fallback = false;
    } else {
        std::cout << "\n"
                  << eigenInfoToString(solver.info()) << "\nUsing Fallback solver" << std::endl;
        use_fallback = true;
        fallback_solver.compute(A_zero_mean);
    }

    geom.unrequireCotanLaplacian();
    geom.unrequireVertexGalerkinMassMatrix();
}

void StreamFunctionSolver::solve_zero_mean(SurfaceMesh &mesh, VertexData<double> &u, const VertexData<double> &f) const {
    size_t n = u.size();
    if (n == 0)
        return;

    Eigen::VectorXd g_vec(n);
    for (Vertex v : mesh.vertices()) {
        g_vec[v.getIndex()] = f[v];
    }

    Eigen::VectorXd rhs = M_full * g_vec;

    Eigen::VectorXd rhs_aug(n + 1);
    rhs_aug.head(n) = rhs;
    rhs_aug[n] = 0.0;

    Eigen::VectorXd solution;
    if (use_fallback) {
        solution = fallback_solver.solve(rhs_aug);
        if (fallback_solver.info() != Eigen::Success)
            throw std::runtime_error(eigenInfoToString(fallback_solver.info()));
    } else {
        solution = solver.solve(rhs_aug);
        if (solver.info() != Eigen::Success)
            throw std::runtime_error(eigenInfoToString(solver.info()));
    }

    for (Vertex v : mesh.vertices()) {
        u[v] = solution[v.getIndex()];
    }
}
} // namespace geometrycentral::surface
