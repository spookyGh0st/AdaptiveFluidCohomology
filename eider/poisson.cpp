#include "poisson.h"

#include <vector>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky> // For SimplicialLDLT

#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>


using namespace geometrycentral;
using namespace geometrycentral::surface;

void solve_poisson_dirichlet_zero_mean(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                              VertexData<double>& f, const VertexData<double>& g) {

    geom.requireCotanLaplacian();
    geom.requireVertexDualAreas(); // Needed for area-weighted mean

    size_t nVertices = mesh.nVertices();
    if (nVertices == 0) return;

    // === 1. Build Laplacian ===
    Eigen::SparseMatrix<double>& L = geom.cotanLaplacian;

    // === 2. Build RHS vector (g) ===
    Eigen::VectorXd rhs(nVertices);
    for (Vertex v : mesh.vertices()) {
        rhs[v.getIndex()] = g[v];
    }

    // === 3. Build area-weight vector for zero mean constraint ===
    Eigen::RowVectorXd areaRow(nVertices);
    for (Vertex v : mesh.vertices()) {
        areaRow[v.getIndex()] = geom.vertexDualAreas[v];
    }

    // === 4. Augment system (L_aug * f = rhs_aug) ===
    Eigen::SparseMatrix<double> L_aug(nVertices + 1, nVertices);
    std::vector<Eigen::Triplet<double>> triplets;

    // Copy L into top rows of L_aug
    for (int k = 0; k < L.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it) {
            triplets.emplace_back(it.row(), it.col(), it.value());
        }
    }

    // Add areaRow as the last row in L_aug
    for (int j = 0; j < nVertices; ++j) {
        if (areaRow(j) != 0.0) {
            triplets.emplace_back(nVertices, j, areaRow(j));
        }
    }

    L_aug.setFromTriplets(triplets.begin(), triplets.end());

    // === 5. Augmented RHS ===
    Eigen::VectorXd rhs_aug(nVertices + 1);
    rhs_aug.head(nVertices) = rhs;
    rhs_aug(nVertices) = 0.0; // Enforce ∑ A_i * f_i = 0

    // === 6. Solve the least-squares problem ===
    Eigen::SparseMatrix<double> L_augT = L_aug.transpose();
    Eigen::SparseMatrix<double> normalMatrix = L_augT * L_aug;
    Eigen::VectorXd normalRHS = L_augT * rhs_aug;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(normalMatrix);
    Eigen::VectorXd f_vec = solver.solve(normalRHS);

    // === 7. Copy result into f ===
    for (Vertex v : mesh.vertices()) {
        f[v] = f_vec[v.getIndex()];
    }

    geom.unrequireCotanLaplacian();
    geom.unrequireVertexDualAreas();
}


void surface::solve_stream_function(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom, VertexData<double>& f,
    const VertexData<double>& g)
{
    if (mesh.hasBoundary()) {
        throw std::runtime_error("Not implemented yet");
    }else {
        solve_poisson_dirichlet_zero_mean(mesh,geom,f,g);
    }
}
