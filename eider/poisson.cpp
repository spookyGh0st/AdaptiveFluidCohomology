#include "poisson.h"

#include <vector>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky> // For SimplicialLDLT

#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>


namespace geometrycentral::surface
{
    void solve_poisson_dirichlet_zero_mean(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                                           VertexData<double>& f, const VertexData<double>& g)
    {
        geom.requireCotanLaplacian();
        geom.requireVertexDualAreas(); // Needed for area-weighted mean

        size_t nVertices = mesh.nVertices();
        if (nVertices == 0) return;

        // === 1. Build Laplacian ===
        Eigen::SparseMatrix<double>& L = geom.cotanLaplacian;

        // === 2. Build RHS vector (g) ===
        Eigen::VectorXd rhs(nVertices);
        for (Vertex v : mesh.vertices())
        {
            rhs[v.getIndex()] = g[v];
        }

        // === 3. Build area-weight vector for zero mean constraint ===
        Eigen::RowVectorXd areaRow(nVertices);
        for (Vertex v : mesh.vertices())
        {
            areaRow[v.getIndex()] = geom.vertexDualAreas[v];
        }

        // === 4. Augment system (L_aug * f = rhs_aug) ===
        Eigen::SparseMatrix<double> L_aug(nVertices + 1, nVertices);
        std::vector<Eigen::Triplet<double>> triplets;

        // Copy L into top rows of L_aug
        for (int k = 0; k < L.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it)
            {
                triplets.emplace_back(it.row(), it.col(), it.value());
            }
        }

        // Add areaRow as the last row in L_aug
        for (int j = 0; j < nVertices; ++j)
        {
            if (areaRow(j) != 0.0)
            {
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
        for (Vertex v : mesh.vertices())
        {
            f[v] = f_vec[v.getIndex()];
        }

        geom.unrequireCotanLaplacian();
        geom.unrequireVertexDualAreas();
    }



    void StreamFunctionSolver::compute(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom)
    {
        geom.requireCotanLaplacian();

        size_t nVertices = mesh.nVertices();
        if (nVertices == 0) return; // Nothing to do on an empty mesh.

        // === 1. Partition vertices and establish index mappings ===
        // We need separate numbering for interior and boundary vertices to build the sub-system.
        globalToInteriorIndex = VertexData(mesh, INVALID_IND);
        globalToBoundaryIndex = VertexData(mesh, INVALID_IND);

        interiorVertices;
        boundaryVertices;
        size_t nInterior = 0;
        size_t nBoundary = 0;

        for (Vertex v : mesh.vertices())
        {
            if (v.isBoundary())
            {
                boundaryVertices.push_back(v);
                globalToBoundaryIndex[v] = nBoundary++;
            }
            else
            {
                interiorVertices.push_back(v);
                globalToInteriorIndex[v] = nInterior++;
            }
        }

        if (nInterior == 0) return; // Solution is just the boundary values already in f.
        if (nBoundary == 0)
        {
            throw std::runtime_error("solve_laplace error: Mesh has no boundary, Dirichlet problem is ill-posed.");
        }

        // === 2. Build the full Cotangent Laplacian Matrix (N x N) ===
        // This operator discretizes the Laplace-Beltrami operator (up to sign/normalization conventions).
        Eigen::SparseMatrix<double>& L_full = geom.cotanLaplacian;

        // === 3. Extract Submatrices L_II and L_IB ===
        // Required because the linear system separates known boundary values from unknown interior values.
        // We build using triplets because direct sparse submatrix slicing is not supported like dense.
        std::vector<Eigen::Triplet<double>> triplets_II, triplets_IB;
        triplets_II.reserve(nInterior * 7); // Pre-allocate assuming ~7 non-zeros per row (typical mesh sparsity)
        triplets_IB.reserve(nInterior * 3); // Guess fewer connections to boundary

        for (int k = 0; k < L_full.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(L_full, k); it; ++it)
            {
                Vertex rowV = mesh.vertex(it.row());
                Vertex colV = mesh.vertex(it.col());

                int interiorRowIdx = globalToInteriorIndex[rowV];
                int interiorColIdx = globalToInteriorIndex[colV];
                int boundaryColIdx = globalToBoundaryIndex[colV];

                // Only consider rows corresponding to interior vertices (equations we need to solve)
                if (interiorRowIdx != INVALID_IND)
                {
                    if (interiorColIdx != INVALID_IND)
                    {
                        // Interaction between two interior vertices -> goes into L_II
                        triplets_II.emplace_back(interiorRowIdx, interiorColIdx, it.value());
                    }
                    else if (boundaryColIdx != INVALID_IND)
                    {
                        // Interaction between an interior and a boundary vertex -> goes into L_IB
                        triplets_IB.emplace_back(interiorRowIdx, boundaryColIdx, it.value());
                    }
                }
            }
        }

        Eigen::SparseMatrix<double> L_II(nInterior, nInterior);
        L_IB = Eigen::SparseMatrix<double> (nInterior, nBoundary);
        L_II.setFromTriplets(triplets_II.begin(), triplets_II.end());
        L_IB.setFromTriplets(triplets_IB.begin(), triplets_IB.end());

        // === 4. Build the boundary value vector x_B ===
        // These are the known values on the right side, scaled by L_IB.

        std::vector<Eigen::Triplet<double>> triplets_MI;
        geom.requireVertexGalerkinMassMatrix();
        for (int k = 0; k < geom.vertexGalerkinMassMatrix.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(geom.vertexGalerkinMassMatrix, k); it; ++it)
            {
                Vertex rowV = mesh.vertex(it.row());
                Vertex colV = mesh.vertex(it.col());

                int interiorRowIdx = globalToInteriorIndex[rowV];
                int interiorColIdx = globalToInteriorIndex[colV];

                // Only consider rows corresponding to interior vertices (equations we need to solve)
                if (interiorRowIdx != -1)
                {
                    if (interiorColIdx != -1)
                    {
                        triplets_MI.emplace_back(interiorRowIdx, interiorColIdx, it.value());
                    }
                }
            }
        }
        M_II = Eigen::SparseMatrix<double> (nInterior, nInterior);
        M_II.setFromTriplets(triplets_MI.begin(), triplets_MI.end());


        solver.compute(L_II);

        geom.unrequireCotanLaplacian();
        geom.unrequireVertexGalerkinMassMatrix();
    }

    void StreamFunctionSolver::solve_dirichlet(VertexData<double>& f,
        const VertexData<double>& g) const
    {
        Eigen::VectorXd x_B(boundaryVertices.size()), B_I(interiorVertices.size());
        for (Vertex v : boundaryVertices) { x_B(globalToBoundaryIndex[v]) = f[v]; }
        for (Vertex v : interiorVertices) { B_I[globalToInteriorIndex[v]] = g[v]; }

        Eigen::VectorXd rhs = (M_II * B_I) - (L_IB * x_B);
        Eigen::VectorXd x_I = solver.solve(rhs);

        for (Vertex v : interiorVertices) {
            f[v] = x_I(globalToInteriorIndex[v]);
        }
    }

    void StreamFunctionSolver::solve(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom, VertexData<double>& f,
        const VertexData<double>& g) const
    {
        if (mesh.hasBoundary())
        {
            solve_dirichlet(f,g);
        } else {
            solve_poisson_dirichlet_zero_mean(mesh, geom, f, g);
            // TODO: make more efficient
        }
    }
}
