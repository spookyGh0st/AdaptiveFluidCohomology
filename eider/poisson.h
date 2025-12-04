#pragma once

#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <Eigen/IterativeLinearSolvers>

namespace geometrycentral::surface {
void solve_poisson_dirichlet_zero_mean(SurfaceMesh &mesh, IntrinsicGeometryInterface &geom, VertexData<double> &f, const VertexData<double> &g);

/// Solves -Δu = f in Ω, u = 0 in ∂Ω
struct StreamFunctionSolver {
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double>> fallback_solver;

    bool use_fallback = false;

    // For bounded mesh
    Eigen::SparseMatrix<double> M_II, L_IB;
    VertexData<std::size_t> globalToInteriorIndex, globalToBoundaryIndex; // global to interior/boundary indices
    std::vector<Vertex> boundaryVertices{}, interiorVertices{};

    // For zero-mean (closed mesh) case
    Eigen::SparseMatrix<double> A_zero_mean;
    Eigen::SparseMatrix<double> L_full, M_full;
    Eigen::VectorXd massVector;

    void compute_dirichlet(SurfaceMesh &mesh, IntrinsicGeometryInterface &geom);
    void compute_zero_mean(SurfaceMesh &mesh, IntrinsicGeometryInterface &geom);
    void compute(SurfaceMesh &mesh, IntrinsicGeometryInterface &geom);
    StreamFunctionSolver() = default;
    StreamFunctionSolver(SurfaceMesh &mesh, IntrinsicGeometryInterface &geom) { compute(mesh, geom); }

    void solve_dirichlet(VertexData<double> &u, const VertexData<double> &f) const;
    void solve_zero_mean(SurfaceMesh &mesh, VertexData<double> &u, const VertexData<double> &f) const;
    /// Solves -Δu = f in Ω, u = 0 in ∂Ω
    void solve(SurfaceMesh &mesh, IntrinsicGeometryInterface &geom, VertexData<double> &u, const VertexData<double> &f) const;
};
} // namespace geometrycentral::surface
