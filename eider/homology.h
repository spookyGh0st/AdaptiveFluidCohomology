#pragma once
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "singular_homology.h"
#include <vector>

namespace geometrycentral::surface {

std::vector<Singular_Circle> singular_homology_basis(ManifoldSurfaceMesh &mesh, std::vector<std::vector<Halfedge>> hom_basis);

EdgeData<double> delta_form(ManifoldSurfaceMesh &mesh, const Singular_Circle &co_loop);
EdgeData<double> delta_form(ManifoldSurfaceMesh &mesh, const std::vector<Halfedge> &co_loop);

struct PressureProjectionSolver {
    void compute(IntrinsicGeometryInterface &geom);
    EdgeData<double> solve(ManifoldSurfaceMesh &mesh, const EdgeData<double> &co_loop) const;
    Eigen::SparseMatrix<double> A, AT;
    // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver {};
    // Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver{};
    // TODO: check wether i need QR to solve pseudoinverse - I hope not, but Im afraid I do!
    Eigen::ConjugateGradient<SparseMatrix<double>> solver;
};

struct AdaptivePressureProjectionSolver {
    void compute(IntrinsicGeometryInterface &geom);
    EdgeData<double> solveWithGuess(ManifoldSurfaceMesh &mesh, const EdgeData<double> &co_loop, VertexData<double>* guess);
    Eigen::SparseMatrix<double> A, AT;
    Eigen::ConjugateGradient<SparseMatrix<double>> solver {};
};

using Harmonic_basis = std::vector<FaceData<Vector2>>;

FaceData<Vector2> whitney_interpolation(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom, EdgeData<double> &h);

Harmonic_basis orthonormalize(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom, const std::vector<FaceData<Vector2>> &h);

Harmonic_basis orthonormal_hom_basis(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom, const std::vector<Singular_Circle> &homology_basis);

Harmonic_basis orthonormal_hom_basis(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom);
} // namespace geometrycentral::surface
