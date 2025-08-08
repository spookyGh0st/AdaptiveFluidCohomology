#pragma once
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "homotopy.h"
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
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver{};
    // TODO: check wether i need QR to solve pseudoinverse - I hope not, but Im afraid I do!
    // Eigen::SparseQR<SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
};

FaceData<Vector2> whitney_interpolation(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom, EdgeData<double> &h);

std::vector<FaceData<Vector2>> orthonormalize(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom, const std::vector<FaceData<Vector2>> &h);

std::vector<FaceData<Vector2>> orthonormal_hom_basis(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom, const std::vector<Singular_Circle> &homology_basis);

std::vector<FaceData<Vector2>> orthonormal_hom_basis(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom);
} // namespace geometrycentral::surface
