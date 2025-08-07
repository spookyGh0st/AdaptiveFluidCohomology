#pragma once
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "singular_homology.h"
#include <vector>

namespace geometrycentral::surface {
enum EdgeType { minimal_st,
                maximal_co_st,
                bridge };

// compare functions for the spanning tree and co spanning tree
using sp_cmp = std::function<bool(Edge, Edge)>;

// TODO: extend to relative homology
// --- Primal Graph MST Computation ---
void computeMinimalSpanningTree(ManifoldSurfaceMesh &mesh, EdgeData<EdgeType> &edgeData, const sp_cmp &fn);

// --- Dual Graph Maximal Spanning Tree ---
Halfedge computePrimalEdgesOfDualMaxST(ManifoldSurfaceMesh &mesh, EdgeData<EdgeType> &edgeData, const sp_cmp &fn);

// --- Extract Edges not in MST or MaxST ---
std::vector<Edge> distinctEdges(ManifoldSurfaceMesh &mesh, EdgeData<EdgeType> &edgeData);

// --- Dual Dijkstra ---
std::pair<FaceData<Halfedge>, FaceData<double>> co_dijkstra(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom, EdgeData<EdgeType> &edgeData, Face orig_face);

// --- Minimal Co-loop ---
std::vector<Halfedge> homotopy_co_loop(FaceData<Halfedge> &prev, Face x, Edge bridge, Halfedge bound_dual_edge);

std::vector<Halfedge> reduce_co_loop(ManifoldSurfaceMesh &mesh, const std::vector<Halfedge> &co_loop);

// --- Homotopy Basis ---
std::vector<std::vector<Halfedge>> homotopy_basis(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom, Face x);

std::vector<std::vector<Halfedge>> greedy_homotopy_basis(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom, Face x);
std::vector<std::vector<Halfedge>> minimal_greedy_homotopy_basis(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom);

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
