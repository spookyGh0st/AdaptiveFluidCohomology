#pragma once

#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"

namespace geometrycentral::surface {
enum EdgeType { minimal_st,
    maximal_co_st,
    bridge };

// compare functions for the spanning tree and co spanning tree
using sp_cmp = std::function<bool(Edge, Edge)>;

// if m has boundary, return boundary face, else, some random face
Face arbitrary_base_face(ManifoldSurfaceMesh& mesh);


// TODO: extend to relative homology
// --- Primal Graph MST Computation ---
void computeMinimalSpanningTree(ManifoldSurfaceMesh &mesh, EdgeData<EdgeType> &edgeData, const sp_cmp &fn);

// --- Dual Graph Maximal Spanning Tree ---
Halfedge computePrimalEdgesOfDualMaxST(ManifoldSurfaceMesh &mesh, EdgeData<EdgeType> &edgeData, const sp_cmp &fn);

// --- Extract Edges not in MST or MaxST ---
std::vector<Edge> distinctEdges(ManifoldSurfaceMesh &mesh, EdgeData<EdgeType> &edgeData);

// --- Dual Dijkstra ---
std::pair<FaceData<Halfedge>, FaceData<double>> co_dijkstra(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom, EdgeData<EdgeType> &edgeData, Face orig_face, bool skip_co);

using Homotopy_cycle = std::vector<Halfedge>;
// --- Minimal Co-loop ---
Homotopy_cycle homotopy_co_loop(FaceData<Halfedge> &prev, Face x, Edge bridge, Halfedge bound_dual_edge);

Homotopy_cycle reduce_co_loop(ManifoldSurfaceMesh &mesh, const std::vector<Halfedge> &co_loop);

// --- Homotopy Basis ---
std::vector<Homotopy_cycle> homotopy_basis(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom, Face x);

std::vector<Homotopy_cycle> greedy_homotopy_basis(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom, Face x);
std::vector<Homotopy_cycle> minimal_greedy_homotopy_basis(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom);
}
