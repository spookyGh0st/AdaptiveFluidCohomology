#pragma once
#include <vector>
#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral::surface
{
    enum EdgeType {
        minimal_st, maximal_co_st, bridge
    };
    // TODO: extend to relative homology
    // --- Primal Graph MST Computation ---
    void computeMinimalSpanningTree(SurfaceMesh& mesh, IntrinsicGeometryInterface& geometry, EdgeData<EdgeType>& edgeData);

    // --- Dual Graph Maximal Spanning Tree ---
    void computePrimalEdgesOfDualMaxST(SurfaceMesh& mesh, IntrinsicGeometryInterface& geometry, EdgeData<EdgeType>& edgeData);

    // --- Extract Edges not in MST or MaxST ---
    std::vector<Edge> distinctEdges(SurfaceMesh& mesh, EdgeData<EdgeType>& edgeData);

    // --- Dual Dijkstra ---
    std::pair<FaceData<Halfedge>,FaceData<double>> co_dijkstra(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom, EdgeData<EdgeType>& edgeData, Face orig_face);

    // --- Minimal Co-loop ---
    std::vector<Halfedge> minimal_co_loop(FaceData<Halfedge>& prev, Face x, Edge bridge);

    // --- Homotopy Basis ---
    std::vector<std::vector<Halfedge>> homotopy_basis(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom, Face x);

    EdgeData<double> delta_form(SurfaceMesh& mesh, const std::vector<Halfedge>& co_loop);

    EdgeData<double> pressure_project(SurfaceMesh& mesh, const EdgeData<double>& co_loop, IntrinsicGeometryInterface& geom);

    FaceData<Vector2> whitney_interpolation(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom, EdgeData<double>& h);

    std::vector<FaceData<Vector2>> orthonormalize(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom, const std::vector<FaceData<Vector2>>& h);

    std::vector<FaceData<Vector2>> orthonormal_hom_basis(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom);
}

