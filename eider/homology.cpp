#include "homology.h"
#include <vector>
#include <algorithm> // For std::sort
#include <tuple>     // For std::tuple
#include <numeric>   // For std::iota (for DSU initialization)

#include "geometrycentral/surface/barycentric_vector.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral::surface
{
    // --- Disjoint Set Union (DSU) ---
    class DSU {
    public:
        std::vector<size_t> parent;
        std::vector<size_t> rank;

        DSU(size_t n) {
            parent.resize(n);
            std::iota(parent.begin(), parent.end(), 0);
            rank.assign(n, 0);
        }

        size_t find(size_t i) {
            if (parent[i] == i)
                return i;
            return parent[i] = find(parent[i]);
        }

        bool unite(size_t i, size_t j) {
            size_t root_i = find(i);
            size_t root_j = find(j);

            if (root_i != root_j) {
                if (rank[root_i] < rank[root_j])
                    std::swap(root_i, root_j);

                parent[root_j] = root_i;

                if (rank[root_i] == rank[root_j])
                    rank[root_i]++;

                return true;
            }
            return false;
        }
    };

    // --- Primal Graph MST Computation ---
    // Computes the Minimal Spanning Tree of the graph of edges, using edge lengths as weights.
    std::vector<Edge> computeMinimalSpanningTree(SurfaceMesh& mesh, IntrinsicGeometryInterface& geometry) {


        geometry.requireEdgeLengths();
        std::vector<std::tuple<double, Edge>> weightedEdges;
        weightedEdges.reserve(mesh.nEdges());

        for (Edge e : mesh.edges()) {
            weightedEdges.emplace_back(geometry.edgeLengths[e], e);
        }

        std::sort(weightedEdges.begin(), weightedEdges.end()); // Ascending for MST

        std::vector<Edge> mstEdges;
        if (mesh.nVertices() == 0) {
            return mstEdges;
        }

        DSU dsu(mesh.nVertices());
        size_t numMstEdges = 0;
        size_t targetNumEdges = mesh.nVertices() > 0 ? mesh.nVertices() - 1 : 0;

        for (const auto& item : weightedEdges) {
            Edge e = std::get<1>(item);

            Halfedge he = e.halfedge();
            Vertex v1 = he.vertex();
            Vertex v2 = he.twin().vertex();

            size_t idx1 = v1.getIndex();
            size_t idx2 = v2.getIndex();

            if (dsu.unite(idx1, idx2)) {
                mstEdges.push_back(e);
                numMstEdges++;
                if (numMstEdges == targetNumEdges && mesh.nVertices() > 0) {
                    break;
                }
            }
        }

        return mstEdges;
    }


    // --- Dual Graph Maximal Spanning Tree - Primal Edge Return ---
    // Computes the Maximal Spanning Tree (MaxST) of the dual graph.
    // Returns the set of primal edges that are "crossed" by the edges of this dual MaxST.
    // The weight of a dual edge is the length of the primal edge it crosses.
    std::vector<Edge> computePrimalEdgesOfDualMaxST(SurfaceMesh& mesh, IntrinsicGeometryInterface& geometry) {

        geometry.requireEdgeLengths();

        std::vector<std::tuple<double, Edge, size_t, size_t>> potentialDualEdges;
        potentialDualEdges.reserve(mesh.nEdges());
        geometry.requireDualEdgeLengths();


        for (Edge primal_e : mesh.edges()) {
            if (!primal_e.isBoundary()) {
                Halfedge h = primal_e.halfedge();
                Face f1 = h.face();
                Face f2 = h.twin().face();
                size_t f1_idx = f1.getIndex();
                size_t f2_idx = f2.getIndex();
                double length = geometry.dualEdgeLengths[primal_e];
                potentialDualEdges.emplace_back(length, primal_e, f1_idx, f2_idx);
            }
        }
        geometry.unrequireDualEdgeLengths();

        std::ranges::reverse(potentialDualEdges);
        std::ranges::sort(potentialDualEdges,
                          [](const auto& a, const auto& b) {
                              return std::get<0>(a) > std::get<0>(b); // Sort descending by length
                          });

        std::vector<Edge> primalEdgesInDualMaxST;

        if (mesh.nFaces() == 0) {
            return primalEdgesInDualMaxST;
        }

        DSU dsu(mesh.nFaces());
        size_t numDualMaxSTEdges = 0;
        size_t targetNumEdges = (mesh.nFaces() > 0) ? (mesh.nFaces() - 1) : 0;

        for (const auto& item : potentialDualEdges) {
            Edge primal_e = std::get<1>(item);
            size_t f1_idx = std::get<2>(item);
            size_t f2_idx = std::get<3>(item);

            if (dsu.unite(f1_idx, f2_idx)) {
                primalEdgesInDualMaxST.push_back(primal_e);
                numDualMaxSTEdges++;
                // This condition is important: an ST (or spanning forest) will have N-C edges
                // where N is number of nodes (faces) and C is number of connected components.
                // For a single component, target is N-1.
                if (numDualMaxSTEdges == targetNumEdges && mesh.nFaces() > 0) {
                    break;
                }
            }
        }

        return primalEdgesInDualMaxST;
    }
}

