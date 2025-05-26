#include "homology.h"
#include <vector>
#include <algorithm> // For std::sort
#include <tuple>     // For std::tuple
#include <numeric>   // For std::iota (for DSU initialization)

#include "geometrycentral/surface/barycentric_vector.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"

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

        /**
         * replaces the set containing i and the set containing y with their union.
         * @param i first element
         * @param j second element
         * @return false, if already in same set, true if union was created
         */
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
    void computeMinimalSpanningTree(SurfaceMesh& mesh, IntrinsicGeometryInterface& geometry, EdgeData<EdgeType>& edgeData) {


        geometry.requireEdgeLengths();
        std::vector<std::tuple<double, Edge>> weightedEdges;
        weightedEdges.reserve(mesh.nEdges());

        for (Edge e : mesh.edges()) {
            weightedEdges.emplace_back(geometry.edgeLengths[e], e);
        }

        std::sort(weightedEdges.begin(), weightedEdges.end()); // Ascending for MST

        std::vector<Edge> mstEdges;
        if (mesh.nVertices() == 0) {
            return ;
        }

        DSU dsu(mesh.nVertices());
        size_t numMstEdges = 0;
        size_t targetNumEdges = mesh.nVertices() > 0 ? mesh.nVertices() - 1 : 0;

        for (const auto& item : weightedEdges) {
            Edge e = std::get<1>(item);
            Halfedge he = e.halfedge();

            size_t idx1 = he.tailVertex().getIndex();
            size_t idx2 = he.tipVertex().getIndex();

            if (dsu.unite(idx1, idx2)) {
                mstEdges.push_back(e);
                numMstEdges++;
                if (numMstEdges == targetNumEdges && mesh.nVertices() > 0) {
                    break;
                }
            }
        }

        for (Edge e: mstEdges)
        {
            edgeData[e] = EdgeType::minimal_st;
        }
    }


    // --- Dual Graph Maximal Spanning Tree - Primal Edge Return ---
    // Computes the Maximal Spanning Tree (MaxST) of the dual graph.
    // Returns the set of primal edges that are "crossed" by the edges of this dual MaxST.
    // The weight of a dual edge is the length of the primal edge it crosses.
    void computePrimalEdgesOfDualMaxST(SurfaceMesh& mesh, IntrinsicGeometryInterface& geometry, EdgeData<EdgeType>& edgeData) {

        geometry.requireEdgeLengths();

        std::vector<Edge> dualEdges;
        dualEdges.reserve(mesh.nEdges());
        geometry.requireDualEdgeLengths();

        for (Edge primal_e : mesh.edges()) {
            if (primal_e.isBoundary()  || edgeData[primal_e] == EdgeType::minimal_st) {
                continue;
            }
            dualEdges.emplace_back(primal_e);
        }

        // Sort them by dual edge lengths. First do so by edge lengths and reserve to ensure distinct sets tho.
        std::ranges::sort(dualEdges, [&geometry](const auto& a, const auto& b) { return geometry.edgeLengths[a] < geometry.edgeLengths[b]; });
        std::ranges::reverse(dualEdges);
        std::ranges::sort(dualEdges, [&geometry](const auto& a, const auto& b) { return geometry.dualEdgeLengths[a] > geometry.dualEdgeLengths[b]; });

        std::vector<Edge> primalEdgesInDualMaxST;
        primalEdgesInDualMaxST.reserve(mesh.nEdges());

        DSU dsu(mesh.nFaces());
        size_t numDualMaxSTEdges = 0;
        size_t targetNumEdges = (mesh.nFaces() > 0) ? (mesh.nFaces() - 1) : 0;

        for (const auto& primal_e : dualEdges) {
            Face f1 = primal_e.halfedge().face();
            Face f2 = primal_e.halfedge().twin().face();

            if (dsu.unite(f1.getIndex(), f2.getIndex())) { // not same set
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

        geometry.unrequireDualEdgeLengths();
        for (Edge e: primalEdgesInDualMaxST)
        {
            edgeData[e] = EdgeType::maximal_co_st;
        }
    }

    std::vector<Edge> distinctEdges(SurfaceMesh& mesh, EdgeData<EdgeType>& edgeData) {
        std::vector<Edge> neither;
        for (Edge e: mesh.edges()) {
            if (edgeData[e] == EdgeType::bridge) neither.push_back(e);
        }
        return neither;
    }

    class MinHeap
    {
        using T = Face;
        std::vector<std::pair<double, T>> heap;
        FaceData<int> pos{}; // pos[vertex] = index in heap, or -1 if not present

        void swapNodes(int i, int j)
        {
            std::swap(heap[i], heap[j]);
            pos[heap[i].second] = i;
            pos[heap[j].second] = j;
        }

        void heapifyUp(int idx)
        {
            while (idx > 0)
            {
                int parent = (idx - 1) / 2;
                if (heap[idx].first < heap[parent].first)
                {
                    swapNodes(idx, parent);
                    idx = parent;
                }
                else
                {
                    break;
                }
            }
        }

        void heapifyDown(int idx)
        {
            int n = static_cast<int>(heap.size());
            while (true)
            {
                int left = 2 * idx + 1;
                int right = 2 * idx + 2;
                int smallest = idx;

                if (left < n && heap[left].first < heap[smallest].first)
                    smallest = left;
                if (right < n && heap[right].first < heap[smallest].first)
                    smallest = right;

                if (smallest != idx)
                {
                    swapNodes(idx, smallest);
                    idx = smallest;
                }
                else
                {
                    break;
                }
            }
        }

    public:
        explicit MinHeap(SurfaceMesh& m) : pos(m,-1)
        {
            heap.reserve(m.nFaces());
        }

        bool empty() const { return heap.empty(); }

        void insert(const T& vertex, double key)
        {
            heap.push_back({key, vertex});
            int idx = static_cast<int>(heap.size()) - 1;
            pos[vertex] = idx;
            heapifyUp(idx);
        }

        std::pair<double, T> extractMin()
        {
            if (empty()) return {std::numeric_limits<double>::infinity(), T()};

            std::pair<double, T> root = heap[0];
            std::pair<double, T> last = heap.back();

            heap[0] = last;
            pos[last.second] = 0;
            heap.pop_back();
            pos[root.second] = -1;

            if (!empty())
                heapifyDown(0);

            return root;
        }

        void decreaseKey(const T& vertex, double newKey)
        {
            int idx = pos[vertex];
            if (idx == -1) return;
            if (heap[idx].first <= newKey) return;

            heap[idx].first = newKey;
            heapifyUp(idx);
        }

        double value(const T& vertex) const
        {
            int idx = pos[vertex];
            if (idx == -1) return std::numeric_limits<double>::max();
            return heap[idx].first;
        }
    };

    std::pair<FaceData<Halfedge>,FaceData<double>> co_dijkstra(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom, EdgeData<EdgeType>& edgeData, Face orig_face)
    {
        geom.requireDualEdgeLengths();
        FaceData<Halfedge> prev(mesh, Halfedge());
        FaceData<double> dist(mesh, std::numeric_limits<double>::max());
        // Halfedge such that he.face is prev and he.edge is edge adjacent to this and prev
        MinHeap Q(mesh);
        for (Face f : mesh.faces())
        {
            Q.insert(f, std::numeric_limits<double>::max());
        }
        Q.decreaseKey(orig_face, 0);
        dist[orig_face] = 0;

        while (!Q.empty())
        {
            auto u = Q.extractMin();
            for (Halfedge he : u.second.adjacentHalfedges())
            {
                if (he.edge().isBoundary() || edgeData[he.edge()] != EdgeType::maximal_co_st) continue;
                Face v = he.twin().face();
                double alt = u.first + geom.dualEdgeLengths[he.edge()];
                if (alt < dist[v])
                {
                    Q.decreaseKey(v, alt);
                    dist[v] = alt; prev[v] = he;
                }
            }
        }
        geom.unrequireDualEdgeLengths();
        return {prev, dist};
    }

    std::vector<Halfedge> minimal_co_loop(FaceData<Halfedge>& prev, Face x, Edge bridge)
    {
        Halfedge he = bridge.halfedge();
        Face f = he.face();
        std::vector<Halfedge> co_loop;
        co_loop.push_back(he);
        while (f != x)
        {
            co_loop.push_back(prev[f]);
            f = prev[f].face();
        }
        f = he.twin().face();
        while (f != x)
        {
            co_loop.push_back(prev[f].twin());
            f = prev[f].face();
        }
        return co_loop;
    }

    std::vector<std::vector<Halfedge>> homotopy_basis(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom, Face x)
    {
        EdgeData<EdgeType> edge_data(mesh, EdgeType::bridge);
        computeMinimalSpanningTree(mesh, geom, edge_data);
        computePrimalEdgesOfDualMaxST(mesh, geom, edge_data);
        auto dist_edges = distinctEdges(mesh, edge_data);
        auto prev = co_dijkstra(mesh, geom, edge_data, x);
        std::vector<std::vector<Halfedge>> co_loops;
        for (Edge e : dist_edges)
        {
            co_loops.push_back(minimal_co_loop(prev.first, x, e));
        }
        return co_loops;
    }
}

