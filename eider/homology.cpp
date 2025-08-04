#include "homology.h"
#include <vector>
#include <algorithm>
#include <tuple>
#include <numeric>
#include <bitset>

#include "geometrycentral/surface/barycentric_vector.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"

namespace geometrycentral::surface
{
    // --- Disjoint Set Union (DSU) ---
    class DSU {
    public:
        std::vector<size_t> parent;
        std::vector<size_t> rank;

        explicit DSU(size_t n) {
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
    void computeMinimalSpanningTree(SurfaceMesh& mesh, EdgeData<EdgeType>& edgeData, const sp_cmp& fn) {


        std::vector<Edge> weightedEdges;
        weightedEdges.reserve(mesh.nEdges());

        for (Edge e : mesh.edges()) {
            if (edgeData[e] == EdgeType::bridge)
                weightedEdges.emplace_back(e);
        }


        std::ranges::sort(weightedEdges, fn);

        std::vector<Edge> mstEdges;

        DSU dsu(mesh.nVertices());
        size_t numMstEdges = 0;
        size_t targetNumEdges = mesh.nVertices() > 0 ? mesh.nVertices() - 1 : 0;

        for (Edge e : weightedEdges) {
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
    Halfedge computePrimalEdgesOfDualMaxST(SurfaceMesh& mesh, EdgeData<EdgeType>& edgeData, const sp_cmp& fn) {

        std::vector<Edge> dualEdges;
        dualEdges.reserve(mesh.nEdges());

        for (Edge primal_e : mesh.edges()) {
            if (primal_e.isBoundary()  || edgeData[primal_e] == EdgeType::minimal_st) {
                continue;
            }
            dualEdges.emplace_back(primal_e);
        }


        std::ranges::sort(dualEdges, fn);

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

        for (Edge e: primalEdgesInDualMaxST) {
            edgeData[e] = EdgeType::maximal_co_st;
        }

        // If manifold with boundary, add single boundary edge
        if (mesh.hasBoundary())
        {
            for (Edge e: mesh.edges()) {
                if (e.isBoundary()) {
                    edgeData[e] = EdgeType::maximal_co_st;
                    if (e.halfedge().isInterior()) return e.halfedge();
                    return e.halfedge().twin();
                }
            }
        }
        return Halfedge();
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
            heap.emplace_back(key, vertex);
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
                // TODO: If using dijkstra for optimal, use if (he.edge().isBoundary()) continue;
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

    // Compute for each edge e in G\T^*  the length of the shortest loop in T^* that contains e
    EdgeData<double> edge_coloop_lengths(SurfaceMesh& mesh, FaceData<double>& dist, const EdgeData<EdgeType>& data) {
        EdgeData<double> loop_length(mesh,NAN);
        for (Edge e: mesh.edges()) {
            if (data[e] == EdgeType::maximal_co_st) continue;
            double d = 0;
            for (Halfedge he: e.adjacentHalfedges()) if (he.isInterior()) d +=dist[he.face()];
            loop_length[e] = d;
        }
        return loop_length;
    }

    // computes the shortest loop of G \ T^* that contains e and x
    std::vector<Halfedge> minimal_coloop(FaceData<Halfedge> prev, Edge bridge, Halfedge x)
    {
        Halfedge start, end;
        std::vector<Halfedge> co_loop;
        if (bridge.isBoundary()) {
            start=x;
            end=bridge.halfedge().isInterior()? bridge.halfedge() : bridge.halfedge().twin();
        }else {
            start = bridge.halfedge();
            end = bridge.halfedge().twin();
        }
        co_loop.push_back(end);
        Face f = end.face();
        while (f != start.face()) {
            co_loop.push_back(prev[f]);
            assert(co_loop.back().isInterior());
            f = co_loop.back().face();
        }
        if (bridge.isBoundary()) co_loop.push_back(start.twin());
        return co_loop;
    }

    std::vector<Halfedge> homotopy_co_loop(FaceData<Halfedge>& prev, Face x, Edge bridge, Halfedge bound_dual_edge)
    {
        Halfedge he = bridge.halfedge();
        assert(he.isInterior());
        Face f = he.face();
        std::vector<Halfedge> co_loop;
        co_loop.push_back(he);
        while (f != x)
        {
            co_loop.push_back(prev[f]);
            f = prev[f].face();
        }

        std::vector<Halfedge> back_co_loop {};
        if (bridge.isBoundary()) { he = bound_dual_edge; back_co_loop.push_back(he.twin()); }
        else {he = bridge.halfedge().twin(); }
        f = he.face();

        while (f != x)
        {
            back_co_loop.push_back(prev[f].twin());
            f = prev[f].face();
        }
        for (int i = back_co_loop.size() - 1; i >= 0; --i) {
            co_loop.push_back(back_co_loop[i]);
        }
        return co_loop;
    }

    // TODO: move to down
    std::vector<Halfedge> reduce_co_loop(SurfaceMesh& mesh, const std::vector<Halfedge>& co_loop) {
        // this is in O(n), could be improved to be in size of co_loop
        std::vector<Halfedge> reduced {}; reduced.reserve(co_loop.size());
        HalfedgeData<bool> in_loop(mesh,false);
        for (Halfedge he : co_loop) { in_loop[he] = true; }
        for (Halfedge he : co_loop) {
            if (!he.edge().isBoundary() && in_loop[he] && in_loop[he.twin()]) continue;
            reduced.push_back(he);
        }
        return reduced;
    }

    double path_length(IntrinsicGeometryInterface& geom, const std::vector<Halfedge>& path) {
        double s = 0;
        for (Halfedge he: path) {
            s += geom.dualEdgeLengths[he.edge()];
        }
        return s;
    }

    Halfedge add_boundary_edge(Face x, EdgeData<EdgeType>& data)
    {
        Halfedge h;
        for (Halfedge he: x.adjacentHalfedges()) {
            if (he.edge().isBoundary()) { h = he; break; }
        }
        if (h == Halfedge()) return Halfedge();
        data[h.edge()] = EdgeType::maximal_co_st;
        return h;
    }

    std::vector<std::vector<Halfedge>> greedy_homotopy_basis(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom, Face x)
    {
        EdgeData<EdgeType> edge_data(mesh, EdgeType::bridge);
        geom.requireEdgeLengths(); geom.requireDualEdgeLengths();
        auto prev_dist = co_dijkstra(mesh, geom, edge_data, x);
        for (Face f: mesh.faces()) {
            if (f == x) continue;
            Edge e = prev_dist.first[f].edge(); assert(!e.isBoundary());
            edge_data[e] = EdgeType::maximal_co_st;
        }
        Halfedge x_he = add_boundary_edge(x,edge_data);
        auto coloop_lengths = edge_coloop_lengths(mesh,prev_dist.second,edge_data);
        computeMinimalSpanningTree(mesh, edge_data, [&l = coloop_lengths] (Edge a, const Edge& b) {return l[a] > l[b]; });

        auto dist_edges = distinctEdges(mesh, edge_data);
        std::vector<std::vector<Halfedge>> co_loops;
        for (Edge e : dist_edges)
        {
            co_loops.push_back(homotopy_co_loop(prev_dist.first,x,e,x_he));
        }
        return co_loops;
    }

    std::vector<std::vector<Halfedge>> minimal_greedy_homotopy_basis(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom)
    {
        // TODO add logic for non-boundary case
        std::vector<std::vector<Halfedge>> basis; double basis_length = std::numeric_limits<double>::max();
        for (Face y: mesh.faces()) {
            bool onBoundary = false;
            for (Edge e: y.adjacentEdges()) if (e.isBoundary()) onBoundary=true;
            if (!onBoundary) continue;

            double d = 0;
            std::vector<std::vector<Halfedge>> b = greedy_homotopy_basis(mesh, geom, y);
            for (const auto& path:b) d+= path_length(geom,path);
            if ( d < basis_length)
            {
                basis = b; basis_length = d;
            }
        }
        return basis;
    }

    std::vector<std::vector<Halfedge>> homotopy_basis(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom, Face x)
    {
        EdgeData<EdgeType> edge_data(mesh, EdgeType::bridge);
        geom.requireEdgeLengths(); geom.requireDualEdgeLengths();

        Halfedge b_halfedge = computePrimalEdgesOfDualMaxST(mesh, edge_data, [&l = geom.dualEdgeLengths] (Edge a, Edge b) { return l[a] < l[b]; });
        computeMinimalSpanningTree(mesh, edge_data, [&l = geom.edgeLengths] (Edge a, const Edge& b) {return l[a] > l[b]; });
        geom.unrequireEdgeLengths(); geom.unrequireDualEdgeLengths();

        auto dist_edges = distinctEdges(mesh, edge_data);
        auto prev = co_dijkstra(mesh, geom, edge_data, x);
        std::vector<std::vector<Halfedge>> co_loops;
        for (Edge e : dist_edges)
        {
            // co_loops.push_back(reduce_co_loop(mesh,minimal_co_loop(prev.first, x, e)));
            co_loops.push_back(homotopy_co_loop(prev.first, x, e,b_halfedge));
        }
        return co_loops;
    }

    inline static constexpr std::size_t n_bits = 64;
    EdgeData<std::bitset<n_bits>> homology_bitvectors(SurfaceMesh& mesh, const std::vector<std::vector<Halfedge>>& homotopy_basis) {
        EdgeData<std::bitset<n_bits>> bitsets(mesh, std::bitset<n_bits>(0));
        if (homotopy_basis.size() > n_bits) throw std::runtime_error("Genus exceeds bitvector size");
        for (std::size_t i = 0; i < homotopy_basis.size(); i++) {
            for (Halfedge he: homotopy_basis[i]) {
                // TODO: Should i set you or flip you - what about basis that go forth and back, i.e. contain an edge twice?
                bitsets[he.edge()][i].flip();
            }
        }
        return  bitsets;
    }

    bool cycle_contractable(EdgeData<std::bitset<n_bits>>& hom_bitvectors, const std::vector<Halfedge>& cycle) {
        std::bitset<n_bits> v {};
        for (Halfedge he : cycle) {
            v ^= hom_bitvectors[he.edge()];
        }
        return v.none();
    }

    void shortest_path_cotree() {
        // TODO
    }


    void co_cut_locus()
    {
        /*
        1. Compute the shortest path tree T from basepoint x (using Dijkstra).
        2. Construct the cut locus graph C:
        - For each edge e = (u,v) ∈ E not in T:
        a. Let P1 = shortest path x → u in T
        b. Let P2 = shortest path x → v in T
        c. Form cycle γ_e = P1 + (u,v) + reverse(P2)
        d. If γ_e is non-contractible (i.e., a valid loop):
        Add edge e to cut locus graph C
        3. Remove all degree-1 vertices from C iteratively.
        (This gives the reduced cut locus Φ)
        4. Return Φ
        */
    }


    EdgeData<double> delta_form(SurfaceMesh& mesh, const std::vector<Halfedge>& co_loop) {
        EdgeData<double> delta(mesh, 0);
        for (Halfedge he : co_loop) {
            delta[he.edge()] = he.orientation() ? 1 : -1;
        }
        return  delta;

    }


    void PressureProjectionSolver::compute(IntrinsicGeometryInterface& geom)
    {
        geom.requireDECOperators();
        A = geom.d0; AT = A.transpose();
        solver.compute(AT*A);
        // TODO: QR - solver.compute(A);
        geom.unrequireDECOperators();
    }

    EdgeData<double> PressureProjectionSolver::solve(SurfaceMesh& mesh, const EdgeData<double>& co_loop) const
    {
        Eigen::VectorXd x = co_loop.toVector();
        // TODO: QR - Eigen::VectorXd c = solver.solve(x);
        Eigen::VectorXd c = solver.solve(AT * x);
        return EdgeData<double>(mesh, x - A*c);
    }

    // Function to compute the Whitney interpolated vector field
    FaceData<Vector2> whitney_interpolation(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom, EdgeData<double>& h) {
        geom.requireHalfedgeVectorsInFace();
        geom.requireFaceAreas();
        FaceData<Vector2> vField(mesh, Vector2::zero());


        for (Face f : mesh.faces()) {
            Vector2 vT = Vector2::zero();
            double A = geom.faceAreas[f];
            assert(A > 0);
            auto grad = [&,A](Halfedge he) -> Vector2 {
                Vector2 ei = geom.halfedgeVectorsInFace[he.next()];
                return 1/(2*A) * ei.rotate90();
            };
            for (Halfedge he: f.adjacentHalfedges()){
                Edge e = he.edge();
                double wij = (he.orientation()) ? h[e] : -h[e];
                double A = geom.faceAreas[f];
                // vT += wij * grad(he) * A/3 - wij * grad(he.next()) * A/3;
                Vector2 ni = geom.halfedgeVectorsInFace[he.next()].rotate90(), nj = geom.halfedgeVectorsInFace[he.next().next()].rotate90();
                vT += wij /6 * (nj - ni);
            }
            assert(vT.isFinite());
            vField[f] = vT;
        }

        return vField;
    }


    using VectorX2d = Eigen::Matrix<Vector2,-1,1>;
    using MatrixX2d = Eigen::Matrix<Vector2,-1,-1>;

    using InnerProductFn = std::function<double(const VectorX2d&, const VectorX2d&)>;

    void modifiedGramSchmidt(const MatrixX2d& A, MatrixX2d& Q, Eigen::MatrixXd& R, const InnerProductFn& innerProduct) {
        const int m = A.rows();
        const int n = A.cols();

        Q = MatrixX2d::Zero(m, n);
        R = Eigen::MatrixXd::Zero(n, n);

        for (int j = 0; j < n; ++j) {
            VectorX2d v = A.col(j);

            for (int i = 0; i < j; ++i) {
                R(i, j) = innerProduct(Q.col(i), A.col(j));
                VectorX2d w = Q.col(i);
                for (int k = 0; k< v.size(); ++k) { v[k] -= w[k] * R(i,j); }
            }

            R(j, j) = std::sqrt(innerProduct(v, v));
            if (R(j, j) == 0.0) {
                throw std::runtime_error("Linearly dependent column detected");
            }
            for (int k = 0; k< v.size(); ++k) { Q.col(j)[k] = v[k] / R(j,j); }
        }
    }


    std::vector<FaceData<Vector2>> orthonormalize(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom,const std::vector<FaceData<Vector2>>& X) {
        MatrixX2d matrix(mesh.nFaces(), X.size());
        // TODO: Think about indexing
        for (std::size_t m = 0; m < X.size(); m++) {
            for (std::size_t n = 0; n < mesh.nFaces(); n++) {
                matrix(n,m) = X[m][mesh.face(n)] * std::sqrt(geom.faceAreas[mesh.face(n)]);
            }
        }
        MatrixX2d Q {}; Eigen::MatrixXd R{};
        auto inner_product = [](const VectorX2d& a, const VectorX2d& b) -> double {
            double s = 0; for (std::size_t i = 0; i < a.size(); i++) s += dot(a[i], b[i]); return s;
        };
        modifiedGramSchmidt(matrix,Q,R, inner_product);
        std::vector<FaceData<Vector2>> h{};
        for (std::size_t m = 0; m < X.size(); m++)
        {
            FaceData<Vector2>& hm = h.emplace_back(mesh);
            for (std::size_t n = 0; n < mesh.nFaces(); n++)
            {
                Face f = mesh.face(n);
                hm[f] = Q(n, m) * (1./ std::sqrt(geom.faceAreas[f]));
            }
        }
        return h;
    }

    std::vector<FaceData<Vector2>> orthonormal_hom_basis(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom, const std::vector<std::vector<Halfedge>>& homotopy_b)
    {
        PressureProjectionSolver pp_solver {};
        pp_solver.compute(geom);
        std::vector<FaceData<Vector2>> h(homotopy_b.size());
        for (std::size_t i = 0; i < homotopy_b.size(); i++){
            auto& basis = homotopy_b[i];
            auto df =delta_form(mesh, reduce_co_loop(mesh, basis));
            assert(df.raw().allFinite());
            EdgeData<double> pf = pp_solver.solve(mesh,df);
            assert(pf.raw().allFinite());
            h[i] = whitney_interpolation(mesh,geom,pf);
            assert(h[i].raw().allFinite());
        }
        return orthonormalize(mesh, geom, h);
    }

    std::vector<FaceData<Vector2>> orthonormal_hom_basis(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom)
    {
        assert(&mesh == &geom.mesh);
        auto homotopy_b= homotopy_basis(mesh, geom,mesh.face(0));
        return orthonormal_hom_basis(mesh,geom,homotopy_b);
    }
}

