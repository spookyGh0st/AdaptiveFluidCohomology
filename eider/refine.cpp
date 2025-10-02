#include "refine.h"

#include <unordered_set>
#include "util.h"

namespace geometrycentral::surface {

CornerData<bool> mark_longest_edge(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom) {
    CornerData<bool> marked_corner(mesh,false);
    geom.requireEdgeLengths();
    for (Face f : mesh.faces()) {
        double maxEl = std::numeric_limits<double>::lowest();
        Corner c;
        for (Halfedge he : f.adjacentHalfedges()) {
            if (geom.edgeLengths[he.edge()] > maxEl) {
                maxEl = geom.edgeLengths[he.edge()];
                c = he.oppositeCorner();
            }
        }
        marked_corner[c] = true;
    }
    geom.unrequireEdgeLengths();
    return marked_corner;
}

CornerData<bool> mark_pattern(ManifoldSurfaceMesh& mesh) {
    CornerData<bool> marked_corner(mesh,false);
    std::vector<Halfedge> se;
    auto isMarked = [&](Face f) {for (Corner c: f.adjacentCorners()) { if (marked_corner[c]) return true; } return false; };
    for (Edge e: mesh.edges()) {
        if (!e.isBoundary()) {
            se.push_back(e.halfedge()); break;
        }
    }
    while (!se.empty()) {
        Halfedge e = se.back(); se.pop_back();
        std::vector<Halfedge> border;
        for (Halfedge he: e.edge().adjacentHalfedges()) {
            if (!he.isInterior()) continue;
            if (isMarked(he.face())) continue;
            marked_corner[he.next().next().corner()] = true;
            if (he.next().twin().isInterior())
                border.push_back(he.next().twin().prevOrbitFace());
            if (he.next().next().twin().isInterior())
                border.push_back(he.next().next().twin().next());
        }
        for (Halfedge he:border) {
            if (isMarked(he.face())) continue;
            se.push_back(he);
        }
    }
    return marked_corner;
}

CornerData<bool> mark_random(ManifoldSurfaceMesh& mesh) {
    CornerData<bool> marked_corner(mesh,false);
    std::mt19937 gen(22);
    std::uniform_int_distribution<int> dist(0, 3);
    for (Face f: mesh.faces()) {
        Halfedge he = f.halfedge();
        for (int i = 0; i < dist(gen); ++i) {
            he = he.next();
        }
        marked_corner[he.corner()] = true;
    }
    return marked_corner;
}

CornerData<bool> mark_faces(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, MARKING_STRATEGY marking) {
    switch (marking) {
    case MARKING_STRATEGY::LONGEST_EDGE:
        return mark_longest_edge(mesh,geom);
    case MARKING_STRATEGY::PATTERN:
        return mark_pattern(mesh);
    case MARKING_STRATEGY::RANDOM:
        return mark_random(mesh);
    default: ;
        throw std::runtime_error("Not Implemented");
    }
}

AdaptiveTriangulation::AdaptiveTriangulation(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, MARKING_STRATEGY marking)
    :tri(mesh,geom), idx(*tri.intrinsicMesh), marked_corner(mark_faces(*tri.intrinsicMesh, tri,marking)) {
}

Halfedge AdaptiveTriangulation::vertex_bisection(Halfedge he, AdaptiveTransfer* transfer) {
    assert (he.isInterior());
    assert(marked_corner[he.oppositeCorner()]);
    if (he.twin().isInterior()) { assert(marked_corner[he.twin().oppositeCorner()]); }

    Halfedge twin_he = he.twin();
    Vertex vi = he.tailVertex(), vj = he.tipVertex();

    he = tri.splitEdge(he, 0.5);

    // If transfer is supplied, update that
    if(transfer) { transfer->refineEdge(vi,vj,he.tailVertex()); }

    // ensure indices are kept in order
    assert(twin_he.tipVertex() == he.tailVertex());
    assert(he.tipVertex() == vj );
    if (he.isInterior()) {
        Face fl = he.prevOrbitFace().twin().face(), fr = he.face();
        if (idx[fl] > idx[fr]) std::swap(idx[fl],idx[fr]);
    }
    assert(twin_he.tailVertex() == vj );
    if (twin_he.isInterior()) {
        Face fl = twin_he.face(), fr = twin_he.next().twin().face();
         if (idx[fl] > idx[fr]) std::swap(idx[fl],idx[fr]);
    }

    // Update refinement edges
    for (Halfedge he_o : he.vertex().outgoingHalfedges()) {
        if (!he_o.isInterior()) continue;
        setRefinementEdge(he_o.next());
    }

    return he;
}

void AdaptiveTriangulation::refine(std::vector<Face> faces, AdaptiveTransfer* transfer) {
    FaceData<bool> marked_faces (mesh(),false);
    std::unordered_set<Edge> start_edges;

    start_edges.reserve(faces.size()*2);

    if(transfer) transfer->startRefine();

    while (!faces.empty()) {
        Face f = faces.back();
        faces.pop_back();
        marked_faces[f] = true;
        Halfedge he = getRefinementEdge(f);
        if (!he.edge().isBoundary() && !marked_faces[he.twin().face()]) {
            faces.push_back(he.twin().face());
        }
        if (he.edge().isBoundary() || getRefinementEdge(he.twin().face()) == he.twin()) {
            start_edges.insert(he.edge());
        }
    }

    while (!start_edges.empty()) {
        Edge e = *start_edges.begin();
        start_edges.erase(start_edges.begin());

        if(e.halfedge().isInterior()) marked_faces[e.halfedge().face()] = false;
        if(e.halfedge().twin().isInterior()) marked_faces[e.halfedge().twin().face()] = false;

        Halfedge split_he = e.halfedge().isInterior() ? e.halfedge() : e.halfedge().twin();
        Vertex new_v = vertex_bisection(split_he,transfer).vertex();


        // Check if we can now refine one of the neighbours, i.e.
        // if the primal edge of one of the neighbours is one
        // of the newly created primal edges
        for (Halfedge he_o : new_v.outgoingHalfedges()) {
            if (!he_o.isInterior())
                continue;
            Halfedge side_he = he_o.next();
            if (side_he.edge().isBoundary())
                continue;
            Face side_f = side_he.twin().face();
            if (!marked_faces[side_f])
                continue;
            if (getRefinementEdge(side_f).edge() == side_he.edge())
                start_edges.insert(side_he.edge());
        }
    }
    assert(!marked_faces.raw().any()); // If this is the case, then we didn't process all faces!

    if(transfer) transfer->endRefine();
    tri.refreshQuantities();
}

Halfedge AdaptiveTriangulation::coarse_halfedge(Vertex v) {
    assert((!v.isBoundary() && v.faceDegree() == 4) || (v.isBoundary() && v.faceDegree() == 2));
    size_t k = 0;
    std::size_t min_idx = std::numeric_limits<std::size_t>::max();
    Halfedge he;
    for (Halfedge ohe : v.outgoingHalfedges()) {
        if (!ohe.isInterior()) continue;
        std::size_t f_idx = idx[ohe.face()];
        if (f_idx < min_idx) {
            min_idx = f_idx;
            he = ohe;
        };
    }
    // TODO: This should not be necessary
    // if (v.isBoundary()) { for (Halfedge ohe: v.outgoingHalfedges()) if (!ohe.edge().isBoundary()) he = ohe; }
    assert(he != Halfedge());
    assert(!he.edge().isBoundary());
    return he.twin().next();
}

Halfedge AdaptiveTriangulation::vertex_biunion(Halfedge he, AdaptiveTransfer* transfer) {
    // TODO: Assert left face has smaller idx then right face
    assert(he.isInterior());
    std::size_t l_idx = idx[he.prevOrbitFace().twin().face()];
    std::size_t r_idx = he.twin().isInterior()? idx[he.twin().face()] : -1;

    // Assert all adjacent faces have the refinement edges opposite to this vertex
    for (Halfedge ohe: he.vertex().outgoingHalfedges()) {
        if (!ohe.isInterior()) continue;
        assert(getRefinementEdge(ohe.face()) == ohe.next());
    }

    Vertex vi = he.prevOrbitFace().twin().next().tipVertex();
    Vertex vj = he.tipVertex();
    Vertex vp = he.tailVertex();
    he = tri.collapseEdgeTriangular(he);
    assert(vi == he.tailVertex());
    assert(vj == he.tipVertex());

    // Update Inverse coarsening map
    if(transfer){ transfer->coarseEdge(vi,vj,vp); }

    // update refinement edges
    for (Halfedge ahe: he.edge().adjacentHalfedges()) { if (ahe.isInterior()) setRefinementEdge(ahe); }

    // update index;
    idx[he.face()] = l_idx;
    if (he.twin().isInterior()) {
        idx[he.twin().face()] = r_idx;
    }

    return he;
}

inline bool vertex_is_good(IntrinsicTriangulation &T, Vertex v) {
    if (v == Vertex())
        return false;
    if (T.vertexLocations[v].type == SurfacePointType::Vertex)
        return false;
    return (v.isBoundary() && v.faceDegree() == 2) || (!v.isBoundary() && v.faceDegree() == 4);
}

inline bool vertexMarked(Vertex v, const FaceData<bool>& marked_faces){
    for(Face f:v.adjacentFaces()) if (!marked_faces[f]) return false;
    return true;
}

void AdaptiveTriangulation::coarse(const std::vector<Face> &faces, AdaptiveTransfer* transfer) {
    if(transfer) transfer->startCoarse();

    // Mark all potential good vertices, i.e. these vertices with only marked faces around it
    FaceData<bool> marked_faces(mesh(),false);
    for (Face f: faces){ marked_faces[f] = true;}

    std::unordered_set<Vertex> good_vertices{};
    for (Vertex v : tri.intrinsicMesh->vertices()) {
        if (vertex_is_good(tri, v) && vertexMarked(v,marked_faces))
            good_vertices.insert(v);
    }

    while (!good_vertices.empty()) {
        // pop element
        Vertex v = *good_vertices.begin();
        good_vertices.erase(v);

        Halfedge he = coarse_halfedge(v);
        he = vertex_biunion(he,transfer);
        assert(he != Halfedge());

        std::array<Vertex, 2> adjacent_v;
        adjacent_v[0] = he.next().tipVertex();
        if (he.twin().isInterior()){
            adjacent_v[1] = he.twin().next().tipVertex();

        }

        for (Vertex nv : adjacent_v) {
            if (vertex_is_good(tri, nv) && vertexMarked(nv,marked_faces))
                good_vertices.insert(nv);
        }

        // unmark adjacent faces
        for(Halfedge he: he.edge().adjacentHalfedges())
            if(he.isInterior()) marked_faces[he.face()] = false;
    }
    if (transfer) transfer->endCoarse();
    tri.refreshQuantities();
}

double etaRSqr(Face T, IntrinsicGeometryInterface &geom, const VertexData<double> &f, const VertexData<double> &u) {
    double f_st = 0, h_t = diameter(geom, T);
    for (Vertex v : T.adjacentVertices())
        f_st += f[v];
    f_st /= 3;

    double lu = 0;
    for (Vertex v : T.adjacentVertices())
        lu += laplacian(geom, v, u);
    lu /= 3;

    double jump_sum = 0;
    for (Edge e : T.adjacentEdges()) {
        if (e.isBoundary())
            continue; // only sum over interior edges
        double h_e = diameter(geom, e), l_e = geom.edgeLengths[e], j = 0;
        for (Halfedge he : e.adjacentHalfedges()) // for both sides of edge
        {
            j += dot(grad(geom, he.face(), u), geom.halfedgeVectorsInFace[he].rotateCW(PI / 2));
        }
        jump_sum += h_e * h_e * j * j; // one_he_e form L_2 norm.
    }
    jump_sum *= 1. / 2.;
    return h_t * h_t * (std::pow(f_st + lu, 2) * geom.faceAreas[T]) + jump_sum;
}

void AdaptiveTriangulation::setRefinementEdge(Halfedge he) {
    for (Halfedge fhe: he.face().adjacentHalfedges())
        marked_corner[fhe.corner()] = false;
    marked_corner[he.oppositeCorner()] = true;
}

Halfedge AdaptiveTriangulation::getRefinementEdge(Face f) {
    Halfedge he;
    for (Halfedge fhe: f.adjacentHalfedges())
        if (marked_corner[fhe.oppositeCorner()]) he= fhe;
    return he;
}

FaceData<double> poisson_residual_error_sqr(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom, const VertexData<double> &u, const VertexData<double> &f) {
    geom.requireHalfedgeVectorsInFace();
    geom.requireEdgeCotanWeights();
    FaceData<double> eta(mesh);
    for (Face T : mesh.faces()) {
        eta[T] = etaRSqr(T, geom, f, u);
    }
    geom.unrequireHalfedgeVectorsInFace();
    geom.unrequireEdgeCotanWeights();
    return eta;
}

std::array<std::vector<Face>,2> select_doerfler(ManifoldSurfaceMesh &mesh, FaceData<double> residual, const DoeflerConf& conf) {
    if (conf.theta_refine == 1) {
        std::vector<Face> faces; faces.reserve(mesh.nFaces());
        for (Face f:mesh.faces()) faces.push_back(f); return {faces,{}};
    }
    if (conf.theta_coarse == 1) {
        std::vector<Face> faces; faces.reserve(mesh.nFaces());
        for (Face f:mesh.faces()) faces.push_back(f); return std::array<std::vector<Face>,2>({ {}, faces });
    }
    if (mesh.nFaces() == 0)
        return {};
    std::vector<Face> faces;
    faces.reserve(mesh.nFaces());
    double total_res = 0;
    for (Face f : mesh.faces()) {
        faces.push_back(f);
        total_res += residual[f];
    }

    std::ranges::sort(faces, [&residual](const Face &a, const Face &b) -> bool { return residual[a] * residual[a] < residual[b] * residual[b]; });

    std::vector<Face> mark_refine, mark_coarse;
    mark_refine.reserve(conf.theta_refine * mesh.nFaces());
    mark_coarse.reserve(conf.theta_refine * mesh.nFaces());
    double accum_res = 0;
    for (int i = faces.size() - 1; i >= 0; --i) {
        if (residual[faces[i]] > conf.threshold_refine && accum_res <= conf.theta_refine * total_res)
            mark_refine.push_back(faces[i]);
        else break;
        accum_res += residual[faces[i]];
    }
    accum_res = 0;
    for (int i = 0; i < faces.size(); ++i) {
        if (residual[faces[i]]< conf.threshold_coarse && accum_res <= conf.theta_coarse*total_res)
            mark_coarse.push_back(faces[i]);
        else break;
        accum_res += residual[faces[i]];
    }
    return { mark_refine, mark_coarse };
}

DoeflerConf DoerflerPreset(DoerflerPresetConf preset) {
    DoeflerConf conf;
    switch (preset) {
    case DoerflerPresetConf::LOW: // Conservative refinement
        conf.theta_refine     = 0.3;
        conf.threshold_refine = 1e-4;
        conf.theta_coarse     = 0.3;
        conf.threshold_coarse = 1e-6;
        break;
    case DoerflerPresetConf::MEDIUM: // Conservative refinement
        conf.theta_refine     = 0.3;
        conf.threshold_refine = 1e-5;
        conf.theta_coarse     = 0.3;
        conf.threshold_coarse = 1e-7;
        break;
    case DoerflerPresetConf::HIGH: // Conservative refinement
        conf.theta_refine     = 0.3;
        conf.threshold_refine = 1e-6;
        conf.theta_coarse     = 0.3;
        conf.threshold_coarse = 1e-8;
        break;
    case DoerflerPresetConf::VERY_HIGH: // Conservative refinement
        conf.theta_refine     = 0.3;
        conf.threshold_refine = 1e-8;
        conf.theta_coarse     = 0.3;
        conf.threshold_coarse = 1e-10;
        break;
    case DoerflerPresetConf::UNIFORM_REFINE:
        conf.theta_refine     = 1.0;   // force refine
        conf.threshold_refine = 0.0;   // no threshold
        conf.theta_coarse     = 0.0;   // disable coarsening
        conf.threshold_coarse = 0.0;   // disable coarsening
        break;
    case DoerflerPresetConf::UNIFORM_COARSE:
        conf.theta_refine     = 0.0;   // disable refining
        conf.threshold_refine = 1e9;   // effectively prevent refining
        conf.theta_coarse     = 1.0;   // force coarsening
        conf.threshold_coarse = 1e9;   // no threshold
        break;
    }
    return conf;
}

IncrementingIndex::IncrementingIndex(ManifoldSurfaceMesh &mesh) : idx(mesh,invalidIdx) {}

std::size_t &IncrementingIndex::operator[](Face f) {
    if (idx[f] == invalidIdx) idx[f] = current++;
    return idx[f];
}

void IncrementingIndex::compress() {
    std::vector<Face> faces;
    faces.reserve(idx.getMesh()->nVertices());
    for(Face f: idx.getMesh()->faces()){
        faces.push_back(f);
    }
    std::ranges::sort(faces,[&](Face f1, Face f2) { return this->operator[](f1) < this->operator[](f2); } );
    for (int i = 0; i < faces.size(); ++i) {
        idx[faces[i]] = i;
    }
}
}