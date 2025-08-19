#include "refine.h"

#include <unordered_set>
#include "util.h"

namespace geometrycentral::surface {

AdaptiveTriangulation::AdaptiveTriangulation(IntrinsicTriangulation &tri): tri(tri), idx(*tri.intrinsicMesh), marked_corner(*tri.intrinsicMesh,false) {
    tri.requireEdgeLengths();
    for (Face f : tri.intrinsicMesh->faces()) {
        double maxEl = std::numeric_limits<double>::lowest();
        Corner c;
        for (Halfedge he : f.adjacentHalfedges()) {
            if (tri.edgeLengths[he.edge()] > maxEl) {
                maxEl = tri.edgeLengths[he.edge()];
                c = he.oppositeCorner();
            }
        }
        marked_corner[c] = true;
    }
    tri.unrequireEdgeLengths();
}

Halfedge AdaptiveTriangulation::vertex_bisection(Halfedge he) {
    assert (he.isInterior());
    assert(marked_corner[he.oppositeCorner()]);
    if (he.twin().isInterior()) { assert(marked_corner[he.twin().oppositeCorner()]); }

    Halfedge twin_he = he.twin();
    Vertex  tip_vertex = he.tipVertex();

    // TODO: does this return the right he?
    he = tri.splitEdge(he, 0.5);
    // ensure indices are kept in order
    assert(twin_he.tipVertex() == he.tailVertex());
    assert(he.tipVertex() == tip_vertex );
    if (he.isInterior()) {
        Face fl = he.prevOrbitFace().twin().face(), fr = he.face();
        if (idx[fl] > idx[fr]) std::swap(idx[fl],idx[fr]);
    }
    assert(twin_he.tailVertex() == tip_vertex );
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

void AdaptiveTriangulation::refine(std::vector<Face> faces) {
    std::unordered_set<Face> marked_faces;
    std::unordered_set<Edge> start_edges;
    // marked_edges.reserve(faces.size()*2);
    while (!faces.empty()) {
        Face f = faces.back();
        faces.pop_back();
        marked_faces.insert(f);
        Halfedge he = getRefinementEdge(f);
        if (!he.edge().isBoundary() && !marked_faces.contains(he.twin().face())) {
            faces.push_back(he.twin().face());
        }
        if (he.edge().isBoundary() || getRefinementEdge(he.twin().face()) == he.twin()) {
            start_edges.insert(he.edge());
        }
    }

    while (!start_edges.empty()) {
        Edge e = *start_edges.begin();
        start_edges.erase(start_edges.begin());

        marked_faces.erase(e.halfedge().face());
        marked_faces.erase(e.halfedge().twin().face());

        Halfedge split_he = e.halfedge().isInterior() ? e.halfedge() : e.halfedge().twin();
        Vertex new_v = vertex_bisection(split_he).vertex();


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
            if (!marked_faces.contains(side_f))
                continue;
            if (getRefinementEdge(side_f).edge() == side_he.edge())
                start_edges.insert(side_he.edge());
        }
    }

    assert(marked_faces.empty());
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

Halfedge AdaptiveTriangulation::vertex_biunion(Halfedge he) {
    // TODO: Assert left face has smaller idx then right face
    std::size_t l_idx = idx[he.prevOrbitFace().twin().face()];
    std::size_t r_idx = he.twin().isInterior()? idx[he.twin().face()] : -1;

    // Assert all adjacent faces have the refinement edges opposite to this vertex
    for (Halfedge ohe: he.vertex().outgoingHalfedges()) {
        if (!ohe.isInterior()) continue;
        assert(getRefinementEdge(ohe.face()) == ohe.next());
    }

    Vertex vj = he.tipVertex();
    he = tri.collapseEdgeTriangular(he);
    assert(vj == he.tipVertex());

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

void AdaptiveTriangulation::coarse(const std::function<bool(Vertex)> &f) {
    std::unordered_set<Vertex> good_vertices{};
    for (Vertex v : tri.intrinsicMesh->vertices()) {
        if (vertex_is_good(tri, v) && f(v))
            good_vertices.insert(v);
    }

    while (!good_vertices.empty()) {
        // pop element
        Vertex v = *good_vertices.begin();
        good_vertices.erase(v);

        Halfedge he = coarse_halfedge(v);

        std::array<Vertex, 2> adjacent_v;
        adjacent_v[0] = he.next().tipVertex();
        if (he.twin().isInterior())
            adjacent_v[1] = he.twin().next().tipVertex();
        Halfedge cv = vertex_biunion(he);
        assert(he != Halfedge());

        /*
        for (Vertex nv : adjacent_v) {
            if (vertex_is_good(m, nv) && f(nv))
                good_vertices.insert(nv);
        }
         */
    }
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

std::vector<Face> select_doerfler(ManifoldSurfaceMesh &mesh, FaceData<double> residual, double theta, double threshold) {
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

    std::vector<Face> result;
    result.reserve(theta * mesh.nFaces());
    double accum_res = 0;
    for (int i = faces.size() - 1; i >= 0; --i) {
        if (residual[faces[i]] < threshold)
            break;
        result.push_back(faces[i]);
        accum_res += residual[faces[i]];
        if (accum_res >= theta * total_res)
            break;
    }
    return result;
}
}