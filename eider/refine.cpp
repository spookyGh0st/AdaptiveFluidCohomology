#include "refine.h"

#include <unordered_set>
#include "util.h"

namespace geometrycentral::surface {
void refine(IntrinsicTriangulation &tri, std::vector<Face> faces) {
    auto &M = tri.intrinsicMesh;

    FaceData<Halfedge> refinement_edges(*M);
    tri.requireEdgeLengths();
    for (Face f : M->faces()) {
        double maxEl = std::numeric_limits<double>::lowest();
        for (Halfedge he : f.adjacentHalfedges()) {
            if (tri.edgeLengths[he.edge()] > maxEl) {
                maxEl = tri.edgeLengths[he.edge()];
                refinement_edges[f] = he;
            }
        }
    }
    tri.unrequireEdgeLengths();

    ManifoldSurfaceMesh &mesh = *tri.intrinsicMesh;
    std::unordered_set<Face> marked_faces;
    std::unordered_set<Edge> start_edges;
    // marked_edges.reserve(faces.size()*2);
    while (!faces.empty()) {
        Face f = faces.back();
        faces.pop_back();
        marked_faces.insert(f);
        Halfedge he = refinement_edges[f];
        if (!he.edge().isBoundary() && !marked_faces.contains(he.twin().face())) {
            faces.push_back(he.twin().face());
        }
        if (he.edge().isBoundary() || refinement_edges[he.twin().face()] == he.twin()) {
            start_edges.insert(he.edge());
        }
    }

    while (!start_edges.empty()) {
        Edge e = *start_edges.begin();
        start_edges.erase(start_edges.begin());
        marked_faces.erase(e.halfedge().face());
        marked_faces.erase(e.halfedge().twin().face());

        Vertex new_v = tri.splitEdge(e.halfedge(), 0.5).vertex();

        // Update refinement edges
        for (Halfedge he_o : new_v.outgoingHalfedges()) {
            if (!he_o.isInterior())
                continue;
            refinement_edges[he_o.face()] = he_o.next();
        }

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
            if (refinement_edges[side_f].edge() == side_he.edge())
                start_edges.insert(side_he.edge());
        }
    }

    assert(marked_faces.empty());
    tri.refreshQuantities();
}

Halfedge coarse_halfedge(Vertex v) {
    assert((!v.isBoundary() && v.faceDegree() == 4) || (v.isBoundary() && v.faceDegree() == 2));
    size_t k = 0;
    std::size_t min_idx = std::numeric_limits<std::size_t>::max();
    Halfedge he;
    for (Halfedge ohe : v.outgoingHalfedges()) {
        if (!ohe.isInterior())
            continue;
        std::size_t f_idx = ohe.face().getIndex();
        if (f_idx < min_idx) {
            min_idx = f_idx;
            he = ohe;
        };
    }
    assert(he != Halfedge());
    return he.twin().next();
}

inline bool vertex_is_good(IntrinsicTriangulation &T, Vertex v) {
    if (v == Vertex())
        return false;
    if (T.vertexLocations[v].type == SurfacePointType::Vertex)
        return false;
    return (v.isBoundary() && v.degree() == 2) || (!v.isBoundary() && v.degree() == 4);
}

void coarse(IntrinsicTriangulation &m, const std::function<bool(Vertex)> &f) {
    std::unordered_set<Vertex> good_vertices{};
    for (Vertex v : m.intrinsicMesh->vertices()) {
        if (vertex_is_good(m, v) && f(v))
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
        Vertex cv = m.intrinsicMesh->collapseEdgeTriangular(he);
        assert(cv != Vertex());

        for (Vertex nv : adjacent_v) {
            if (vertex_is_good(m, nv) && f(nv))
                good_vertices.insert(nv);
        }
    }
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