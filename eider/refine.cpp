#include "refine.h"

#include <geometrycentral/surface/surface_mesh.h>
#include <set>

#include "util.h"

namespace gcs=geometrycentral::surface;

void refine(gcs::IntrinsicTriangulation& tri, std::vector<gcs::Face> faces) {
    auto& M = tri.intrinsicMesh;

    gcs::FaceData<gcs::Halfedge> refinement_edges (*M);
    tri.requireEdgeLengths();
    for (gcs::Face f: M->faces()) {
        double maxEl = std::numeric_limits<double>::lowest();
        for (gcs::Halfedge he: f.adjacentHalfedges()) {
            if (tri.edgeLengths[he.edge()] > maxEl) {
                maxEl = tri.edgeLengths[he.edge()];
                refinement_edges[f] = he;
            }
        }
    }
    tri.unrequireEdgeLengths();

    std::set<gcs::Edge> marked_edges;
    std::set<gcs::Edge> start_edges;
    // marked_edges.reserve(faces.size()*2);

    while (!faces.empty())
    {
        gcs::Face f = faces.back(); faces.pop_back();
        gcs::Halfedge he = refinement_edges[f];
        marked_edges.insert(he.edge());
        if (he.edge().isBoundary() || refinement_edges[he.twin().face()] == he.twin()) {
            start_edges.insert(he.edge());
        }else if (!he.edge().isBoundary() && !marked_edges.contains(refinement_edges[he.twin().face()].edge())) {
            faces.push_back(he.twin().face());
        }
    }

    while (!start_edges.empty())
    {
        gcs::Edge e = *start_edges.begin();
        start_edges.erase(start_edges.begin());

        gcs::Halfedge he = e.halfedge();
        gcs::Halfedge n_he = tri.splitEdge(he,0.5);
        refinement_edges[n_he.face()] = n_he.next();
        refinement_edges[he.face()] = he.prevOrbitFace();
        if (!n_he.edge().isBoundary())
        {
            refinement_edges[n_he.twin().face()] = n_he.twin().prevOrbitFace();
            refinement_edges[he.twin().face()] = he.twin().next();
        }

        auto he_r = n_he.next(), he_l = n_he.next().next().twin().next();
        if (he_r.isInterior() && refinement_edges[he_r.twin().face()] == he_r.twin() && marked_edges.contains(he_r.edge())) {
            start_edges.insert(he_r.edge());
        }
        if (he_l.isInterior() && refinement_edges[he_l.twin().face()] == he_l.twin() && marked_edges.contains(he_l.edge())) {
            start_edges.insert(he_l.edge());
        }
        marked_edges.erase(e);
    }

    assert(marked_edges.empty());
    tri.refreshQuantities();
}

using namespace geometrycentral::surface;
using namespace geometrycentral;

double etaR(Face T, IntrinsicGeometryInterface& geom, const VertexData<double>& f, const VertexData<double>& u)
{
    double f_st = 0, h_t = geom.faceAreas[T];
    for (Vertex v: T.adjacentVertices()) f_st = f[v];
    f_st /= 3;

    double lu = 0;
    for (Vertex v: T.adjacentVertices()) lu += laplacian(geom,v,u);
    lu /= 3;

    double jump_sum = 0;
    for (Edge e : T.adjacentEdges())
    {
        if (e.isBoundary()) continue; // only sum over interior edges
        double h_e = geom.edgeLengths[e], j = 0;
        for (Halfedge he : e.adjacentHalfedges()) // for both sides of edge
        {
            j += dot(grad(geom, he.face(), u), geom.halfedgeVectorsInFace[he].rotateCW(PI/2));
        }
        jump_sum += h_e * h_e * j * j;
    }

    return h_t * h_t * std::pow(f_st + lu,2) + jump_sum;
}

FaceData<double> poisson_residual_error(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, const VertexData<double>& f, const VertexData<double>& u)
{
    geom.requireHalfedgeVectorsInFace(); geom.requireEdgeCotanWeights();
    FaceData<double> eta(mesh);
    for (Face T: mesh.faces()) {
        eta[T] = etaR(T,geom,f,u);
    }
    geom.unrequireHalfedgeVectorsInFace(); geom.requireEdgeCotanWeights();
    return eta;
}