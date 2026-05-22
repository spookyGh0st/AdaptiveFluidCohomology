#pragma once

#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/utilities/vector2.h>

namespace geometrycentral::surface {

inline double L2Norm(FaceData<Vector2> d, IntrinsicGeometryInterface &geom) {
    if (!d.toVector().allFinite())
        return std::numeric_limits<double>::quiet_NaN();
    double s = 0;
    for (Face f : geom.mesh.faces()) {
        s += d[f].norm2() * geom.faceAreas[f];
    }
    return std::sqrt(s);
}
inline double L2NormSqr(VertexData<double> d, IntrinsicGeometryInterface &geom) {
    if (!d.toVector().allFinite())
        return std::numeric_limits<double>::quiet_NaN();
    double s = 0;
    geom.requireVertexDualAreas();
    for (Vertex v : geom.mesh.vertices()) {
        double area = geom.vertexDualAreas[v]; // dual area around vertex
        s += d[v] * d[v] * area;
    }
    geom.unrequireVertexDualAreas();
    return s;
}
inline double L2Norm(VertexData<double> d, IntrinsicGeometryInterface &geom) {
    if (!d.toVector().allFinite())
        return std::numeric_limits<double>::quiet_NaN();
    return std::sqrt(L2NormSqr(d, geom));
}

inline double integral(VertexData<double> d, IntrinsicGeometryInterface &geom) {
    if (!d.toVector().allFinite())
        return std::numeric_limits<double>::quiet_NaN();
    double s = 0;
    for (Face f : geom.mesh.faces()) {
        double fs = 0;
        for (Vertex v : f.adjacentVertices()) {
            fs += d[v];
        }
        s += fs / 3 * geom.faceAreas[f];
    }
    return s;
}

/**
 * Compute the gradient of a scalar function f over a triangle face by
 * \( \nabla f = \frac{1}{2A} \sum_{i=0}^2 f_i \cdot R_{90}(e_i^*),\)
 * where \( f_i \) is the value at vertex \( i \), \( e_i^* \) is the edge vector
 * opposite vertex \( i \), \( R_{90} \) is a 90° counter-clockwise rotation,
 * and \( A \) is the triangle area.
 *
 * @param geom Intrinsic geometry data, requires halfedgeVectorsInFace and faceAreas
 * @param face The triangle face.
 * @param f Scalar function on vertices.
 * @return Gradient of f over the face.
 */
inline Vector2 grad(IntrinsicGeometryInterface &geom, Face face, const VertexData<double> &f) {
    assert(face.isTriangle());
    Vector2 g = Vector2::zero();
    for (Halfedge he : face.adjacentHalfedges()) {
        Vector2 efp = geom.halfedgeVectorsInFace[he.next()];
        g += (efp.rotate90()) * f[he.vertex()];
    }
    return g / (2 * geom.faceAreas[face]);
}

/**
 * Compute the scalar curl of a vector field u around a vertex v.
 *
 * The scalar curl is approximated as:
 *   curl(v) = (1 / A) * sum_over_faces ( u_f · t_e ) * |e|,
 * where:
 *   - u_f is the vector field in face f,
 *   - t_e is the unit tangent vector along edge e in counter-clockwise order around v,
 *   - |e| is the length of the edge,
 *   - A is the dual area at vertex v (e.g. circumcentric or barycentric).
 *
 * This mimics integrating circulation around the vertex and normalizing by area.
 *
 * @param geom Intrinsic geometry (requires halfedgeVectorsInFace, edgeLengths, vertexDualAreas)
 * @param v The vertex at which to compute the curl.
 * @param u Vector field defined per face.
 * @return Scalar curl at the vertex.
 */
inline double curl(IntrinsicGeometryInterface &geom, Vertex v, const FaceData<Vector2> &u) {
    double sum = 0.0;

    for (Halfedge he : v.outgoingHalfedges()) {
        if (!he.isInterior())
            continue;

        Face f = he.face();
        Vector2 u_f = u[f]; // vector field in face

        Vector2 edgeVec = geom.halfedgeVectorsInFace[he]; // edge from v to neighbor
        double len = edgeVec.norm();

        Vector2 tangent = edgeVec.unit().rotate90(); // CCW tangent vector along edge (boundary of face)
        sum += dot(u_f, tangent) * len;
    }

    double area = geom.vertexDualAreas[v]; // should be precomputed, e.g., barycentric or Voronoi dual
    return sum / area;
}

inline double laplacian(IntrinsicGeometryInterface &geom, Vertex v, const VertexData<double> &f) {
    double sum = 0.0;

    for (Halfedge he : v.outgoingHalfedges()) {
        Vertex vi = he.tailVertex();
        Vertex vj = he.tipVertex();
        double wij = geom.edgeCotanWeights[he.edge()];
        sum += wij * (f[vi] - f[vj]);
    }
    return sum;
}

inline double derive(Vertex p, FaceData<Vector2> &u, IntrinsicGeometryInterface &geom, const VertexData<double> &f) {
    double a = 0, s = 0;
    for (Face face : p.adjacentFaces()) {
        s += geom.faceAreas[face] * dot(u[face], grad(geom, face, f));
        a += geom.faceAreas[face];
    }
    return (1 / a) * s;
}

inline Vector2 Lamb(Face face, const VertexData<double> &w, const FaceData<Vector2> &u) {
    double s = 0;
    for (Vertex v : face.adjacentVertices()) {
        s += w[v];
    }
    s = s / 3;

    return -s * u[face].rotate90();
};

inline double diameter(IntrinsicGeometryInterface &geom, Face f) {
    assert(f.isTriangle());
    double d = 0;
    for (Edge e : f.adjacentEdges())
        d = std::max(d, geom.edgeLengths[e]);
    return d;
}

inline double diameter(IntrinsicGeometryInterface &geom, Edge e) {
    return geom.edgeLengths[e];
}

inline double max_diameter(ManifoldSurfaceMesh &mesh, IntrinsicGeometryInterface &geom) {
    double d = 0;
    for (Face f : mesh.faces())
        d = std::max(d, diameter(geom, f));
    return d;
}
} // namespace geometrycentral::surface
