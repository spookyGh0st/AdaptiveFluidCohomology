#pragma once

#include <geometrycentral/utilities/vector2.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/manifold_surface_mesh.h>

namespace geometrycentral::surface
{
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
    inline Vector2 grad(IntrinsicGeometryInterface& geom, Face face, const VertexData<double>& f)
    {
        assert(face.isTriangle());
        Vector2 g = Vector2::zero();
        for (Halfedge he : face.adjacentHalfedges())
        {
            Vector2 efp = geom.halfedgeVectorsInFace[he.next()];
            g += (efp.rotate90()) * f[he.vertex()];
        }
        return g / (2 * geom.faceAreas[face]);
    }

    inline double laplacian(IntrinsicGeometryInterface& geom, Vertex v, const VertexData<double>& f)
    {
        double sum = 0.0;

        for (Halfedge he : v.outgoingHalfedges()) {
            Vertex vi = he.tailVertex();
            Vertex vj = he.tipVertex();
            double wij =geom.edgeCotanWeights[he.edge()];
            sum += wij * (f[vi] - f[vj]);
        }
        return sum;
    }

    inline double diameter(IntrinsicGeometryInterface& geom, Face f) {
        assert(f.isTriangle());
        double d = 0;
        for (Edge e: f.adjacentEdges())
            d = std::max(d,geom.edgeLengths[e]);
        return d;
    }

    inline double diameter(IntrinsicGeometryInterface& geom, Edge e) {
        return geom.edgeLengths[e];
    }

    inline double max_diameter(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom) {
        double d = 0;
        for (Face f: mesh.faces())
            d = std::max(d,diameter(geom, f));
        return d;
    }
}
