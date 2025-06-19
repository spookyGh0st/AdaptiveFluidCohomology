#pragma once

#include <geometrycentral/utilities/vector2.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/manifold_surface_mesh.h>

namespace geometrycentral::surface
{
    inline Vector2 grad(IntrinsicGeometryInterface& geom, Face face, const VertexData<double>& f)
    {
        Vector2 g = Vector2::zero();
        for (Halfedge he : face.adjacentHalfedges())
        {
            Vector2 efp = geom.halfedgeVectorsInFace[he.next()];
            g += (-efp.rotate90()) * f[he.vertex()];
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
