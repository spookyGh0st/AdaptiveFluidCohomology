#pragma once

#include <geometrycentral/utilities/vector2.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>

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
            if (!he.isInterior()) continue;
            Vertex vi = he.vertex();
            Vertex vj = he.twin().vertex();

            double cotAlpha = geom.halfedgeCotanWeights[he];
            sum += cotAlpha * (f[vj] - f[vi]);
        }

        return sum / geom.vertexDualAreas[v];
    }
}
