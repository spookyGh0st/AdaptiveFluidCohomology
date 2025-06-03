#pragma once

#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>


namespace geometrycentral::surface {
    void solve_poisson_dirichlet_zero_mean(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                                           VertexData<double>& f, const VertexData<double>& g);
    void solve_stream_function(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                                           VertexData<double>& f, const VertexData<double>& g);
}
