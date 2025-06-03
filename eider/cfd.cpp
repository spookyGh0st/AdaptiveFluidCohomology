
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/utilities/vector2.h"
#include "poisson.h"
using namespace geometrycentral::surface;
using namespace geometrycentral;

FaceData<Vector2> velocity(
    SurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
    const VertexData<double>& w, const std::vector<double>& c,const std::vector<FaceData<Vector2>>& h)
{
    VertexData<double> f(mesh,0);
    solve_poisson_dirichlet_zero_mean(mesh,geom,f,w);
    FaceData<Vector2> u;

    throw std::runtime_error("Not implemented");

    assert(c.size() == h.size());
    for (std::size_t i = 0; i < h.size(); i++) {
        for (Face f: mesh.faces()) {
            u[f] = u[f] + c[i] * h[i][f];
        }
    }
    return u;
}