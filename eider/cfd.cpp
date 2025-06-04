#include "cfd.h"

#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/utilities/vector2.h"
#include "poisson.h"

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

    inline double derive(Vertex p, FaceData<Vector2>& u, IntrinsicGeometryInterface& geom, const VertexData<double>& f)
    {
        double a = 0, s = 0;
        for (Face face : p.adjacentFaces())
        {
            s += geom.faceAreas[face] * dot(u[face], grad(geom, face, f));
            a += geom.faceAreas[face];
        }
        return (1 / a) * s;
    }


    inline Vector2 Lamb(Face face, const VertexData<double>& w, const FaceData<Vector2>& u)
    {
        double s = 0;
        for (Vertex v : face.adjacentVertices())
        {
            s += w[v];
        }
        s = s / 3;

        return -s * u[face].rotate90();
    };


    wc_wrapper operator+(const wc_wrapper& a, const wc_wrapper& b)
    {
        assert(a.c.size() == b.c.size());

        wc_wrapper result;
        result.w = a.w + b.w;
        result.c.resize(a.c.size());
        for (size_t i = 0; i < a.c.size(); ++i)
            result.c[i] = a.c[i] + b.c[i];
        return result;
    }

    wc_wrapper operator*(double s, const wc_wrapper& wc)
    {
        wc_wrapper result;
        result.w = wc.w * s;
        result.c.resize(wc.c.size());
        for (size_t i = 0; i < wc.c.size(); ++i)
            result.c[i] = wc.c[i] * s;
        return result;
    }

    FaceData<Vector2> velocity(
        SurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
        const wc_wrapper& wc, const std::vector<FaceData<Vector2>>& h)
    {
        VertexData<double> f(mesh, 0);
        solve_stream_function(mesh, geom, f, wc.w);
        FaceData<Vector2> u(mesh, Vector2::zero());
        for (Face face : mesh.faces())
        {
            u[face] = -grad(geom, face, f).rotate90();
        }

        assert(wc.c.size() == h.size());
        for (std::size_t i = 0; i < h.size(); i++)
            for (Face face : mesh.faces())
                u[face] = u[face] + wc.c[i] * h[i][face];
        return u;
    }


    wc_wrapper evalRHS(
        SurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
        const wc_wrapper& wc, const std::vector<FaceData<Vector2>>& h)
    {
        geom.requireFaceAreas(); geom.requireHalfedgeVectorsInFace();
        auto u = velocity(mesh, geom, wc, h);

        auto l = FaceData<Vector2>(mesh, Vector2::zero());
        for (Face f : mesh.faces()) { l[f] = Lamb(f, wc.w, u); }

        VertexData<double> dw(mesh, 0);
        for (Vertex v : mesh.vertices())
        {
            dw[v] = -derive(v, u, geom, wc.w);
        }
        std::vector<double> dc(wc.c.size(), 0);
        for (std::size_t i = 0; i < wc.c.size(); i++)
            for (Face f : mesh.faces())
                dc[i] += dot(l[f], h[i][f]) * geom.faceAreas[f];

        geom.unrequireFaceAreas(); geom.unrequireHalfedgeVectorsInFace();
        return {dw, dc};
    }


    wc_wrapper RK4Step(
        SurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
        const std::vector<FaceData<Vector2>>& h,
        const wc_wrapper& x, double dt
    )
    {
        auto F = [&mesh, &geom, &h](const wc_wrapper& wc) -> wc_wrapper { return evalRHS(mesh, geom, wc, h); };
        wc_wrapper k1 = F(x), k2 = F(x + dt / 2 * k1);
        wc_wrapper k3 = F(x + dt / 2 * k2), k4 = F(x + dt * k3);
        return x + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4);
    }


}


