#pragma once

#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <vector>

namespace geometrycentral::surface {
struct Singular_Circle {
    Singular_Circle(ManifoldSurfaceMesh &mesh, const std::vector<Halfedge> &circle);
    Singular_Circle() = default;
    EdgeData<Halfedge> next;
    Edge start_e;
};
// Iterator definition
struct Singular_Circle_Iterator {
    const Singular_Circle &circle;
    Halfedge current;
    bool started = false;

    Singular_Circle_Iterator(const Singular_Circle &c, Halfedge start);

    Halfedge operator*() const;

    Singular_Circle_Iterator &operator++();

    bool operator!=(const Singular_Circle_Iterator &other);
};

Singular_Circle_Iterator begin(const Singular_Circle &circle);

Singular_Circle_Iterator end(const Singular_Circle &circle);
} // namespace geometrycentral::surface
