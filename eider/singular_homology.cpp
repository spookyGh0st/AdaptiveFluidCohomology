#include "singular_homology.h"
namespace geometrycentral::surface {
Singular_Circle::Singular_Circle(ManifoldSurfaceMesh &mesh, const std::vector<Halfedge> &circle) : next(mesh, Halfedge()), start_e(circle.begin()->edge()) {
    for (std::size_t i = 0; i < circle.size(); ++i) {
        next[circle[i].edge()] = circle[(i + 1) % circle.size()];
    }
}
Singular_Circle_Iterator::Singular_Circle_Iterator(const Singular_Circle &c,
                                                   Halfedge start)
    : circle(c), current(start) {}
Halfedge Singular_Circle_Iterator::operator*() const {
    return current;
}
Singular_Circle_Iterator &Singular_Circle_Iterator::operator++() {
    current = circle.next[current.edge()];
    return *this;
}
bool Singular_Circle_Iterator::operator!=(
    const Singular_Circle_Iterator &other) {
    if (current == other.current && started)
        return false;
    started = true;
    return true;
}
Singular_Circle_Iterator begin(const Singular_Circle &circle) {
    Halfedge start = circle.next[circle.start_e];
    return Singular_Circle_Iterator(circle, start);
}
Singular_Circle_Iterator end(const Singular_Circle &circle) {
    Halfedge start = circle.next[circle.start_e];
    return Singular_Circle_Iterator(circle, start);
}
}; // namespace geometrycentral::surface
