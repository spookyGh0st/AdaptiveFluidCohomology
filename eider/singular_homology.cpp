#include "singular_homology.h"
namespace geometrycentral::surface {
Singular_Circle::Singular_Circle(ManifoldSurfaceMesh &mesh, const std::vector<Halfedge> &circle)
    : nextLeft(mesh) {
    for (std::size_t i = 0; i < circle.size(); ++i) {
        Halfedge he = circle[i];
        Halfedge nhe = circle[(i + 1) % circle.size()];
        // if we are at the boundary, we set this to an arbitrary value
        if (!he.isInterior()) {
            nextLeft[he] = true;
            continue;
        }
        // otherwise to either left or false
        if (he.next().twin() == nhe) {
            nextLeft[he] = false;
        } else {
            assert(he.next().next().twin() == nhe);
            nextLeft[he] = true;
        }
    }
}
Singular_Circle_Iterator::Singular_Circle_Iterator(const Singular_Circle &c,
                                                   Halfedge start)
    : circle(c), current(start), start(start) {}
Halfedge Singular_Circle_Iterator::operator*() const {
    return current;
}
Singular_Circle_Iterator &Singular_Circle_Iterator::operator++() {
    auto nl = circle.nextLeft[current];
    assert(nl.has_value());
    if (!current.isInterior())
        current = start;
    else {
        if (*nl)
            current = current.next().next().twin();
        else
            current = current.next().twin();
    }
    assert(current != Halfedge());
    return *this;
}
bool Singular_Circle_Iterator::operator!=(
    const Singular_Circle_Iterator &other) {
    if (current == other.current && started)
        return false;
    started = true;
    return true;
}
Halfedge startHe(const Singular_Circle &circle) {
    const auto &m = circle.nextLeft.getMesh();
    if (m->hasBoundary()) {
        for (const auto &b : m->boundaryLoops()) {
            for (Halfedge he : b.adjacentHalfedges()) {
                assert(he.twin().isInterior());
                if (circle.nextLeft[he.twin()].has_value())
                    return he.twin();
            }
        }
    }

    Halfedge start_e;
    for (Halfedge he : circle.nextLeft.getMesh()->halfedges()) {
        if (circle.nextLeft[he].has_value())
            return he;
    }
    throw std::runtime_error("This should never be reached");
}
Singular_Circle_Iterator begin(const Singular_Circle &circle) {
    return Singular_Circle_Iterator(circle, startHe(circle));
}
Singular_Circle_Iterator end(const Singular_Circle &circle) {
    return Singular_Circle_Iterator(circle, startHe(circle));
}
}; // namespace geometrycentral::surface
