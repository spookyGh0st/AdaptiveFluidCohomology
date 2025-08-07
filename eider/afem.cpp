#include "afem.h"

void onSplit_boundary_e(EdgeData<Halfedge> &next, Edge e, Halfedge he, Edge *start_e) {
    if (!he.isInterior())
        return;
    Halfedge out_he = next[e];
    next[e] = Halfedge();
    if (next[he.next().edge()] != Halfedge()) {
        next[he.edge()] = out_he;
    } else {
        Halfedge new_he = he.next().next().twin();
        Halfedge side_he = new_he.next();
        if (next[side_he.edge()].edge().isBoundary()) {
            // Path leads to boundary
            next[side_he.edge()] = new_he.twin();
            next[new_he.edge()] = he.twin();
            next[he.edge()] = out_he;
        } else {
            // TODO: If boundary, connect other end to new he;
            next[he.edge()] = new_he;
            next[new_he.edge()] = out_he;
            assert(out_he == new_he.next().twin());
        }
    }

    // The path now always goes over he.
    // If we started at e, we now start at he.edge()
    if (*start_e == e)
        *start_e = he.edge();
}

void onSplit_marked_e(EdgeData<Halfedge> &next, Edge e, Halfedge he1, Halfedge he2, Edge *start_e) {
    if (e.isBoundary()) {
        onSplit_boundary_e(next, e, he1, start_e);
        onSplit_boundary_e(next, e, he2, start_e);
        return;
    }
    Halfedge out_e = next[e];
    next[e] = Halfedge(); // clear data;
    std::array<Halfedge, 4> cyc_edges{
        he1.next(),
        he1.next().next().twin().next(),
        he2.next(),
        he2.next().next().twin().next(),
    };
    Halfedge in, out;
    for (Halfedge be : cyc_edges) {
        if (next[be.edge()] == Halfedge())
            continue;
        if (next[be.edge()].edge() == e)
            in = be;
        else if (be == out_e.twin())
            out = be;
    }
    assert(in != Halfedge());
    assert(out != Halfedge());
    assert(in != out);

    Halfedge inner_he;
    // 4 cases:
    if (in.tailVertex() == out.tipVertex()) {
        inner_he = in.next().next().twin();
        next[in.edge()] = inner_he;
        next[inner_he.edge()] = out_e;
    } else if (in.tipVertex() == out.tailVertex()) {
        inner_he = in.next().twin();
        next[in.edge()] = inner_he;
        next[inner_he.edge()] = out_e;
    } else {
        inner_he = in.next().twin();
        Halfedge nn = in.next().twin().next().next();
        next[in.edge()] = inner_he;
        next[inner_he.edge()] = nn.twin();
        next[nn.edge()] = out_e;
    }

    // In all cases the path goes over inner_he
    // Update Start HE, in case the split edge was it
    if (e == *start_e) {
        *start_e = inner_he.edge();
    }
}

void onSplit_straight_face(EdgeData<Halfedge> &next, Halfedge he) {
    if (!he.isInterior())
        return;
    Halfedge bhe1 = he.next();
    Halfedge bhe2 = he.next().next().twin().next();
    if (next[bhe1.edge()] == Halfedge()) {
        assert(next[bhe2.edge()] == Halfedge());
    } else if (next[bhe1.edge()] == bhe2.twin()) {
        next[bhe1.edge()] = bhe1.next().twin();
        next[bhe1.next().edge()] = bhe2.twin();
    } else {
        assert(next[bhe2.edge()] == bhe1.twin());
        next[bhe2.edge()] = bhe2.prevOrbitFace().twin();
        next[bhe2.prevOrbitFace().edge()] = bhe1.twin();
    }
}

void onSplit_unmarked_e(Edge e, Halfedge he1, Halfedge he2, EdgeData<Halfedge> &next) {
    onSplit_straight_face(next, he1);
    onSplit_straight_face(next, he2);
}

void onSplit(Edge e, Halfedge he1, Halfedge he2, EdgeData<Halfedge> &next, Edge *start_e) {
    if (next[e] == Halfedge())
        onSplit_unmarked_e(e, he1, he2, next);
    else
        onSplit_marked_e(next, e, he1, he2, start_e);
}
