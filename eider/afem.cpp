#include "afem.h"

void onSplit_boundary_e(EdgeData<Halfedge>& next, Edge e, Halfedge he) {
    if (!he.isInterior()) return;
    Halfedge out_he = next[e];
    next[e] = Halfedge();
    if (next[he.next().edge()] != Halfedge()) {
        next[he.edge()] = out_he;
    }else {
        Halfedge new_he = he.next().next().twin();
        Halfedge side_he = new_he.next();
        if (next[side_he.edge()].edge().isBoundary()) {
            // Path leads to boundary
            next[side_he.edge()] = new_he.twin();
            next[new_he.edge()] = he.twin();
            next[he.edge()] = out_he;
        }else {
            // TODO: If boundary, connect other end to new he;
            next[he.edge()] = new_he;
            next[new_he.edge()] = out_he;
            assert(out_he == new_he.next().twin());
        }
    }
}

void onSplit_marked_e(EdgeData<Halfedge>& next, Edge e, Halfedge he1, Halfedge he2)
{
    if (e.isBoundary()) {
        onSplit_boundary_e(next,e,he1);
        onSplit_boundary_e(next,e,he2);
        return;
    }
    Halfedge out_e = next[e];
    next[e] = Halfedge(); // clear data;
    std::array<Halfedge, 4> cyc_edges {
        he1.next(),
        he1.next().next().twin().next(),
        he2.next(),
        he2.next().next().twin().next(),
      };
    Halfedge in, out;
    for (Halfedge be: cyc_edges) {
        if (next[be.edge()] == Halfedge()) continue;
        if (next[be.edge()].edge() == e) in = be;
        else if (be == out_e.twin()) out = be;
    }
    assert(in != Halfedge()); assert(out != Halfedge()); assert(in != out);

    // 4 cases:
    if (in.tailVertex() == out.tipVertex())
    {
        next[in.edge()] = in.next().next().twin();
        next[in.next().next().edge()] = out_e;
    } else if (in.tipVertex() == out.tailVertex()) {
        next[in.edge()] = in.next().twin();
        next[in.next().edge()] = out_e;
    } else {
        next[in.edge()] = in.next().twin();
        const Halfedge nn = in.next().twin().next().next();
        next[in.next().edge()] = nn.twin();
        next[nn.edge()] =  out_e;
    }
}

void onSplit_straight_face(EdgeData<Halfedge>& next, Halfedge he)
{
    if (!he.isInterior()) return;
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

void onSplit_unmarked_e(EdgeData<Halfedge>& next, Edge e, Halfedge he1, Halfedge he2) {
    onSplit_straight_face(next,he1);
    onSplit_straight_face(next,he2);
}

void onSplit(EdgeData<Halfedge>& next, Edge e, Halfedge he1, Halfedge he2)
{
    if (next[e] == Halfedge())
        onSplit_unmarked_e(next,e,he1,he2);
    else
        onSplit_marked_e(next,e,he1,he2);
}
