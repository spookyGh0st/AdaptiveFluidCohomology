#include "afem.h"
#include <span>

void connect_split_path(std::span<Halfedge> b_halfedges, HalfedgeData<std::optional<bool>>&nextLeft) {
    Halfedge in, out;
    for (Halfedge b_he: b_halfedges){
        GC_SAFETY_ASSERT(!b_he.isDead(), "all adjacent hes should be alive");
        GC_SAFETY_ASSERT(b_he.isInterior(), "all adjacent hes should be interior");
        if (nextLeft[b_he].has_value()) in = b_he;
        else if (nextLeft[b_he.twin()].has_value()) out = b_he;
    }
    if (in == Halfedge() && out == Halfedge()) return;
    GC_SAFETY_ASSERT(in != Halfedge(), "in  must be set")
    GC_SAFETY_ASSERT(out != Halfedge(), "out  must be set")

    // Three cases, this should also work for the boundary case
    if(in.nextLeft().face() == out.face()){
        nextLeft[in] = true;
        nextLeft[in.nextLeft()] = true;
    } else if (in.nextRight().face() == out.face()){
        nextLeft[in] = false;
        nextLeft[in.nextRight()] = false;
    } else {
        nextLeft[in] = true;
        nextLeft[in.nextLeft()] = false;
        nextLeft[in.nextLeft().nextRight()] = true;
        GC_SAFETY_ASSERT(in.nextLeft().nextRight().nextLeft() == out.twin(), "Path malformed");
    }
}

void onSplit(Edge e, Halfedge he1, Halfedge he2, HalfedgeData<std::optional<bool>> &nextLeft) {
    if(!e.isBoundary()){
        for (Halfedge he: he1.vertex().outgoingHalfedges()) {
            nextLeft[he] = {};
            nextLeft[he.twin()] = {};
        }
        Halfedge b1 = he1.next();
        Halfedge b2 = b1.next().twin().next();
        Halfedge b3 = b2.next().twin().next();
        Halfedge b4 = b3.next().twin().next();
        assert(b4.next().twin().next() == b1);
        std::array<Halfedge, 4> cyc_edges{ b1,b2,b3,b4 };
        connect_split_path(cyc_edges,nextLeft);
    } else{
        // TODO: clear previos data
        Halfedge he = (e == he1.edge())? he1 : he2; // Abuse bug in splitBoundaryEdge
        assert(he.isInterior() && he.twin().isInterior());
        std::array<Halfedge, 2> cyc_edges{
            he.twin().prevOrbitFace(),
            he.next(),
        };
        connect_split_path(cyc_edges,nextLeft);
    }
}


void connect_collapsed_path(std::span<Halfedge> b_halfedges, HalfedgeData<std::optional<bool>>&nextLeft) {
    Halfedge in, out;
    for (Halfedge b_he: b_halfedges){
        GC_SAFETY_ASSERT(!b_he.isDead(), "all adjacent hes should be alive");
        GC_SAFETY_ASSERT(b_he.isInterior(), "all adjacent hes should be interior");
        if (nextLeft[b_he].has_value()) in = b_he;
        else if (nextLeft[b_he.twin()].has_value()) out = b_he;
    }
    if (in == Halfedge() && out == Halfedge()) return;
    GC_SAFETY_ASSERT(in != Halfedge(), "in  must be set")
    GC_SAFETY_ASSERT(out != Halfedge(), "out  must be set")

    if(in.nextLeft() == out.twin()){
        nextLeft[in] = true;
    } else if (in.nextRight() == out.twin()){
        nextLeft[in] = false;
    } else {
        auto twin = out.twin();
        if (in.nextLeft().nextRight() == twin) {
            nextLeft[in] = true;
            nextLeft[in.nextLeft()] = false;
        } else if (in.nextRight().nextLeft() == twin) {
            nextLeft[in] = false;
            nextLeft[in.nextRight()] = true;
        } else if (in.nextLeft().nextLeft() == twin) {
            nextLeft[in] = true;
            nextLeft[in.nextLeft()] = true;
        } else if (in.nextRight().nextRight() == twin) {
            nextLeft[in] = false;
            nextLeft[in.nextRight()] = false;
        }
    }
}

void onCollapse(Halfedge he, HalfedgeData<std::optional<bool>> &next) {
    GC_SAFETY_ASSERT(!he.isDead(), "he must be allive");
    if (he.edge().isBoundary()){
        if (!he.isInterior()) he = he.twin();
        std::array<Halfedge,2> b_halfedges = { he.next(), he.next().next(), };
        connect_collapsed_path(b_halfedges,next);
    } else {
        std::array<Halfedge,4> b_halfedges = {
            he.next(), he.next().next(),
            he.twin().next(), he.twin().next().next()
        };
        connect_collapsed_path(b_halfedges,next);
    }
}
