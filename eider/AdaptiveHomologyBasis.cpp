#include "AdaptiveHomologyBasis.h"
#include "homotopy.h"
#include <span>

void connect_split_path(Halfedge in, Halfedge out, HalfedgeData<std::optional<bool>>&nextLeft) {
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

void onsplit_boundary(Halfedge he, HalfedgeData<std::optional<bool>> &nextLeft) {
    GC_SAFETY_ASSERT(he.edge().isBoundary(),"Always have he be the inner halfedge on the boundary to be split");
    std::array<Halfedge, 4> b_halfedges{
        he,
        he.next(),
        he.next().next().twin().next(),
        he.next().next().twin().next().next(),
    };
    Halfedge in, out;
    for (Halfedge b_he: b_halfedges){
        GC_SAFETY_ASSERT(!b_he.isDead(), "all adjacent hes should be alive");
        GC_SAFETY_ASSERT(b_he.isInterior(), "all adjacent hes should be interior");
        if (nextLeft[b_he].has_value()) in = b_he;
        else if (nextLeft[b_he.twin()].has_value()) out = b_he;
    }
    if (in == Halfedge() && out == Halfedge()) return;
    if (in != Halfedge() && out != Halfedge()) {
        return connect_split_path(in,out,nextLeft);
    }
    if (out == Halfedge()) {
        if (in.next().edge().isBoundary()) {
            nextLeft[in] = false;
            nextLeft[in.next().twin()] = true;
        }else {
            GC_SAFETY_ASSERT(in.prevOrbitFace().edge().isBoundary(),"other way")
            nextLeft[in] = true;
            nextLeft[in.prevOrbitFace().twin()] = true;
        }
    } else {
        if (out.next().edge().isBoundary()) {
            nextLeft[out.next()] = true;
        }else {
            GC_SAFETY_ASSERT(out.prevOrbitFace().edge().isBoundary(),"other way")
            nextLeft[out.prevOrbitFace()] = false;
        }
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
        Halfedge in, out;
        for (Halfedge b_he: cyc_edges){
            GC_SAFETY_ASSERT(!b_he.isDead(), "all adjacent hes should be alive");
            GC_SAFETY_ASSERT(b_he.isInterior(), "all adjacent hes should be interior");
            if (nextLeft[b_he].has_value()) in = b_he;
            else if (nextLeft[b_he.twin()].has_value()) out = b_he;
        }
        if (in == Halfedge() && out == Halfedge()) return;
        GC_SAFETY_ASSERT(in != Halfedge() && out != Halfedge(), "Both in and out need to be set");
        connect_split_path(in, out,nextLeft);
    } else{
        Halfedge he = he1.isInterior()? he1 : he2;
        for (Halfedge he: he.vertex().outgoingHalfedges()) {
            nextLeft[he] = {};
            nextLeft[he.twin()] = {};
        }
        onsplit_boundary(he,nextLeft);
    }
}


void onCollapse_boundary(Halfedge he,HalfedgeData<std::optional<bool>>&nextLeft ) {
    GC_SAFETY_ASSERT(he.edge().isBoundary(),"Always have he be the inner halfedge on the boundary to be split");
    std::array<Halfedge, 3> b_halfedges{ he, he.next(), he.next().next() };
    Halfedge in, out;
    for (Halfedge b_he: b_halfedges){
        GC_SAFETY_ASSERT(!b_he.isDead(), "all adjacent hes should be alive");
        GC_SAFETY_ASSERT(b_he.isInterior(), "all adjacent hes should be interior");
        if (nextLeft[b_he].has_value()) in = b_he;
        else if (nextLeft[b_he.twin()].has_value()) out = b_he;
    }
    if (in == Halfedge() && out == Halfedge()) return;
    if (in != Halfedge() && out != Halfedge()) { // path runs along split edge;
        if (in.nextLeft() == out.twin()) {
            nextLeft[in] = true;
        } else {
            assert(in.nextRight() == out.twin());
            nextLeft[in] = false;
        }
        return;
    }
    if (out == Halfedge()) {
        if (!in.nextLeft().isInterior()) { nextLeft[in] = true; nextLeft[in.nextLeft()] = true; } // boundary is on left
        else { assert(!in.nextRight().isInterior()); nextLeft[in] = false; nextLeft[in.nextRight()] = true;} // on right
    } else {
        if (out.prevOrbitFace().edge().isBoundary()) { nextLeft[out.prevOrbitFace()] = false;; }
        else { assert(out.next().edge().isBoundary()); nextLeft[out.next()] = true; }
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
        onCollapse_boundary(he,next);
    } else {
        std::array<Halfedge,4> b_halfedges = {
            he.next(), he.next().next(),
            he.twin().next(), he.twin().next().next()
        };
        connect_collapsed_path(b_halfedges,next);
    }
}
AdaptiveHomologyBasis::AdaptiveHomologyBasis(IntrinsicTriangulation &icit) : mesh(*icit.intrinsicMesh), geom(icit) {
    Face b_face = arbitrary_base_face(mesh);
    auto homotopy = greedy_homotopy_basis(mesh, geom,b_face);
    homologyB = singular_homology_basis(mesh, homotopy);

    for (auto &b : homologyB) {
        icit.edgeSplitCallbackList.emplace_back([&b](Edge e, Halfedge he1, Halfedge he2) {
          onSplit(e, he1, he2, b.nextLeft); });
    }
    for (auto &b : homologyB) {
        icit.edgeCollapseCallbackList.emplace_back([&b](Halfedge he) { onCollapse(he, b.nextLeft); });
    }
}

Harmonic_basis AdaptiveHomologyBasis::harmonicBasis() const {
    return orthonormal_hom_basis(mesh,geom);
}
