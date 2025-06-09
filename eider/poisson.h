#pragma once

#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>


namespace geometrycentral::surface {
    void solve_poisson_dirichlet_zero_mean(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                                           VertexData<double>& f, const VertexData<double>& g);

    struct StreamFunctionSolver {

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        Eigen::SparseMatrix<double> M_II, L_IB;
        VertexData<std::size_t> globalToInteriorIndex, globalToBoundaryIndex; // global to interior/boundary indices
        std::vector<Vertex> boundaryVertices {}, interiorVertices {};

        void compute(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom);

        void solve_dirichlet(VertexData<double> & f, const VertexData<double>& g) const;
        void solve(SurfaceMesh & mesh, IntrinsicGeometryInterface & geom, VertexData<double> & f, const VertexData<double>& g) const;
    };
}
