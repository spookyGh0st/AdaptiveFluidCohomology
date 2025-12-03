#pragma  once
#include <eider/AdaptiveFluidSolver.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

namespace geometrycentral::surface{

using MeshP = std::unique_ptr<ManifoldSurfaceMesh>;
using GeomP = std::unique_ptr<VertexPositionGeometry>;

VertexData<Vector3> intrinsic_geom(IntrinsicTriangulation& Tri, VertexPositionGeometry& inputG);

std::pair<MeshP, GeomP > uniform_refine(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, int refine,MARKING_STRATEGY strategy = MARKING_STRATEGY::PATTERN);

std::pair<MeshP, GeomP > delauny_refine(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, double degree = 25, double h_size = std::numeric_limits<double>::infinity());

extern AdaptiveFluidSolverData static_solver_data;
}
