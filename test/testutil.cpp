#include "testutil.h"

namespace geometrycentral::surface {
AdaptiveFluidSolverData static_solver_data = AdaptiveFluidSolverData(
    DOPRI5PresetConf::MEDIUM,
    DoerflerPresetConf::LOW,
    0.01,
    false,
    false,
    MARKING_STRATEGY::PATTERN,
    false,
    false);

std::pair<MeshP, GeomP> uniform_refine(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geom, int refine, MARKING_STRATEGY strategy) {
    AdaptiveTriangulation atri(mesh, geom, strategy);
    for (int i = 0; i < refine; ++i) {
        std::vector<Face> faces;
        for (Face f : atri.mesh().faces())
            faces.push_back(f);
        atri.refine(faces);
    }
    atri.mesh().compress();
    atri.intrinsicTriangulation().refreshQuantities();
    auto nMesh = atri.mesh().copy();
    auto nGeom = std::make_unique<VertexPositionGeometry>(*nMesh, intrinsic_geom(atri.intrinsicTriangulation(), geom).reinterpretTo(*nMesh));
    return {std::move(nMesh), std::move(nGeom)};
}
std::pair<MeshP, GeomP> delauny_refine(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geom, double degree, double h_size) {
    IntegerCoordinatesIntrinsicTriangulation atri(mesh, geom);
    atri.delaunayRefine(degree);
    atri.intrinsicMesh->compress();
    atri.refreshQuantities();
    auto nMesh = atri.intrinsicMesh->copy();
    auto nGeom = std::make_unique<VertexPositionGeometry>(*nMesh, intrinsic_geom(atri, geom).reinterpretTo(*nMesh));
    return {std::move(nMesh), std::move(nGeom)};
}
VertexData<Vector3> intrinsic_geom(IntrinsicTriangulation &Tri, VertexPositionGeometry &inputG) {
    auto &mesh = *Tri.intrinsicMesh;
    VertexData<Vector3> int_positions(mesh);
    for (Vertex v : mesh.vertices()) {
        int_positions[v] = Tri.vertexLocations[v].interpolate(inputG.vertexPositions);
    }
    return int_positions;
}
} // namespace geometrycentral::surface
