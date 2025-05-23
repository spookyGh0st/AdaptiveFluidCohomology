#pragma once
#include <vector>
#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral::surface
{
    std::vector<Edge> computeMinimalSpanningTree(SurfaceMesh& mesh, IntrinsicGeometryInterface& geometry);
    std::vector<Edge> computePrimalEdgesOfDualMaxST(SurfaceMesh& mesh, IntrinsicGeometryInterface& geometry);
}

