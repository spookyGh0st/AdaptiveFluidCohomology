#pragma once

#include "geometrycentral/surface/common_subdivision.h"
#include "geometrycentral/surface/transfer_functions.h"
#include "poisson.h"
#include "refine.h"
#include "cfd.h"
#include "homology.h"

using namespace geometrycentral::surface;

void onSplit(Edge e, Halfedge he1, Halfedge he2, HalfedgeData<std::optional<bool>> &nextLeft);

void onCollapse(Halfedge he, HalfedgeData<std::optional<bool>> &nextLeft);
