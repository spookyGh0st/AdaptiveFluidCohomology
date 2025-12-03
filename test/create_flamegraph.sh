#!/usr/bin/env bash
set -euo pipefail


# Get absolute path of this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Set working directory relative to script location

cd "$SCRIPT_DIR"
cd ..
PATH="$(pwd)/external/FlameGraph:$PATH"
cd "cmake-build-relwithdebinfo/test"

perf record -q -F 997 --call-graph dwarf ./eider-test --gtest_filter=EvaluatorTest.EvaluatePerformance1
perf script > out.perf
stackcollapse-perf.pl out.perf > out.folded

grep 'AdaptiveFluidSolver::step' out.folded > tmp.folded
sed -i 's|eider-test;geometrycentral::surface::EvaluatorTest_EvaluatePerformance1_Test::TestBody;||g' tmp.folded
sed -i 's/geometrycentral::surface:://g' tmp.folded
sed -i 's/AdaptiveFluidSolver:://g' tmp.folded
sed -i 's/AdaptiveHomologyBasis:://g' tmp.folded
sed -i 's/AdaptiveTriangulation:://g' tmp.folded
sed -i 's/AdaptiveVertexTransfer:://g' tmp.folded
sed -i 's/SurfaceMesh:://g' tmp.folded
sed -i 's/BaseGeometryInterface:://g' tmp.folded
sed -i 's/IntrinsicGeometryInterface:://g' tmp.folded
sed -i 's/Eigen:://g' tmp.folded
sed -i 's/PressureProjectionSolver:://g' tmp.folded

max_frames="8"
tmp="$(mktemp)"
awk -v N="$max_frames" '
{
  cnt = $NF
  $NF = ""
  sub(/[ \t]+$/, "")
  stack = $0
  n = split(stack, f, ";")
  out = f[1]
  for (i = 2; i <= N && i <= n; ++i) out = out ";" f[i]
  print out " " cnt
}' tmp.folded > "$tmp"

flamegraph.pl --title "Adaptive, recomputed h" --bgcolors "#ffffff" "$tmp" > out.svg
inkscape -D out.svg  -o "../../../tex/thesis/figures/svg/flamegraph-recomputed.pdf" --export-dpi=600 --export-background=white --export-background-opacity=1 

