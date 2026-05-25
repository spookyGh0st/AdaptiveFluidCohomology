## Project Structure

The code is organized into several key components that correspond to the major theoretical sections:

```
.
├── CMakeLists.txt              # Build configuration
├── cfd.h/.cpp                  # Core fluid dynamics logic (RHS evaluation, time-stepping)
├── homology.h/.cpp             # Homology computation pipeline (delta forms, projection, interpolation)
├── homotopy.h/.cpp             # Homotopy basis computation (Tree-Cotree algorithm)
├── poisson.h/.cpp              # Poisson solver for the stream-function
├── refine.h/.cpp               # Adaptive mesh refinement and coarsening logic
├── transfer.h/.cpp             # Data transfer between different mesh resolutions
├── singular_homology.h/.cpp    # Data structure for representing homology cycles
├── AdaptiveHomologyBasis.h/.cpp# Manages the adaptive homology basis
├── AdaptiveFluidSolver.h/.cpp  # Main class tying all components together
└── util.h                      # Utility functions (norms, gradients, etc.)
```

## Core Concepts and Implementation

### 1. Homology and Harmonic Bases

The fluid solver requires a basis of harmonic 1-forms, which are derived from the topology of the mesh.
The library implements the full pipeline described in the paper to compute these forms.

*   **Theory**: The computation follows the chain of isomorphisms:
    `Homotopy Basis` → `Homology Representatives` → `Closed 1-Form (Delta Form)` → `Harmonic 1-Form (Projection)` → `Orthonormal Basis`.

*   **Implementation**:
    *   **Homotopy Basis**: The process starts in `homotopy.h/.cpp`. The function `greedy_homotopy_basis()` implements the greedy Tree-Cotree algorithm to find a basis for the fundamental group. It uses a Disjoint Set Union (DSU) data structure and Dijkstra's algorithm on the dual mesh.
    *   **Homology Representation**: The resulting homotopy cycles are converted into `Singular_Circle` objects, defined in `singular_homology.h/.cpp`. A `Singular_Circle` cleverly represents a cycle not as a list of mesh elements (which would be invalidated by mesh mutations) but as a `HalfedgeData<std::optional<bool>>` map. This map stores the turning direction at each half-edge, making the cycle representation robust to refinement.
    *   **Harmonic Form Computation**: This pipeline is implemented in `homology.h/.cpp`.
        *   `delta_form()` converts a `Singular_Circle` into a discrete, closed 1-form (`EdgeData<double>`).
        *   `PressureProjectionSolver` solves the least-squares problem to project this closed form onto the space of harmonic forms, ensuring it is also co-closed (`Δη = 0`).
        *   `whitney_interpolation()` converts the resulting harmonic 1-form (on edges) into a piecewise-constant vector field on faces (`FaceData<Vector2>`).
        *   `orthonormalize()` uses a modified Gram-Schmidt process to produce an orthonormal basis.

### 2. Fluid Solver

*   **Theory**: The core of the solver involves two steps: reconstructing the velocity from the vorticity and harmonic components, and then evolving the vorticity and harmonic coefficients over time. 
* These correspond to `alg:velocity` and `alg:evalRHS` in the paper.

*   **Implementation**: The logic resides in `cfd.h/.cpp`.
    *   **Velocity Reconstruction**: The `velocity()` function implements `alg:velocity`. It first solves the Poisson equation `−Δψ = ω` using the `StreamFunctionSolver` (from `poisson.h/.cpp`), then computes the velocity `u = ∇ψ + ∑cᵢhᵢ`.
    *   **Time Evolution**: The `evalRHS()` function implements `alg:evalRHS`. It computes the time derivatives `dω/dt` and `dcᵢ/dt` using the discrete operators for the directional derivative (`derive()`) and the Lamb vector (`Lamb()`), which are defined in `util.h`.

### 3. Adaptive Mesh and Time Stepping

*   **Theory**: The method uses an AFEM loop (solve-estimate-mark-refine) for spatial adaptivity and an embedded Runge-Kutta method (DOPRI5) for temporal adaptivity.

*   **Implementation**:
    *   **Spatial Adaptivity (`refine.h/.cpp`)**:
        *   The `AdaptiveTriangulation` class manages the intrinsic triangulation.
        *   `poisson_residual_error_sqr()` computes the a posteriori error estimator `η_R` for the stream-function Poisson solve.
        *   `select_doerfler()` implements the Dörfler marking strategy to select faces for refinement and coarsening based on the residual.
        *   `refine()` and `coarse()` orchestrate the mesh mutations by calling `vertex_bisection()` and `vertex_biunion()`.
    *   **Temporal Adaptivity (`cfd.h/.cpp`)**:
        *   `adaptive_step()` implements the Dormand-Prince 5(4) method (`alg:Dopri5Step`). It computes two solutions of different orders to estimate the local error and adjusts the next time step `dt` accordingly.

### 4. Adaptive Data Transfer & Topology Preservation

When the mesh is refined or coarsened, all simulation data must be transferred. 

*   **Theory**: Vorticity is transferred using an L2-optimal projection. The harmonic basis is not interpolated directly.
* Instead, the underlying singular cycles that generate the basis are updated during mesh mutations. 
* The new harmonic basis is then recomputed from these updated cycles, guaranteeing it remains a valid discretization of the same continuous harmonic form.

*   **Implementation**:
    *   **Scalar/Vector Field Transfer (`transfer.h/.cpp`)**: The classes `AdaptiveVertexTransfer` and `AdaptiveFaceTransfer` build interpolation matrices (`P_A`, `P_B`) on-the-fly to solve the L2-minimization problem for transferring data between the old (`A`) and new (`B`) meshes via a common subdivision (`S`).
    *   **Harmonic Basis Adaptivity (`AdaptiveHomologyBasis.h/.cpp`)**: 
        *   The `AdaptiveHomologyBasis` class holds the `Homology_basis` (a vector of `Singular_Circle`s).
        *   It registers callbacks with the `IntrinsicTriangulation`. The `onSplit()` and `onCollapse()` functions are invoked during mesh mutations to locally update the `Singular_Circle`'s half-edge map, preserving its homology class.
        *   After the mesh is adapted, `harmonicBasis()` is called. This function re-runs the pipeline from `delta_form()` onwards to generate a new, high-quality harmonic basis on the new mesh.

### 5. Main Solver Loop

*   **Theory**: The final adaptive algorithm combines all the above components into a single time-stepping loop.

*   **Implementation**: The `AdaptiveFluidSolver` class in `AdaptiveFluidSolver.h/.cpp` orchestrates the entire simulation. Its `step()` method performs one time step, which involves:
    1.  Calling `adapt()` to refine/coarsen the mesh and update all data.
    2.  Calling `adaptive_step()` (or `RK4Step` for static time stepping) to advance the simulation state `(ω, c)`.
    3.  Updating the internal state `wc` and time `dt`.

## Building

The project uses submodules for dependencies and lfs for results and videos.
Use the following commands to clone it (replace with https url if needed)

```
git lfs install
git clone --recurse-submodules git@gitlab.informatik.uni-bonn.de:abtb0/master-thesis.git
```

The project uses a standard CMake configuration. 

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build --parallel 8 
```

The source code includes the library 'eider' implementing AFEM for two dimensions 
as well as a test framework for it.
The code itself is of scientific nature and not production ready.


## Run

To see all registered test use
```shell
./build/test/eider-test --gtest_list_tests
```

Most of the tests are visual tests, asserting that specific parts work as expected.
Perhaps of most interesting are the following tests:
- `AdaptiveFluidCohomology.Main` The Main application
- `homologyTest.*` Each test tests one step of the homology pipeline
- `EvaluatorTest.*` Creates the evaluation results

A specific test can be run similar to  
```shell
./build/test/eider-test --gtest_filter=AdaptiveFluidCohomology.Main
```

## Generating Results

This project includes a suite of tests designed to generate video outputs from different fluid solver configurations. 
A custom workflow is available to run these tests and then combine the resulting videos into a side-by-side comparison using FFmpeg.

The project uses a custom CMake target to automate running the specific tests that produce video files. This target will:
1.  Ensure the test executable (`eider-test`) is compiled and up-to-date.
2.  Run the test suite with a filter for `VideoTest.*`.

To generate the individual videos (`evaluation_s.mp4`, `evaluation_c.mp4`, etc.), run the following command from your build directory:

```bash
cmake --build build --target generate_videos
```
After the command completes successfully, the individual video files will be created in the videos directory

Similarly, to generate the results, execute 
```bash
cmake --build build --target generate_evaluation
```
