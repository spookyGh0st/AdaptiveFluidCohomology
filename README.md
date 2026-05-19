Project overview

The source code includes the library 'eider' implementing AFEM for two dimensions 
as well as a test framework for it.


The code itself is of scientific nature and not production ready.
// list missing features.

## Obtain the project

```bash

```

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=DebugWithRelease
cmake --build build
```

## Run

To see all registered test use
```shell
./build/test/eider-test --gtest_list_tests
```

Most of the tests are visual tests, asserting that specific parts work as expected.
Perhaps of most interesting are the following tests:
-  afemTest, AdaptiveFluidCohomologyHo

A specific test can be run similar to  
```shell
./build/test/eider-test --gtest_filter=InteractiveTestName.*
```

TODO: First refactor, rename files,


eider
AdaptiveFluidSolver
poisson.h

transfer:

homotopy.h
singular_homology.h
homology.h
AdaptiveHomologyBasis.h

- fjdsk
- jfdlsk


test