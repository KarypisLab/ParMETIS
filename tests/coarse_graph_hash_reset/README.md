# Coarse-graph hash reset reproducer

`repro.c` calls `ParMETIS_V3_PartMeshKway` on a three-way split of a tetrahedral
mesh. `newtonian.connectivity` contains 11,878 four-node elements with
zero-based node IDs. It was extracted from
`SimVascular/svMultiPhysics/tests/cases/fluid/newtonian/mesh/mesh-complete.mesh.vtu`.
The reproducer and fixture are distributed under the BSD-3-Clause license in
`LICENSE.txt`; this directory contains no generated binaries.

Build ParMETIS at the revision being tested. From the ParMETIS repository root,
with MPI, METIS, and GKlib already available:

```sh
cmake -S . -B build-repro -DCMAKE_C_COMPILER=mpicc \
  -DGKLIB_PATH=/path/to/GKlib -DMETIS_PATH=/path/to/METIS
cmake --build build-repro --target parmetis
mpicc -std=gnu99 -Iinclude -Ilibparmetis \
  -I/path/to/METIS/include -I/path/to/GKlib/include \
  tests/coarse_graph_hash_reset/repro.c \
  build-repro/libparmetis/libparmetis.a \
  /path/to/METIS/lib/libmetis.a /path/to/GKlib/lib/libGKlib.a \
  -lm -o build-repro/coarse-graph-hash-reset
mpirun -np 3 build-repro/coarse-graph-hash-reset \
  tests/coarse_graph_hash_reset/newtonian.connectivity
```

The precise METIS/GKlib library names and locations depend on how those
dependencies were built. The reported before/after run used upstream ParMETIS
source with the METIS and GKlib static libraries vendored by svMultiPhysics.
On macOS, `mpirun --oversubscribe -np 3` may be needed.

Before the fix, METIS reports unequal weights for opposite directed edges
(for example, `16 -> 53 = 7` versus `53 -> 16 = 6`), then asserts in
`libmetis__SetupGraph`. With the fix, all three ranks return status 1
(`METIS_OK`) and edgecut 268. The reproducer uses valid zero-based node IDs;
no solver executable, Python package, or VTK file is needed to run it.
