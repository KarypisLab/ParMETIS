/* SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others. */
/* SPDX-License-Identifier: BSD-3-Clause */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "parmetislib.h"

static void fail(const char *message, MPI_Comm comm)
{
  fprintf(stderr, "%s\n", message);
  MPI_Abort(comm, 1);
}

int main(int argc, char **argv)
{
  int rank, nprocs, nelem = 0, nodes_per_elem = 0;
  int *counts, *displacements, *root_connectivity = NULL;
  idx_t *elmdist, *eptr, *eind, *part;
  idx_t wgtflag = 0, numflag = 0, ncon = 1, ncommonnodes = 3;
  idx_t nparts, edgecut, options[10] = {0};
  real_t *tpwgts, ubvec[1] = {1.05};
  MPI_Comm comm;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if (argc != 2 || nprocs < 2 || sizeof(idx_t) != sizeof(int))
    fail("usage: mpirun -np N repro connectivity.txt (32-bit idx_t)", MPI_COMM_WORLD);

  if (rank == 0) {
    FILE *input = fopen(argv[1], "r");
    if (input == NULL || fscanf(input, "%d %d", &nelem, &nodes_per_elem) != 2)
      fail("cannot read mesh dimensions", MPI_COMM_WORLD);
    root_connectivity = malloc((size_t)nelem * nodes_per_elem * sizeof(int));
    if (root_connectivity == NULL)
      fail("cannot allocate mesh connectivity", MPI_COMM_WORLD);
    for (int i = 0; i < nelem * nodes_per_elem; ++i) {
      if (fscanf(input, "%d", &root_connectivity[i]) != 1)
        fail("cannot read mesh connectivity", MPI_COMM_WORLD);
    }
    fclose(input);
  }

  MPI_Bcast(&nelem, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nodes_per_elem, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (nelem < nprocs || nodes_per_elem != 4)
    fail("expected at least one tetrahedron per rank", MPI_COMM_WORLD);

  elmdist = malloc((size_t)(nprocs + 1) * sizeof(idx_t));
  counts = malloc((size_t)nprocs * sizeof(int));
  displacements = malloc((size_t)nprocs * sizeof(int));
  tpwgts = malloc((size_t)nprocs * sizeof(real_t));
  if (!elmdist || !counts || !displacements || !tpwgts)
    fail("cannot allocate partition arrays", MPI_COMM_WORLD);

  for (int i = 0; i < nprocs; ++i) {
    elmdist[i] = (idx_t)round((double)i * nelem / nprocs);
    tpwgts[i] = (real_t)1.0 / nprocs;
  }
  elmdist[nprocs] = nelem;
  for (int i = 0; i < nprocs; ++i) {
    displacements[i] = elmdist[i] * nodes_per_elem;
    counts[i] = (elmdist[i + 1] - elmdist[i]) * nodes_per_elem;
  }

  int local_elements = elmdist[rank + 1] - elmdist[rank];
  eptr = malloc((size_t)(local_elements + 1) * sizeof(idx_t));
  eind = malloc((size_t)counts[rank] * sizeof(idx_t));
  part = malloc((size_t)local_elements * sizeof(idx_t));
  if (!eptr || !eind || !part)
    fail("cannot allocate local mesh", MPI_COMM_WORLD);

  MPI_Scatterv(root_connectivity, counts, displacements, MPI_INT,
               eind, counts[rank], MPI_INT, 0, MPI_COMM_WORLD);
  for (int i = 0; i <= local_elements; ++i)
    eptr[i] = i * nodes_per_elem;

  options[0] = 1;
  options[PMV3_OPTION_DBGLVL] = 0;
  options[PMV3_OPTION_SEED] = 10;
  nparts = nprocs;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);
  int status = ParMETIS_V3_PartMeshKway(elmdist, eptr, eind, NULL, &wgtflag,
      &numflag, &ncon, &ncommonnodes, &nparts, tpwgts, ubvec,
      options, &edgecut, part, &comm);
  MPI_Comm_free(&comm);
  printf("rank %d: status=%d edgecut=%d local_elements=%d\n",
         rank, status, (int)edgecut, local_elements);

  free(root_connectivity);
  free(counts);
  free(displacements);
  free(elmdist);
  free(eptr);
  free(eind);
  free(part);
  free(tpwgts);
  MPI_Finalize();
  return status == METIS_OK ? 0 : 1;
}
