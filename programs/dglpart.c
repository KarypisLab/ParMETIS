/*
 * dglpart.c
 * 
 * This file partitions a graph for use in DistDGL
 *
 * Started 12/31/2020
 * George
 *
 */

#include <parmetisbin.h>

#define NCON    5

/*************************************************************************/
/*! Entry point of the partitioning code */
/*************************************************************************/
int main(int argc, char *argv[])
{
  idx_t mype, npes;
  MPI_Comm comm;

  MPI_Init(&argc, &argv);
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);
  gkMPI_Comm_size(comm, &npes);
  gkMPI_Comm_rank(comm, &mype);

  if (argc != 3) {
    if (mype == 0)
      printf("Usage: %s <fstem> <nparts-per-pe>\n", argv[0]);

    MPI_Finalize();
    exit(0);
  }

  DistDGL_GPart(argv[1], atoi(argv[2]), comm);

  gkMPI_Comm_free(&comm);

  MPI_Finalize();

  return 0;
}



/***********************************************************************************/
/*! Partition, move, and save the local graphs */
/***********************************************************************************/
void DistDGL_GPart(char *fstem, idx_t nparts_per_pe, MPI_Comm comm)
{
  idx_t i, npes, mype, nparts;
  graph_t *graph, *mgraph;
  idx_t *part, *cpart, *mpart;
  idx_t numflag=0, wgtflag=0, options[10], edgecut, ndims;
  real_t ipc2redist, *xyz=NULL, *tpwgts=NULL, *ubvec=NULL;

  gkMPI_Comm_size(comm, &npes);
  gkMPI_Comm_rank(comm, &mype);

  /* read and create the graph */
  graph = DistDGL_ReadGraph(fstem, comm);
  gkMPI_Barrier(comm);


  /*======================================================================
  / Partition the graph 
  /=======================================================================*/
  options[0] = 1;
  options[1] = 15;
  options[2] = 1;
  wgtflag = 2;
  numflag = 0;
  edgecut = 0;

  nparts = npes*nparts_per_pe;
  part   = imalloc(graph->nvtxs, "DistDGL_GPart: part");
  tpwgts = rsmalloc(nparts*graph->ncon, 1.0/(real_t)nparts, "DistDGL_GPart: tpwgts");
  ubvec  = rsmalloc(graph->ncon, 1.05, "DistDGL_GPart: unvec");

  if (mype == 0)
    printf("\nDistDGL partitioning, ncon: %"PRIDX", nparts: %"PRIDX"\n", graph->ncon, nparts);

  ParMETIS_V3_PartKway(graph->vtxdist, graph->xadj, graph->adjncy, graph->vwgt, NULL, 
      &wgtflag, &numflag, &(graph->ncon), &nparts, tpwgts, ubvec, options, &edgecut, 
      part, &comm);


  /*======================================================================
  / Move the graph based on the partitioning
  /=======================================================================*/

  /* "coarsen" to create an npes-level partitioning and save the original
     partitioning information in the vwgt */
  cpart = imalloc(graph->nvtxs, "DistDGL_GPart: cpart");
  for (i=0; i<graph->nvtxs; i++) {
    graph->vwgt[i] = part[i];  
    cpart[i] = part[i]/nparts_per_pe; 
  }
  graph->ncon = 1; /* so it will just move the part[] info */

  mgraph = DistDGL_MoveGraph(graph, cpart, comm);

  /* "uncoarsen" the partition into mpart */
  mpart = imalloc(mgraph->nvtxs, "DistDGL_GPart: mpart");
  icopy(mgraph->nvtxs, mgraph->vwgt, mpart);

  gk_free((void **)&tpwgts, &ubvec, &part, &cpart, &mpart, LTERM);
}





/******************************************************************************
* This function takes a graph and its partition vector and creates a new
* graph corresponding to the one after the movement
*******************************************************************************/
graph_t *DistDGL_MoveGraph(graph_t *ograph, idx_t *part, MPI_Comm comm)
{
  idx_t npes, mype;
  ctrl_t *ctrl;
  graph_t *graph, *mgraph;
  idx_t options[5] = {0, 0, 1, 0, 0};

  gkMPI_Comm_size(comm, &npes);
  ctrl = SetupCtrl(PARMETIS_OP_KMETIS, NULL, 1, npes, NULL, NULL, comm); 
  mype = ctrl->mype;

  ctrl->CoarsenTo = 1;  /* Needed by SetUpGraph, otherwise we can FP errors */
  graph = SetUpGraph(ctrl, 1, ograph->vtxdist, ograph->xadj, ograph->vwgt, ograph->vsize,
              ograph->adjncy, ograph->adjwgt, 0);
  AllocateWSpace(ctrl, 0);

  CommSetup(ctrl, graph);
  graph->where = part;
  graph->ncon  = 1;

  mgraph = MoveGraph(ctrl, graph);

  graph->where = NULL;
  FreeInitialGraphAndRemap(graph);
  FreeCtrl(&ctrl);

  return mgraph;
}


