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

#define CHUNKSIZE (1<<15)

/* The following is to perform a cyclic distribution of the input vertex IDs
   in order to balance the adjancency lists during the partitioning computations */
/* (u%npes)*lnvtxs + u/npes */
#define ToCyclicMap(id, nbuckets, bucketsize) \
            (((id)%(nbuckets))*(bucketsize) + (id)/(nbuckets))
/*
#define ToCyclicMap(id, nbuckets, bucketdist) \
            ((bucketdist)[(id)%(nbuckets)] + (id)/(nbuckets))
*/
/* (u%lnvtxs)*npes + u/lnvtxs */
#define FromCyclicMap(id, nbuckets, bucketsize) \
            (((id)%(bucketsize))*(nbuckets) + (id)/(bucketsize))

/*! The structure that contains size-information of types of moved data */
typedef struct mvinfo_t {
  idx_t nvtxs;       /*!< The number of vertices to move */ 
  idx_t nedges;      /*!< The number of edges to move */      
  idx_t nvmdata;     /*!< The #of idx_t-equivalent of vertex metadata */
  idx_t nemdata;     /*!< The #of idx_t-equivalent of edge metadata */
} mvinfo_t;


int DistDGL_GPart(char *fstem, idx_t nparts_per_pe, MPI_Comm comm);
graph_t *DistDGL_ReadGraph(char *fstem, MPI_Comm comm);
graph_t *DistDGL_MoveGraph(graph_t *ograph, idx_t *part, idx_t nparts_per_pe, MPI_Comm comm);
void DistDGL_CheckMGraph(ctrl_t *ctrl, graph_t *graph, idx_t nparts_per_pe);
void i2kvsorti(size_t n, i2kv_t *base);
void i2kvsortii(size_t n, i2kv_t *base);
void DistDGL_WriteGraphs(char *fstem, graph_t *graph, idx_t nparts_per_pe, MPI_Comm comm);
idx_t DistDGL_mapFromCyclic(idx_t u, idx_t npes, idx_t *vtxdist);
idx_t DistDGL_mapToCyclic(idx_t u, idx_t npes, idx_t *vtxdist);


/*************************************************************************/
/*! Entry point of the partitioning code */
/*************************************************************************/
int main(int argc, char *argv[])
{
  idx_t mype, npes;
  MPI_Comm comm;
  int retval;

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

  retval = DistDGL_GPart(argv[1], atoi(argv[2]), comm);

  gkMPI_Comm_free(&comm);

  MPI_Finalize();

  return retval;
}


/*************************************************************************/
/*! Partition, move, and save the local graphs */
/*************************************************************************/
int DistDGL_GPart(char *fstem, idx_t nparts_per_pe, MPI_Comm comm)
{
  idx_t i, npes, mype, nparts;
  graph_t *graph, *mgraph;
  idx_t *part;
  idx_t numflag=0, wgtflag=0, options[10], edgecut, ndims;
  real_t *tpwgts=NULL, *ubvec=NULL;

  gkMPI_Comm_size(comm, &npes);
  gkMPI_Comm_rank(comm, &mype);

  /* read and create the graph */
  graph = DistDGL_ReadGraph(fstem, comm);
  gkMPI_Barrier(comm);
  if (graph == NULL)
    return EXIT_FAILURE;

  /* Report peak memory use before partitioning */
  for (i=npes-1; i>=0; i--) {
    if (i == mype && mype == npes-1)
      printf("\n-------------------------------------------------------\n");
    if (i == mype) 
      printf("[%03"PRIDX"] proc/self/stat/VmPeak:     %.2f MB\n", mype, (float)gk_GetProcVmPeak()/(1024.0*1024.0));
    if (i == 0 && mype == 0)
      printf("-------------------------------------------------------\n");
    fflush(stdout);
    gkMPI_Barrier(comm);
  }

  /*======================================================================
  / Partition the graph 
  /=======================================================================*/
  options[0] = 1;
  options[1] = 15 + (PARMETIS_DBGLVL_TWOHOP|PARMETIS_DBGLVL_FAST|PARMETIS_DBGLVL_DROPEDGES|PARMETIS_DBGLVL_ONDISK);
  options[2] = 1;
  wgtflag = 2;
  numflag = 0;
  edgecut = 0;

  nparts = npes*nparts_per_pe;
  part   = imalloc(graph->nvtxs, "DistDGL_GPart: part");
  tpwgts = rsmalloc(nparts*graph->ncon, 1.0/(real_t)nparts, "DistDGL_GPart: tpwgts");
  ubvec  = rsmalloc(graph->ncon, 1.02, "DistDGL_GPart: unvec");

  if (mype == 0)
    printf("\nDistDGL partitioning, ncon: %"PRIDX", nparts: %"PRIDX" [%s, %s, MPI %d.%d]\n", 
        graph->ncon, nparts, 
        (sizeof(idx_t)==8 ? "i64" : "i32"), (sizeof(real_t)==8 ? "r64" : "r32"),
        MPI_VERSION, MPI_SUBVERSION);

  ParMETIS_V3_PartKway(graph->vtxdist, graph->xadj, graph->adjncy, graph->vwgt, NULL, 
      &wgtflag, &numflag, &(graph->ncon), &nparts, tpwgts, ubvec, options, &edgecut, 
      part, &comm);

  /*======================================================================
  / Move the graph based on the partitioning
  /=======================================================================*/
  mgraph = DistDGL_MoveGraph(graph, part, nparts_per_pe, comm);

  /*======================================================================
  / Write the different partitions to disk 
  /=======================================================================*/
  DistDGL_WriteGraphs(fstem, mgraph, nparts_per_pe, comm);

  /* Report peak memory use after partitioning */
  for (i=npes-1; i>=0; i--) {
    if (i == mype && mype == npes-1)
      printf("\n-------------------------------------------------------\n");
    if (i == mype) 
      printf("[%03"PRIDX"] proc/self/stat/VmPeak:     %.2f MB\n", mype, (float)gk_GetProcVmPeak()/(1024.0*1024.0));
    if (i == 0 && mype == 0)
      printf("-------------------------------------------------------\n");
    fflush(stdout);
    gkMPI_Barrier(comm);
  }

  return EXIT_SUCCESS;
}


/*************************************************************************/
/*! This function takes a graph and its partition vector and creates a new
     graph corresponding to the one after the movement */
/*************************************************************************/
graph_t *DistDGL_MoveGraph(graph_t *ograph, idx_t *part, idx_t nparts_per_pe, 
             MPI_Comm comm)
{
  idx_t npes, mype, nparts, idxwidth;
  ctrl_t *ctrl;
  idx_t h, i, ii, j, jj, k, nvtxs, nsnbrs, nrnbrs;
  idx_t *xadj, *adjncy, *mvtxdist;
  idx_t *where, *newlabel, *vtype, *lpwgts, *gpwgts;
  idx_t *sgraph, *rgraph;
  mvinfo_t *sinfo, *rinfo;
  graph_t *graph, *mgraph;
  idx_t *vmptr, *emptr;
  char *vmdata, *emdata;
  idx_t *iptr, ilen; 

  idxwidth = sizeof(idx_t);

  gkMPI_Comm_size(comm, &npes);
  ctrl = SetupCtrl(PARMETIS_OP_KMETIS, NULL, 1, npes, NULL, NULL, comm); 
  npes = ctrl->npes;
  mype = ctrl->mype;

  ctrl->CoarsenTo = 1;  /* Needed by SetUpGraph, otherwise we can FP errors */
  graph = SetupGraph(ctrl, 1, ograph->vtxdist, ograph->xadj, ograph->vwgt, ograph->vsize,
              ograph->adjncy, ograph->adjwgt, 0);
  AllocateWSpace(ctrl, 0);

  CommSetup(ctrl, graph);

  nparts = npes*nparts_per_pe;

  WCOREPUSH;

  /* read the metadata from the disk */
  {
    size_t size;
    char filein[256];

    sprintf(filein, "emdata-%d-%"PRIDX".bin", (int)getpid(), mype);
    ograph->emdata = gk_creadfilebin(filein, &size);
    if (size != ograph->emdata_size)
      printf("[%4"PRIDX"] size: %zu != %zu\n", mype, size, ograph->emdata_size);
    gk_rmpath(filein);
    
    sprintf(filein, "vmdata-%d-%"PRIDX".bin", (int)getpid(), mype);
    ograph->vmdata = gk_creadfilebin(filein, &size);
    if (size != ograph->vmdata_size)
      printf("[%4"PRIDX"] size: %zu != %zu\n", mype, size, ograph->vmdata_size);
    gk_rmpath(filein);
  }

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;
  where  = part;

  vmptr  = ograph->vmptr;
  emptr  = ograph->emptr;
  vmdata = ograph->vmdata;
  emdata = ograph->emdata;
  vtype  = ograph->vtype;

  mvtxdist = imalloc(nparts+1, "DistDGL_MoveGraph: mvtxdist");

  /* Let's do a prefix scan to determine the labeling of the nodes given */
  lpwgts = iwspacemalloc(ctrl, nparts+1);
  gpwgts = iwspacemalloc(ctrl, nparts+1);
  sinfo  = (mvinfo_t *)iwspacemalloc(ctrl, nparts*sizeof(mvinfo_t)/idxwidth);
  rinfo  = (mvinfo_t *)iwspacemalloc(ctrl, nparts*sizeof(mvinfo_t)/idxwidth);

  for (i=0; i<nparts; i++)
    sinfo[i].nvtxs = sinfo[i].nedges = sinfo[i].nvmdata = sinfo[i].nemdata = 0;

  for (i=0; i<nvtxs; i++) {
    sinfo[where[i]].nvtxs   += 1;
    sinfo[where[i]].nedges  += xadj[i+1]-xadj[i];  
    sinfo[where[i]].nvmdata += (vmptr[i+1]-vmptr[i])/idxwidth; /* vmdata */
    for (j=xadj[i]; j<xadj[i+1]; j++) 
      sinfo[where[i]].nemdata += (emptr[j+1]-emptr[j])/idxwidth; /* emdata */
  }
  for (i=0; i<nparts; i++)
    lpwgts[i] = sinfo[i].nvtxs;

  gkMPI_Scan((void *)lpwgts, (void *)gpwgts, nparts, IDX_T, MPI_SUM, ctrl->comm);
  gkMPI_Allreduce((void *)lpwgts, (void *)mvtxdist, nparts, IDX_T, MPI_SUM, ctrl->comm);
  MAKECSR(i, nparts, mvtxdist);


  /* gpwgts[i] will store the label of the first vertex for each domain 
     in each processor */
  for (i=0; i<nparts; i++) 
    /* We were interested in an exclusive scan */
    gpwgts[i] = mvtxdist[i] + gpwgts[i] - lpwgts[i];

  newlabel = iwspacemalloc(ctrl, nvtxs+graph->nrecv);
  for (i=0; i<nvtxs; i++) 
    newlabel[i] = gpwgts[where[i]]++;

  /* Send the newlabel info to processors storing adjacent interface nodes */
  CommInterfaceData(ctrl, graph, newlabel, newlabel+nvtxs);

  /* Tell everybody what and from where they will get it. */
  gkMPI_Alltoall((void *)sinfo, nparts_per_pe*(sizeof(mvinfo_t)/idxwidth), IDX_T, 
                 (void *)rinfo, nparts_per_pe*(sizeof(mvinfo_t)/idxwidth), IDX_T, 
                 ctrl->comm);

  /* Use lpwgts/gpwgts as pointers to where data will be sent/received, respectively */
  for (nsnbrs=0, i=0; i<nparts; i++) {
    lpwgts[i] = 3*sinfo[i].nvtxs 
                + sinfo[i].nedges 
                + sinfo[i].nvtxs 
                + sinfo[i].nvmdata 
                + sinfo[i].nedges 
                + sinfo[i].nemdata
                ;
    if (sinfo[i].nvtxs > 0)
      nsnbrs++;
  }
  MAKECSR(i, nparts, lpwgts);
  //myprintf(ctrl, "lpwgts: %d %d %d\n", lpwgts[0], lpwgts[1], lpwgts[2]);

  /* The target locations are designed to pack in consecutive memory locations 
     the different chunks of the subpartitions that are coming from the 
     different processors. */
  for (nrnbrs=0, k=0; k<nparts_per_pe; k++) {
    for (i=0; i<npes; i++) {
      gpwgts[k*npes+i] = 3*rinfo[i*nparts_per_pe+k].nvtxs 
                         + rinfo[i*nparts_per_pe+k].nedges 
                         + rinfo[i*nparts_per_pe+k].nvtxs 
                         + rinfo[i*nparts_per_pe+k].nvmdata 
                         + rinfo[i*nparts_per_pe+k].nedges 
                         + rinfo[i*nparts_per_pe+k].nemdata
                         ;
      if (rinfo[i*nparts_per_pe+k].nvtxs > 0)
        nrnbrs++;
    }
  }
  MAKECSR(i, nparts, gpwgts);
  //myprintf(ctrl, "gpwgts: %d %d %d\n", gpwgts[0], gpwgts[1], gpwgts[2]);

  /* Update the max # of sreq/rreq/statuses */
  CommUpdateNnbrs(ctrl, gk_max(nsnbrs, nrnbrs));

  rgraph = iwspacemalloc(ctrl, gpwgts[nparts]);
  WCOREPUSH;  /* for freeing the send part early */
  sgraph = iwspacemalloc(ctrl, lpwgts[nparts]);

  /* Issue the receives first */
  for (j=0, i=0; i<nparts; i++) {
    if (rinfo[i].nvtxs > 0) {
      //myprintf(ctrl, "[%"PRIDX"]Irecv from: %"PRIDX" tag: %"PRIDX"\n", i, i%npes, 1+i/nparts_per_pe);
      gkMPI_Irecv((void *)(rgraph+gpwgts[i]), gpwgts[i+1]-gpwgts[i], IDX_T, 
          i%npes, 1+i/npes, ctrl->comm, ctrl->rreq+j++);
    }
    else 
      PASSERT(ctrl, gpwgts[i+1]-gpwgts[i] == 0);
  }

  /* Assemble the graph to be sent and send it */
  for (i=0; i<nvtxs; i++) {
    PASSERT(ctrl, where[i] >= 0 && where[i] < nparts);
    ii = lpwgts[where[i]];
    sgraph[ii++] = xadj[i+1]-xadj[i];
    sgraph[ii++] = where[i];
    sgraph[ii++] = vtype[i];
    for (j=xadj[i]; j<xadj[i+1]; j++)
      sgraph[ii++] = newlabel[adjncy[j]];

    /* vertex metadata */
    sgraph[ii++] = vmptr[i+1]-vmptr[i];
    iptr = (idx_t *)(vmdata+vmptr[i]);
    ilen = (vmptr[i+1]-vmptr[i])/idxwidth;
    for (h=0; h<ilen; h++)
      sgraph[ii++] = iptr[h];

    /* edge metadata */
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      sgraph[ii++] = emptr[j+1]-emptr[j];
      iptr = (idx_t *)(emdata+emptr[j]);
      ilen = (emptr[j+1]-emptr[j])/idxwidth;
      for (h=0; h<ilen; h++)
        sgraph[ii++] = iptr[h];
    }

    lpwgts[where[i]] = ii;
  }
  SHIFTCSR(i, nparts, lpwgts);
  //myprintf(ctrl, "lpwgts: %d %d %d\n", lpwgts[0], lpwgts[1], lpwgts[2]);

  for (j=0, i=0; i<nparts; i++) {
    if (sinfo[i].nvtxs > 0) {
      //myprintf(ctrl, "[%"PRIDX"]Send to: %"PRIDX" tag: %"PRIDX"\n", i, i/nparts_per_pe, 1+i%nparts_per_pe);
      gkMPI_Isend((void *)(sgraph+lpwgts[i]), lpwgts[i+1]-lpwgts[i], IDX_T, 
          i/nparts_per_pe, 1+i%nparts_per_pe, ctrl->comm, ctrl->sreq+j++);
    }
    else 
      PASSERT(ctrl, lpwgts[i+1]-lpwgts[i] == 0);
  }

  /* Wait for the send/recv to finish */
  gkMPI_Waitall(nrnbrs, ctrl->rreq, ctrl->statuses);
  gkMPI_Waitall(nsnbrs, ctrl->sreq, ctrl->statuses);

  WCOREPOP;  /* frees sgraph */

  /* OK, now go and put the graph into graph_t Format */
  mgraph = CreateGraph();
  
  mgraph->vtxdist = mvtxdist;
  mgraph->gnvtxs  = graph->gnvtxs;
  mgraph->ncon    = 1;
  mgraph->level   = 0;
  mgraph->nvtxs   = mgraph->nedges = 0;

  idx_t nvmdata=0, nemdata=0;
  for (i=0; i<nparts; i++) {
    mgraph->nvtxs  += rinfo[i].nvtxs;
    mgraph->nedges += rinfo[i].nedges;
    nvmdata        += rinfo[i].nvmdata;
    nemdata        += rinfo[i].nemdata;
  }

  nvtxs  = mgraph->nvtxs;
  xadj   = mgraph->xadj   = imalloc(nvtxs+1, "MMG: mgraph->xadj");
  adjncy = mgraph->adjncy = imalloc(mgraph->nedges, "MMG: mgraph->adjncy");
  vmptr  = mgraph->vmptr  = imalloc(nvtxs+1, "MMG: mgraph->vmptr");
  emptr  = mgraph->emptr  = imalloc(mgraph->nedges+1, "MMG: mgraph->emptr");
  vmdata = mgraph->vmdata = gk_cmalloc(nvmdata*idxwidth, "MMG: mgraph->vmdata");
  emdata = mgraph->emdata = gk_cmalloc(nemdata*idxwidth, "MMG: mgraph->emdata");
  where  = mgraph->where  = imalloc(nvtxs, "MMG: mgraph->where");
  vtype  = mgraph->vtype  = imalloc(nvtxs, "MMG: mgraph->vtype");

  idx_t *ivptr=(idx_t *)vmdata, *ieptr=(idx_t *)emdata;
  for (jj=ii=i=0; i<nvtxs; i++) {
    xadj[i]  = rgraph[ii++];
    where[i] = rgraph[ii++];
    vtype[i] = rgraph[ii++];
    for (j=0; j<xadj[i]; j++)
      adjncy[jj+j] = rgraph[ii++];

    /* vertex metadata */
    vmptr[i] = rgraph[ii++];
    ilen = vmptr[i]/idxwidth;
    for (h=0; h<ilen; h++, ivptr++)
      *ivptr = rgraph[ii++];

    /* edge metadata */
    for (j=0; j<xadj[i]; j++) {
      emptr[jj+j] = rgraph[ii++];
      ilen = emptr[jj+j]/idxwidth;
      for (h=0; h<ilen; h++, ieptr++)
        *ieptr = rgraph[ii++];
    }
    jj += xadj[i];

    ASSERT(where[i] >= mype*nparts_per_pe && where[i] < (mype+1)*nparts_per_pe);
  }
  MAKECSR(i, nvtxs, xadj);
  MAKECSR(i, nvtxs, vmptr);
  MAKECSR(i, mgraph->nedges, emptr);

  PASSERT(ctrl, jj == mgraph->nedges);
  PASSERT(ctrl, ii == gpwgts[nparts]);
  PASSERTP(ctrl, jj == mgraph->nedges, (ctrl, "%"PRIDX" %"PRIDX"\n", jj, mgraph->nedges));
  PASSERTP(ctrl, ii == gpwgts[nparts], (ctrl, "%"PRIDX" %"PRIDX" %"PRIDX" %"PRIDX" %"PRIDX"\n", 
      ii, gpwgts[nparts], jj, mgraph->nedges, nvtxs));


#ifdef DEBUG
  IFSET(ctrl->dbglvl, DBG_INFO, rprintf(ctrl, "Checking moved graph...\n"));
  DistDGL_CheckMGraph(ctrl, mgraph, nparts_per_pe);
  IFSET(ctrl->dbglvl, DBG_INFO, rprintf(ctrl, "Moved graph is consistent.\n"));
#endif

#ifdef XXX
  /* write the graph in stdout */
  {
    idx_t u, v, i, j, pe;

    for (pe=0; pe<npes; pe++) {
      if (mype == pe) {
        for (i=0; i<mgraph->nvtxs; i++) {
          for (j=mgraph->xadj[i]; j<mgraph->xadj[i+1]; j++) {
            u = mgraph->vtxdist[mype*nparts_per_pe]+i;
            v = mgraph->adjncy[j];

            printf("XXX%"PRIDX" [%"PRIDX" %"PRIDX"] vmdata: [%s]  emdata: [%s]\n", 
                mype,
                u, v, 
                mgraph->vmdata+mgraph->vmptr[i],
                (mgraph->emptr[j+1]-mgraph->emptr[j] == 0 ? "NULL" : mgraph->emdata+mgraph->emptr[j]));
          }
        }
        fflush(stdout);
      }
      gkMPI_Barrier(comm);
    }
  }
#endif

  WCOREPOP;

  graph->where = NULL;
  FreeInitialGraphAndRemap(graph);
  FreeCtrl(&ctrl);

  mgraph->where = where;
  return mgraph;
}


#ifdef XXXX
/*************************************************************************/
/*! This function takes a graph and its partition vector and creates a new
     graph corresponding to the one after the movement */
/*************************************************************************/
DistDGL_TypePermute(graph_t *ograph, idx_t nparts_per_pe, MPI_Comm comm)
{
  idx_t npes, mype, nparts;
  ctrl_t *ctrl;
  idx_t h, i, ii, j, jj, k, nvtxs;
  idx_t *xadj, *adjncy;
  idx_t *where, *newlabel, *vtype;
  graph_t *graph, *mgraph;
  idx_t *vmptr, *emptr;
  char *vmdata, *emdata;
  idx_t *iptr, ilen; 
  i2kv_t *cand;
  idx_t *cvtxdist, *fvtxdist;

  fvtxdist = ograph->vtxdist;
  vmptr    = ograph->vmptr;
  emptr    = ograph->emptr;
  vmdata   = ograph->vmdata;
  emdata   = ograph->emdata;
  vtype    = ograph->vtype;
  where    = ograph->where;

  gkMPI_Comm_size(comm, &npes);
  ctrl = SetupCtrl(PARMETIS_OP_KMETIS, NULL, 1, npes, NULL, NULL, comm); 
  mype = ctrl->mype;

  ctrl->CoarsenTo = 1;  /* Needed by SetUpGraph, otherwise we can FP errors */
  vtxdist = imalloc(npes+1, "DistDGL_TypePermute: vtxdist");
  for (i=0; i<npes; i++)
    cvtxdist[i] = fvtxdist[(i+1)*nparts_per_pe]-fvtxdist[i*nparts_per_pe];
  MAKECSR(i, npes, cvtxdist);

  graph = SetupGraph(ctrl, 1, cvtxdist, ograph->xadj, NULL, NULL, ograph->adjncy, NULL, 0);
  AllocateWSpace(ctrl, 0);

  CommSetup(ctrl, graph);

  nparts = npes*nparts_per_pe;

  WCOREPUSH;

  nvtxs   = graph->nvtxs;
  xadj    = graph->xadj;
  adjncy  = graph->adjncy;

  cand = (i2kv_t *)iwspacemalloc(ctrl, nvtxs*sizeof(i2kv_t)/sizeof(idx_t));
  for (i=0; i<nvtxs; i++) {
    cand[i].key1 = where[i];
    cand[i].key2 = vtype[i];
    cand[i].val  = i;
  }
  i2kvsortii(nvtxs, cand);

  newlabel = iwspacemalloc(ctrl, nvtxs+graph->nrecv);
  for (i=0; i<nvtxs; i++) 
    newlabel[cand[i].key3] = cvtxdist[mype]+i;

  /* Send the newlabel info to processors storing adjacent interface nodes */
  CommInterfaceData(ctrl, graph, newlabel, newlabel+nvtxs);





  WCOREPOP;

  graph->where = NULL;
  FreeInitialGraphAndRemap(graph);
  FreeCtrl(&ctrl);

}
#endif

/*************************************************************************/
/*! Checks the local consistency of moved graph. */
/*************************************************************************/
void DistDGL_CheckMGraph(ctrl_t *ctrl, graph_t *graph, idx_t nparts_per_pe)
{
  idx_t i, j, jj, k, nvtxs, firstvtx, lastvtx;
  idx_t *xadj, *adjncy, *vtxdist;

  nvtxs   = graph->nvtxs;
  xadj    = graph->xadj;
  adjncy  = graph->adjncy;
  vtxdist = graph->vtxdist;

  firstvtx = vtxdist[nparts_per_pe*ctrl->mype];
  lastvtx  = vtxdist[nparts_per_pe*ctrl->mype+1];

  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (firstvtx+i == adjncy[j])
        myprintf(ctrl, "(%"PRIDX" %"PRIDX") diagonal entry\n", i, i);

      if (adjncy[j] >= firstvtx && adjncy[j] < lastvtx) {
        k = adjncy[j]-firstvtx;
        for (jj=xadj[k]; jj<xadj[k+1]; jj++) {
          if (adjncy[jj] == firstvtx+i)
            break;
        }
        if (jj == xadj[k+1])
          myprintf(ctrl, "(%"PRIDX" %"PRIDX") but not (%"PRIDX" %"PRIDX") [%"PRIDX" %"PRIDX"] [%"PRIDX" %"PRIDX"]\n", 
              i, k, k, i, firstvtx+i, firstvtx+k, 
              xadj[i+1]-xadj[i], xadj[k+1]-xadj[k]);
      }
    }
  }
}


/*************************************************************************/
/*! Reads, distributes, and pre-processes the DistDGL's input files to 
    create the graph used for partitioning */
/*************************************************************************/
graph_t *DistDGL_ReadGraph(char *fstem, MPI_Comm comm)
{
  idx_t i, j, pe, idxwidth;
  idx_t npes, mype, ier;
  graph_t *graph=NULL;
  idx_t gnvtxs, gnedges, nvtxs, ncon; 
  idx_t *vtxdist, *xadj, *adjncy, *vwgt, *vtype;
  MPI_Status stat;
  ssize_t rlen, fsize;
  size_t lnlen=0;
  char *filename=NULL, *line=NULL;
  FILE *fpin=NULL;

  idxwidth = sizeof(idx_t);

  gkMPI_Comm_size(comm, &npes);
  gkMPI_Comm_rank(comm, &mype);

  /* check to see if the required files exist */ 
  ier = 0;
  if (mype == 0) {
    filename = gk_malloc(100+strlen(fstem), "DistDGL_ReadGraph: filename");

    sprintf(filename, "%s_stats.txt", fstem);
    if (!gk_fexists(filename)) {
      printf("ERROR: File '%s' does not exists.\n", filename);
      ier++;
    }

    sprintf(filename, "%s_edges.txt", fstem);
    if (!gk_fexists(filename)) {
      printf("ERROR: File '%s' does not exists.\n", filename);
      ier++;
    }

    sprintf(filename, "%s_nodes.txt", fstem);
    if (!gk_fexists(filename)) {
      printf("ERROR: File '%s' does not exists.\n", filename);
      ier++;
    }
  }
  if (GlobalSEMaxComm(comm, ier) > 0)
    goto ERROR_EXIT;

  
  /* get the basic statistics */
  ier = 0;
  if (mype == 0) {
    sprintf(filename, "%s_stats.txt", fstem);
    fpin = gk_fopen(filename, "r", "DistDGL_ReadGraph: stats.txt");
    if ((i = fscanf(fpin, "%"SCIDX" %"SCIDX" %"SCIDX, &gnvtxs, &gnedges, &ncon)) != 3) {
      printf("ERROR: File '%s' contains %"PRIDX"/3 required information.\n", filename, i);
      ier = 1;
    }
  }
  if (GlobalSEMaxComm(comm, ier) > 0)
    goto ERROR_EXIT;

  ncon++;  /* the +1 is for the node type */
  gkMPI_Bcast(&gnvtxs, 1, IDX_T, 0, comm);
  gkMPI_Bcast(&gnedges, 1, IDX_T, 0, comm);
  gkMPI_Bcast(&ncon, 1, IDX_T, 0, comm);

  printf("[%03"PRIDX"] gnvtxs: %"PRIDX", gnedges: %"PRIDX", ncon: %"PRIDX"\n", 
      mype, gnvtxs, gnedges, ncon);


  /* ======================================================= */
  /* setup the graph structure                               */
  /* ======================================================= */
  graph = CreateGraph();
  graph->gnvtxs = gnvtxs;
  graph->ncon = ncon-1;

  vtxdist = graph->vtxdist = imalloc(npes+1, "DistDGL_ReadGraph: vtxdist");
  for (pe=0; pe<npes; pe++)
    vtxdist[pe] = gnvtxs/npes + (pe < gnvtxs%npes ? 1 : 0);
  MAKECSR(i, npes, vtxdist);
  ASSERT(gnvtxs == vtxdist[npes]);

  nvtxs = graph->nvtxs = vtxdist[mype+1]-vtxdist[mype];

  //printf("[%03"PRIDX"] nvtxs: %"PRIDX", vtxdist: %"PRIDX" %"PRIDX" %"PRIDX" %"PRIDX"\n", 
  //    mype, nvtxs, vtxdist[0], vtxdist[1], vtxdist[2], vtxdist[3]);


  /* ======================================================= */
  /* read and distribute the edges and their metadata */
  /* ======================================================= */
  {
    idx_t u, v, uu, vv, nlinesread, nchunks, chunk, chunksize, lnedges, lnmeta;
    idx_t *coo_buffers_cpos=NULL, *coo_chunks_len=NULL; 
    i2kv_t **coo_buffers=NULL, **coo_chunks=NULL, *lcoo=NULL;
    idx_t *meta_buffers_cpos=NULL, *meta_buffers_len=NULL, *meta_chunks_len=NULL;
    char **meta_buffers=NULL, **meta_chunks=NULL, *lmeta;
    idx_t *emptr;
    char *emdata;
    idx_t firstvtx, lastvtx;

    chunksize = CHUNKSIZE;
    nchunks = 2*gnedges/chunksize;
    nlinesread = 0;
    if (mype == 0) {
      sprintf(filename, "%s_edges.txt", fstem);
      fsize = 2*gk_getfsize(filename)/nchunks;  /* give it a 2x xtra space */
      fpin = gk_fopen(filename, "r", "DistDGL_ReadGraph: edges.txt");
  
      coo_buffers_cpos  = imalloc(npes, "coo_buffers_cpos");
      meta_buffers_cpos = imalloc(npes, "meta_buflen_cpos");
      meta_buffers_len  = ismalloc(npes, fsize, "meta_buflen_len");
  
      coo_buffers  = (i2kv_t **)gk_malloc(npes*sizeof(i2kv_t *), "coo_buffers");
      meta_buffers = (char **)gk_malloc(npes*sizeof(char *), "meta_buffers");
      for (pe=0; pe<npes; pe++) {
        coo_buffers[pe]  = (i2kv_t *)gk_malloc(chunksize*sizeof(i2kv_t), "coo_buffers[pe]");
        meta_buffers[pe] = gk_cmalloc(meta_buffers_len[pe], "meta_buffers[pe]");
      }
    }
  
    /* allocate memory for the chunks that will be collected by each PE */
    coo_chunks_len  = imalloc(nchunks, "coo_chunks_len");
    meta_chunks_len = imalloc(nchunks, "meta_chunks_len");
    coo_chunks  = (i2kv_t **)gk_malloc(nchunks*sizeof(i2kv_t *), "coo_chunks");
    meta_chunks = (char **)gk_malloc(nchunks*sizeof(char *), "meta_chunks");
  
    /*
    if (mype == 0) {
      for (u=0; u<gnvtxs; u++)
        printf("%6d => %6d => %6d [ %2d %5d ]\n", u, 
            DistDGL_mapToCyclic(u, npes, vtxdist), 
            DistDGL_mapFromCyclic(DistDGL_mapToCyclic(u, npes, vtxdist), npes, vtxdist), 
            u%npes, u/npes);
    }
    */

    /* start reading the edge file */
    for (chunk=0;;chunk++) {
      if (mype  == 0) {
        iset(npes, 0, coo_buffers_cpos);
        iset(npes, 0, meta_buffers_cpos);
        nlinesread = 0;
        while (gk_getline(&line, &lnlen, fpin) != -1) {
          rlen = strlen(gk_strtprune(line, "\n\r"));
          nlinesread++;

          /* sscanf(line, "%"SCIDX" %"SCIDX, &uu, &vv); */
          /* NOTE: The edges are read as (dest-id, source-id), as they 
             are assigned to the partitions based on the dest-id; 
             i.e., each partition stores the incoming edges */
          sscanf(line, "%"SCIDX" %"SCIDX, &vv, &uu);
          u = vtxdist[uu%npes] + uu/npes;
          v = vtxdist[vv%npes] + vv/npes;

          ASSERT2(u < gnvtxs);
          ASSERT2(v < gnvtxs);

          /* record the edge in its input direction */
          pe = uu%npes;
          coo_buffers[pe][coo_buffers_cpos[pe]].key1 = u;
          coo_buffers[pe][coo_buffers_cpos[pe]].key2 = v;
          coo_buffers[pe][coo_buffers_cpos[pe]].val  = meta_buffers_cpos[pe];
          coo_buffers_cpos[pe]++;

          /* see if you need to realloc the metadata buffer */
          if (meta_buffers_cpos[pe]+rlen+1 >= meta_buffers_len[pe]) {
            meta_buffers_len[pe] += meta_buffers_len[pe] + rlen + 1;
            meta_buffers[pe] = gk_crealloc(meta_buffers[pe], meta_buffers_len[pe], "meta_buffers[pe]");
          }

          gk_ccopy(rlen+1, line, meta_buffers[pe]+meta_buffers_cpos[pe]); 
          meta_buffers_cpos[pe] += rlen + 1;
  
          /* record the edge in its oppositive direction */
          pe = vv%npes;
          coo_buffers[pe][coo_buffers_cpos[pe]].key1 = v;
          coo_buffers[pe][coo_buffers_cpos[pe]].key2 = u;
          coo_buffers[pe][coo_buffers_cpos[pe]].val  = -1;
          coo_buffers_cpos[pe]++;

          if (coo_buffers_cpos[uu%npes] >= chunksize || coo_buffers_cpos[vv%npes] >= chunksize) 
            break;
        }
      }
  
      /* distributed termination detection */
      if (GlobalSESumComm(comm, nlinesread) == 0)
        break;
  
      /* adjust memory if needed */
      if (chunk >= nchunks) {
        //printf("[%03"PRIDX"] Readjusting nchunks: %"PRIDX"\n", mype, nchunks);
        nchunks *= 2;
        coo_chunks_len  = irealloc(coo_chunks_len, nchunks, "coo_chunks_len");
        meta_chunks_len = irealloc(meta_chunks_len, nchunks, "meta_chunks_len");
        coo_chunks  = (i2kv_t **)gk_realloc(coo_chunks, nchunks*sizeof(i2kv_t *), "coo_chunks");
        meta_chunks = (char **)gk_realloc(meta_chunks, nchunks*sizeof(char *), "meta_chunks");
      }
  
      if (mype == 0) {
        for (pe=1; pe<npes; pe++) {
          gkMPI_Send((void *)&(coo_buffers_cpos[pe]), 1, IDX_T, pe, 0, comm);
          gkMPI_Send((void *)&(meta_buffers_cpos[pe]), 1, IDX_T, pe, 0, comm);
          gkMPI_Send((void *)coo_buffers[pe], coo_buffers_cpos[pe]*sizeof(i2kv_t), MPI_BYTE, pe, 0, comm);
          gkMPI_Send((void *)meta_buffers[pe], meta_buffers_cpos[pe], MPI_CHAR, pe, 0, comm);
        }
        coo_chunks_len[chunk]  = coo_buffers_cpos[0];
        meta_chunks_len[chunk] = meta_buffers_cpos[0];
        coo_chunks[chunk]  = (i2kv_t *)gk_malloc(coo_chunks_len[chunk]*sizeof(i2kv_t), "coo_chunks[chunk]");
        meta_chunks[chunk] = gk_cmalloc(meta_chunks_len[chunk], "meta_chunks[chunk]");
        gk_ccopy(coo_buffers_cpos[0]*sizeof(i2kv_t), (char *)coo_buffers[0], (char *)coo_chunks[chunk]);
        gk_ccopy(meta_buffers_cpos[0], meta_buffers[0], meta_chunks[chunk]);

        //printf("[%03"PRIDX"] chunk: %"PRIDX", u:%"PRIDX", v:%"PRIDX", val:%s\n",
        //    mype, chunk, coo_chunks[chunk][0].key1, coo_chunks[chunk][0].key2, 
        //    (coo_chunks[chunk][0].val == -1 ? "-1" : meta_chunks[chunk]+coo_chunks[chunk][0].val));
      }
      else {
        gkMPI_Recv((void *)&(coo_chunks_len[chunk]), 1, IDX_T, 0, 0, comm, &stat);
        gkMPI_Recv((void *)&(meta_chunks_len[chunk]), 1, IDX_T, 0, 0, comm, &stat);
        coo_chunks[chunk]  = (i2kv_t *)gk_malloc(coo_chunks_len[chunk]*sizeof(i2kv_t), "coo_chunks[chunk]");
        meta_chunks[chunk] = gk_cmalloc(meta_chunks_len[chunk], "meta_chunks[chunk]");
        gkMPI_Recv((void *)coo_chunks[chunk], coo_chunks_len[chunk]*sizeof(i2kv_t), MPI_BYTE, 0, 0, comm, &stat);
        gkMPI_Recv((void *)meta_chunks[chunk], meta_chunks_len[chunk], MPI_CHAR, 0, 0, comm, &stat);

        //printf("[%03"PRIDX"] chunk: %"PRIDX", u:%"PRIDX", v:%"PRIDX", val:%s\n",
        //    mype, chunk, coo_chunks[chunk][0].key1, coo_chunks[chunk][0].key2, 
        //    (coo_chunks[chunk][0].val == -1 ? "-1" : meta_chunks[chunk]+coo_chunks[chunk][0].val));
      }
    }
    nchunks = chunk;

    //printf("[%03"PRIDX"] Final nchunks: %"PRIDX"\n", mype, nchunks);
  
    /* done reading the edge file */
    if (mype == 0) {
      gk_fclose(fpin);
      
      for (pe=0; pe<npes; pe++) 
        gk_free((void **)&coo_buffers[pe], &meta_buffers[pe], LTERM);
  
      gk_free((void **)&coo_buffers_cpos, &meta_buffers_cpos, &meta_buffers_len, 
          &coo_buffers, &meta_buffers, LTERM);
    }
  
    /* consolidate the chunks into lcoo/lmeta lnedges/lnmeta */
    lnedges = isum(nchunks, coo_chunks_len, 1);
    lnmeta  = isum(nchunks, meta_chunks_len, 1);
  
    //printf("[%03"PRIDX"] lnedges: %"PRIDX", lnmeta: %"PRIDX"\n", mype, lnedges, lnmeta);

    lcoo  = (i2kv_t *)gk_malloc(sizeof(i2kv_t)*lnedges, "lcoo");
    lmeta = gk_cmalloc(lnmeta, "lmeta");
  
    lnedges = lnmeta = 0;
    for (chunk=0; chunk<nchunks; chunk++) {
      for (i=0; i<coo_chunks_len[chunk]; i++, lnedges++) {
        lcoo[lnedges] = coo_chunks[chunk][i];
        if (lcoo[lnedges].val != -1)
          lcoo[lnedges].val += lnmeta;
      }
  
      gk_ccopy(meta_chunks_len[chunk], meta_chunks[chunk], lmeta+lnmeta);
      lnmeta += meta_chunks_len[chunk];
  
      gk_free((void **)&coo_chunks[chunk], &meta_chunks[chunk], LTERM);
    }
    gk_free((void **)&coo_chunks_len, &meta_chunks_len, LTERM);
  
    /* sort and remove duplicates */
    i2kvsorti(lnedges, lcoo);
    for (j=0, i=1; i<lnedges; i++) {
      if (lcoo[i].key1 == lcoo[j].key1 && lcoo[i].key2 == lcoo[j].key2) {
        if (lcoo[i].val != -1)
          printf("[%03"PRIDX"]Duplicate edges with metadata: %"PRIDX" [%s][%s]\n",
              mype, i, lmeta+lcoo[i].val, (lcoo[j].val==-1 ? "NULL" : lmeta+lcoo[j].val));
      }
      else {
        lcoo[++j] = lcoo[i];
      }
    }
    lnedges = j+1;

    //printf("[%03"PRIDX"] Done with sorting and de-duplication.\n", mype);
    gkMPI_Barrier(comm);
  
    /* convert the coo into the csr version */
    graph->nvtxs = nvtxs;
    xadj   = graph->xadj   = ismalloc(nvtxs+1, 0, "DistDGL_ReadGraph: xadj");
    adjncy = graph->adjncy = imalloc(lnedges, "DistDGL_ReadGraph: adjncy");
  
    firstvtx = vtxdist[mype];
    lastvtx  = vtxdist[mype+1];
    for (i=0; i<lnedges; i++) {
      ASSERT2(firstvtx <= lcoo[i].key1 && lcoo[i].key1 < lastvtx);
      xadj[lcoo[i].key1-firstvtx]++;
    }
    MAKECSR(i, nvtxs, xadj);
  
    for (i=0; i<lnedges; i++) 
      adjncy[xadj[lcoo[i].key1-firstvtx]++] = lcoo[i].key2;
    SHIFTCSR(i, nvtxs, xadj);
  
    //printf("[%03"PRIDX"] Done with csr conversion.\n", mype);
    gkMPI_Barrier(comm);
  
    /* convert the lmeta into the (emptr, emdata) arrays */
    graph->emdata_size = lnmeta+lnedges*idxwidth;
    emptr  = graph->emptr  = ismalloc(lnedges+1, 0, "DistDGL_ReadGraph: emptr");
    emdata = graph->emdata = gk_cmalloc(graph->emdata_size, "DistDGL_ReadGraph: emdata");
    for (i=0; i<lnedges; i++) {
      if (lcoo[i].val == -1) {
        emptr[i+1] = emptr[i];
      }
      else { 
        j = strlen(lmeta+lcoo[i].val)+1;
        gk_ccopy(j, lmeta+lcoo[i].val, emdata+emptr[i]);
        emptr[i+1] = emptr[i] + ((j+idxwidth-1)/idxwidth)*idxwidth; /* pad them to idxwidth boundaries */
      }
    }

    /* save the emdata into a file for now */
    {
      char fileout[256];
      sprintf(fileout, "emdata-%d-%"PRIDX".bin", (int)getpid(), mype);
      gk_cwritefilebin(fileout, graph->emdata_size, emdata);
      gk_free((void **)&graph->emdata, LTERM);
    }

    gk_free((void **)&lcoo, &lmeta, LTERM);
  }

  //printf("[%03"PRIDX"] Done with edges.\n", mype);
  gkMPI_Barrier(comm);



  /* ======================================================= */
  /* read and distribute the node weights and their metadata */
  /* ======================================================= */
  {
    idx_t u, v, vv, nlinesread, nchunks, chunk, chunksize, lnmeta;
    idx_t *con_buffers_cpos=NULL, **con_buffers=NULL;
    idx_t *con_chunks_len=NULL, **con_chunks=NULL;
    idx_t *meta_buffers_cpos=NULL, *meta_buffers_len=NULL;
    char **meta_buffers=NULL;
    idx_t *meta_chunks_len=NULL;
    char **meta_chunks=NULL;
    char *curstr, *newstr;
    idx_t *vmptr;
    char *vmdata;

    chunksize = CHUNKSIZE;
    nlinesread = 0;
    if (mype == 0) {
      sprintf(filename, "%s_nodes.txt", fstem);
      fsize = 3*gk_getfsize(filename)/(2*npes);  /* give it a 1.5x xtra space */
      fpin = gk_fopen(filename, "r", "DistDGL_ReadGraph: nodes.txt");
  
      con_buffers_cpos  = imalloc(npes, "con_buffers_cpos");
      meta_buffers_cpos = imalloc(npes, "meta_buflen_cpos");
      meta_buffers_len  = ismalloc(npes, fsize, "meta_buflen_len");
  
      con_buffers  = (idx_t **)gk_malloc(npes*sizeof(idx_t *), "con_buffers");
      meta_buffers = (char **)gk_malloc(npes*sizeof(char *), "meta_buffers");
      for (pe=0; pe<npes; pe++) {
        con_buffers[pe]  = imalloc((ncon+1)*chunksize, "con_buffers[pe]");
        meta_buffers[pe] = gk_cmalloc(meta_buffers_len[pe], "meta_buffers[pe]");
      }

      u = 0;
    }
  
    /* allocate memory for the chunks that will be collected by each PE */
    nchunks = 2*gnvtxs/chunksize;
    con_chunks_len  = imalloc(nchunks, "con_chunks_len");
    meta_chunks_len = imalloc(nchunks, "meta_chunks_len");
    con_chunks  = (idx_t **)gk_malloc(nchunks*sizeof(idx_t *), "con_chunks");
    meta_chunks = (char **)gk_malloc(nchunks*sizeof(char *), "meta_chunks");
  
    /* start reading the node file */
    for (chunk=0;;chunk++) {
      if (mype == 0) {
        iset(npes, 0, con_buffers_cpos);
        iset(npes, 0, meta_buffers_cpos);
        nlinesread = 0;
        while (gk_getline(&line, &lnlen, fpin) != -1) {
          rlen = strlen(gk_strtprune(line, "\n\r"));
          nlinesread++;
          pe = u%npes;
          v = vtxdist[u%npes] + u/npes;
          u++;
  
          if (meta_buffers_cpos[pe]+rlen+1 >= meta_buffers_len[pe]) {
            meta_buffers_len[pe] += meta_buffers_len[pe] + rlen + 1;
            meta_buffers[pe] = gk_crealloc(meta_buffers[pe], meta_buffers_len[pe], "meta_buffers[pe]");
          }

          curstr = line;
          newstr = NULL;
          for (i=0; i<ncon; i++) {
            con_buffers[pe][(ncon+1)*con_buffers_cpos[pe]+i] = strtoidx(curstr, &newstr, 10);
            curstr = newstr;
          }
  
          con_buffers[pe][(ncon+1)*con_buffers_cpos[pe]+ncon] = meta_buffers_cpos[pe];
          gk_ccopy(rlen+1, line, meta_buffers[pe]+meta_buffers_cpos[pe]); 
          meta_buffers_cpos[pe] += rlen + 1;
          con_buffers_cpos[pe]++;
  
          if (con_buffers_cpos[pe] >= chunksize) 
            break;
        }
      }
  
      /* distributed termination detection */
      if (GlobalSESumComm(comm, nlinesread) == 0)
        break;
  
      /* adjust memory if needed */
      if (chunk >= nchunks) {
        //printf("[%03"PRIDX"] Readjusting nchunks: %"PRIDX"\n", mype, nchunks);
        nchunks *= 2;
        con_chunks_len  = irealloc(con_chunks_len, nchunks, "con_chunks_len");
        meta_chunks_len = irealloc(meta_chunks_len, nchunks, "meta_chunks_len");
        con_chunks  = (idx_t **)gk_realloc(con_chunks, nchunks*sizeof(idx_t *), "con_chunks");
        meta_chunks = (char **)gk_realloc(meta_chunks, nchunks*sizeof(char *), "meta_chunks");
      }
  
      if (mype == 0) {
        for (pe=1; pe<npes; pe++) {
          gkMPI_Send((void *)&con_buffers_cpos[pe], 1, IDX_T, pe, 0, comm);
          gkMPI_Send((void *)&meta_buffers_cpos[pe], 1, IDX_T, pe, 0, comm);
          gkMPI_Send((void *)con_buffers[pe], (ncon+1)*con_buffers_cpos[pe], IDX_T, pe, 0, comm);
          gkMPI_Send((void *)meta_buffers[pe], meta_buffers_cpos[pe], MPI_CHAR, pe, 0, comm);
        }
        con_chunks_len[chunk]  = con_buffers_cpos[0];
        meta_chunks_len[chunk] = meta_buffers_cpos[0];
        con_chunks[chunk]  = imalloc((ncon+1)*con_chunks_len[chunk], "con_chunks[chunk]");
        meta_chunks[chunk] = gk_cmalloc(meta_chunks_len[chunk], "meta_chunks[chunk]");
        icopy((ncon+1)*con_chunks_len[chunk], con_buffers[0], con_chunks[chunk]);
        gk_ccopy(meta_chunks_len[chunk], meta_buffers[0], meta_chunks[chunk]);
      }
      else {
        gkMPI_Recv((void *)&con_chunks_len[chunk], 1, IDX_T, 0, 0, comm, &stat);
        gkMPI_Recv((void *)&meta_chunks_len[chunk], 1, IDX_T, 0, 0, comm, &stat);
        con_chunks[chunk]  = imalloc((ncon+1)*con_chunks_len[chunk], "con_chunks[chunk]");
        meta_chunks[chunk] = gk_cmalloc(meta_chunks_len[chunk], "meta_chunks[chunk]");
        gkMPI_Recv((void *)con_chunks[chunk], (ncon+1)*con_chunks_len[chunk], IDX_T, 0, 0, comm, &stat);
        gkMPI_Recv((void *)meta_chunks[chunk], meta_chunks_len[chunk], MPI_CHAR, 0, 0, comm, &stat);
      }
    }
    nchunks = chunk;
  
    //printf("[%03"PRIDX"] Final nchunks: %"PRIDX"\n", mype, nchunks);

    /* done reading the node file */
    if (mype == 0) {
      gk_fclose(fpin);
      
      for (pe=0; pe<npes; pe++) 
        gk_free((void **)&con_buffers[pe], &meta_buffers[pe], LTERM);
  
      gk_free((void **)&con_buffers_cpos, &meta_buffers_cpos, &meta_buffers_len, 
          &con_buffers, &meta_buffers, LTERM);
    }
  
    /* populate vwgt and create (vmptr, vmdata) */
    ASSERT2(nvtxs == isum(nchunks, con_chunks_len, 1));
    lnmeta = isum(nchunks, meta_chunks_len, 1);

    //printf("[%03"PRIDX"] nvtxs: %"PRIDX", lnmeta: %"PRIDX"\n", 
    //    mype, isum(nchunks, con_chunks_len, 1), lnmeta);
  
    graph->vmdata_size = lnmeta+nvtxs*idxwidth;
    vmptr  = graph->vmptr  = imalloc(nvtxs+1, "DistDGL_ReadGraph: vmptr");
    vmdata = graph->vmdata = gk_cmalloc(graph->vmdata_size, "DistDGL_ReadGraph: vmdata");

    vwgt  = graph->vwgt  = imalloc(nvtxs*(ncon-1), "DistDGL_ReadGraph: vwgt");
    vtype = graph->vtype = imalloc(nvtxs, "DistDGL_ReadGraph: vwgt");

    nvtxs = 0;
    vmptr[0] = 0;
    for (chunk=0; chunk<nchunks; chunk++) {
      for (i=0; i<con_chunks_len[chunk]; i++, nvtxs++) {
        vtype[nvtxs] = con_chunks[chunk][(ncon+1)*i+0];
        for (j=1; j<ncon; j++) /* the 1st constraint is the vertex type */
          vwgt[(ncon-1)*nvtxs+j-1] = con_chunks[chunk][(ncon+1)*i+j];

        j = strlen(meta_chunks[chunk]+con_chunks[chunk][(ncon+1)*i+ncon])+1;
        gk_ccopy(j+1, meta_chunks[chunk]+con_chunks[chunk][(ncon+1)*i+ncon], vmdata+vmptr[nvtxs]);
        vmptr[nvtxs+1] = vmptr[nvtxs] + ((j+idxwidth-1)/idxwidth)*idxwidth;  /* pad it to idxwidth muptliples */
      }

      gk_free((void **)&con_chunks[chunk], &meta_chunks[chunk], LTERM);
    }
    ASSERT2(nvtxs == vtxdist[mype+1]-vtxdist[mype]);

    /* save the vmdata into a file for now */
    {
      char fileout[256];
      sprintf(fileout, "vmdata-%d-%"PRIDX".bin", (int)getpid(), mype);
      gk_cwritefilebin(fileout, graph->vmdata_size, vmdata);
      gk_free((void **)&graph->vmdata, LTERM);
    }

    gk_free((void **)&con_chunks_len, &meta_chunks_len, LTERM);
  }

#ifdef XXX
  /* write the graph in stdout */
  {
    idx_t u, v, i, j;

    for (pe=0; pe<npes; pe++) {
      if (mype == pe) {
        for (i=0; i<graph->nvtxs; i++) {
          for (j=graph->xadj[i]; j<graph->xadj[i+1]; j++) {
            u = graph->vtxdist[mype]+i;
            v = graph->adjncy[j];

            printf("YYY%"PRIDX" [%"PRIDX" %"PRIDX"] [%"PRIDX" %"PRIDX"] vmdata: [%s]  emdata: [%s]\n", 
                mype,
                u, v, 
                DistDGL_mapFromCyclic(u, npes, vtxdist), 
                DistDGL_mapFromCyclic(v, npes, vtxdist),
                graph->vmdata+graph->vmptr[i],
                graph->emdata+graph->emptr[j]);
          }
        }
        fflush(stdout);
      }
      gkMPI_Barrier(comm);
    }
  }
#endif

ERROR_EXIT:
  gk_free((void **)&filename, &line, LTERM);

  return graph;
}


/*************************************************************************/
/*! Checks the local consistency of moved graph. */
/*************************************************************************/
void DistDGL_WriteGraphs(char *fstem, graph_t *graph, idx_t nparts_per_pe, 
         MPI_Comm comm)
{
  idx_t i, j, ii, jj, k, nvtxs, pnum, firstvtx;
  idx_t *xadj, *adjncy, *vtxdist, *vmptr, *emptr, *where, *vtype, *newlabel;
  char *filename, *vmdata, *emdata;
  FILE *nodefps[nparts_per_pe], *edgefps[nparts_per_pe], *statfps[nparts_per_pe];
  i2kv_t *cand;

  idx_t npes, mype;
  gkMPI_Comm_size(comm, &npes);
  gkMPI_Comm_rank(comm, &mype);

  nvtxs   = graph->nvtxs;
  xadj    = graph->xadj;
  adjncy  = graph->adjncy;
  vtxdist = graph->vtxdist;
  vmptr   = graph->vmptr;
  emptr   = graph->emptr;
  vmdata  = graph->vmdata;
  emdata  = graph->emdata;
  where   = graph->where;
  vtype   = graph->vtype;

  firstvtx = vtxdist[mype*nparts_per_pe];

  {
    ctrl_t *ctrl;
    graph_t *ngraph;
    idx_t *cvtxdist;
    idx_t *nadjncy;

    ctrl = SetupCtrl(PARMETIS_OP_KMETIS, NULL, 1, npes, NULL, NULL, comm); 
    ctrl->CoarsenTo = 1;  /* Needed by SetUpGraph, otherwise we can FP errors */
    cvtxdist = imalloc(npes+1, "cvtxdist");
    for (i=0; i<npes; i++)
      cvtxdist[i] = vtxdist[(i+1)*nparts_per_pe]-vtxdist[i*nparts_per_pe];
    MAKECSR(i, npes, cvtxdist);

    /* allocate memory and copy the adjncy into nadjncy, so that the 
       renumbering will not change the newlabel mapping that you are
       doing bellow */
    nadjncy = imalloc(xadj[nvtxs], "nadjncy");
    icopy(xadj[nvtxs], adjncy, nadjncy);

    ngraph = SetupGraph(ctrl, 1, cvtxdist, xadj, NULL, NULL, nadjncy, NULL, 0);
    AllocateWSpace(ctrl, 0);

    CommSetup(ctrl, ngraph);

    WCOREPUSH;

    cand = (i2kv_t *)gk_malloc(nvtxs*sizeof(i2kv_t), "cand");
    for (i=0; i<nvtxs; i++) {
      cand[i].key1 = where[i];
      cand[i].key2 = vtype[i];
      cand[i].val  = i;
    }
    i2kvsortii(nvtxs, cand);

    newlabel = imalloc(nvtxs+ngraph->nrecv, "newlabel");
    for (i=0; i<nvtxs; i++) 
      newlabel[cand[i].val] = firstvtx+i;

    /* Send the newlabel info to processors storing adjacent interface nodes */
    CommInterfaceData(ctrl, ngraph, newlabel, newlabel+nvtxs);

    WCOREPOP;

    /* copy the nadjncy as this is renumbered */
    icopy(xadj[nvtxs], ngraph->adjncy, adjncy);

    FreeInitialGraphAndRemap(ngraph);
    FreeCtrl(&ctrl);
    gk_free((void **)&cvtxdist, &nadjncy, LTERM);
  }


  filename = gk_malloc(100+strlen(fstem), "DistDGL_WriteGraphs: filename");
  for (k=0; k<nparts_per_pe; k++) {
    sprintf(filename, "p%03"PRIDX"-%s_nodes.txt", mype*nparts_per_pe+k, fstem);
    nodefps[k] = gk_fopen(filename, "w", "DistDGL_ReadGraph: nodes.txt");
    sprintf(filename, "p%03"PRIDX"-%s_edges.txt", mype*nparts_per_pe+k, fstem);
    edgefps[k] = gk_fopen(filename, "w", "DistDGL_ReadGraph: edges.txt");
  }

  for (ii=0; ii<nvtxs; ii++) {
    i = cand[ii].val;
    ASSERT(where[i] >= mype*nparts_per_pe && where[i] < (mype+1)*nparts_per_pe);
    pnum = where[i]%nparts_per_pe;

    fprintf(nodefps[pnum], "%"PRIDX" %s\n", newlabel[i], vmdata+vmptr[i]);

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (emptr[j] < emptr[j+1]) { /* real edge */
        fprintf(edgefps[pnum], "%"PRIDX" %"PRIDX" %s\n", 
            newlabel[i], newlabel[adjncy[j]], emdata+emptr[j]);
      }
    }
  }

  for (k=0; k<nparts_per_pe; k++) {
    gk_fclose(nodefps[k]);
    gk_fclose(edgefps[k]);
  }

  gk_free((void **)&filename, &cand, &newlabel, LTERM);

}


/*************************************************************************/
/*! To and From Cyclic node-ID mapping */
/*************************************************************************/
idx_t DistDGL_mapFromCyclic(idx_t u, idx_t npes, idx_t *vtxdist)
{
  idx_t i;

  for (i=0; i<npes; i++) {
    if (u < vtxdist[i+1])
      break;
  }
  return (u-vtxdist[i])*npes + i; 
}

idx_t DistDGL_mapToCyclic(idx_t u, idx_t npes, idx_t *vtxdist)
{
  return vtxdist[u%npes] + u/npes;
}
  

/*************************************************************************/
/*! Sorts based on increasing <key1, key2>, and decreasing <val> */
/*************************************************************************/
void i2kvsorti(size_t n, i2kv_t *base)
{
#define ikeyval_lt(a, b) \
  ((a)->key1 < (b)->key1 || \
   ((a)->key1 == (b)->key1 && (a)->key2 < (b)->key2) || \
    ((a)->key1 == (b)->key1 && (a)->key2 == (b)->key2 && (a)->val > (b)->val))
  GK_MKQSORT(i2kv_t, base, n, ikeyval_lt);
#undef ikeyval_lt
}

/*************************************************************************/
/*! Sorts based on increasing <key1, key2>, and decreasing <val> */
/*************************************************************************/
void i2kvsortii(size_t n, i2kv_t *base)
{
#define ikeyval_lt(a, b) \
  ((a)->key1 < (b)->key1 || \
   ((a)->key1 == (b)->key1 && (a)->key2 < (b)->key2) || \
    ((a)->key1 == (b)->key1 && (a)->key2 == (b)->key2 && (a)->val < (b)->val))
  GK_MKQSORT(i2kv_t, base, n, ikeyval_lt);
#undef ikeyval_lt
}


