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

#define CHUNKSIZE (1<<16)

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

int DGLPart_GPart(char *fstem, idx_t nparts_per_pe, MPI_Comm comm);
graph_t *DGLPart_ReadGraph(char *fstem, MPI_Comm comm);
void i2kvsorti(size_t n, i2kv_t *base);
void i2kvsortii(size_t n, i2kv_t *base);
void DGLPart_WritePartition(char *fstem, graph_t *graph, idx_t nparts, idx_t *mypart, MPI_Comm comm);
idx_t DGLPart_mapFromCyclic(idx_t u, idx_t npes, idx_t *vtxdist);
idx_t DGLPart_mapToCyclic(idx_t u, idx_t npes, idx_t *vtxdist);


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
      printf("Usage: %s <fstem> <nparts>\n", argv[0]);

    MPI_Finalize();
    exit(0);
  }

  retval = DGLPart_GPart(argv[1], atoi(argv[2]), comm);

  gkMPI_Comm_free(&comm);

  MPI_Finalize();

  return retval;
}


/*************************************************************************/
/*! Partition, move, and save the local graphs */
/*************************************************************************/
int DGLPart_GPart(char *fstem, idx_t nparts, MPI_Comm comm)
{
  idx_t i, npes, mype;
  graph_t *graph, *mgraph;
  idx_t *part;
  idx_t numflag=0, wgtflag=0, options[10], edgecut, ndims;
  real_t *tpwgts=NULL, *ubvec=NULL;

  gkMPI_Comm_size(comm, &npes);
  gkMPI_Comm_rank(comm, &mype);

  /* read and create the graph */
  graph = DGLPart_ReadGraph(fstem, comm);
  gkMPI_Barrier(comm);
  if (graph == NULL)
    return EXIT_FAILURE;

  /* Report peak memory use before partitioning */
  {
    idx_t max = GlobalSEMaxComm(comm, gk_GetProcVmPeak());
    idx_t sum = GlobalSESumComm(comm, gk_GetProcVmPeak());
    if (mype == 0) {
      printf("\nproc/self/stat/VmPeak:\n\tMax: %7.2fMB, Sum: %7.2fMB, Balance: %4.3f\n", 
          (float)max/(1024.0*1024.0), (float)sum/(1024.0*1024.0), 1.0*npes*max/sum);
      fflush(stdout);
    }
    gkMPI_Barrier(comm);
  }

  /*======================================================================
  / Partition the graph 
  /=======================================================================*/
  options[0] = 1;
  options[1] = 15 + (PARMETIS_DBGLVL_TWOHOP|PARMETIS_DBGLVL_FAST|PARMETIS_DBGLVL_DROPEDGES|PARMETIS_DBGLVL_ONDISK);
  //options[1] = 15 + (PARMETIS_DBGLVL_TWOHOP|PARMETIS_DBGLVL_FAST|PARMETIS_DBGLVL_ONDISK);
  options[2] = 1;
  wgtflag = 2;
  numflag = 0;
  edgecut = 0;

  part   = imalloc(graph->nvtxs, "DGLPart_GPart: part");
  tpwgts = rsmalloc(nparts*graph->ncon, 1.0/(real_t)nparts, "DGLPart_GPart: tpwgts");
  ubvec  = rsmalloc(graph->ncon, 1.02, "DGLPart_GPart: unvec");

  if (mype == 0)
    printf("\nDistDGL partitioning, ncon: %"PRIDX", nparts: %"PRIDX" [%s, %s, MPI %d.%d]\n", 
        graph->ncon, nparts, 
        (sizeof(idx_t)==8 ? "i64" : "i32"), (sizeof(real_t)==8 ? "r64" : "r32"),
        MPI_VERSION, MPI_SUBVERSION);

  ParMETIS_V3_PartKway(graph->vtxdist, graph->xadj, graph->adjncy, graph->vwgt, NULL, 
      &wgtflag, &numflag, &(graph->ncon), &nparts, tpwgts, ubvec, options, &edgecut, 
      part, &comm);


  /* Report peak memory use after partitioning */
  {
    idx_t max = GlobalSEMaxComm(comm, gk_GetProcVmPeak());
    idx_t sum = GlobalSESumComm(comm, gk_GetProcVmPeak());
    if (mype == 0) {
      printf("\nproc/self/stat/VmPeak:\n\tMax: %7.2fMB, Sum: %7.2fMB, Balance: %4.3f\n", 
          (float)max/(1024.0*1024.0), (float)sum/(1024.0*1024.0), 1.0*npes*max/sum);
      fflush(stdout);
    }
    gkMPI_Barrier(comm);
  }

  /*======================================================================
  / Write the partitioning vector to disk 
  /=======================================================================*/
  DGLPart_WritePartition(fstem, graph, nparts, part, comm);

  return EXIT_SUCCESS;
}



/*************************************************************************/
/*! Reads, distributes, and pre-processes the DistDGL's input files to 
    create the graph used for partitioning */
/*************************************************************************/
graph_t *DGLPart_ReadGraph(char *fstem, MPI_Comm comm)
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
    filename = gk_malloc(100+strlen(fstem), "DGLPart_ReadGraph: filename");

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
    fpin = gk_fopen(filename, "r", "DGLPart_ReadGraph: stats.txt");
    if ((i = fscanf(fpin, "%"SCIDX" %"SCIDX" %"SCIDX, &gnvtxs, &gnedges, &ncon)) != 3) {
      printf("ERROR: File '%s' contains %"PRIDX"/3 required information.\n", filename, i);
      ier = 1;
    }
    gk_fclose(fpin);
  }
  if (GlobalSEMaxComm(comm, ier) > 0)
    goto ERROR_EXIT;

  ncon++;  /* the +1 is for the node type */
  gkMPI_Bcast(&gnvtxs, 1, IDX_T, 0, comm);
  gkMPI_Bcast(&gnedges, 1, IDX_T, 0, comm);
  gkMPI_Bcast(&ncon, 1, IDX_T, 0, comm);

  if (mype == 0) {
    printf("gnvtxs: %"PRIDX", gnedges: %"PRIDX", ncon: %"PRIDX"\n", gnvtxs, gnedges, ncon-1);
    fflush(stdout);
  }
  gkMPI_Barrier(comm);


  /* ======================================================= */
  /* setup the graph structure                               */
  /* ======================================================= */
  graph = CreateGraph();
  graph->gnvtxs = gnvtxs;
  graph->ncon = ncon-1;

  vtxdist = graph->vtxdist = imalloc(npes+1, "DGLPart_ReadGraph: vtxdist");
  for (pe=0; pe<npes; pe++)
    vtxdist[pe] = gnvtxs/npes + (pe < gnvtxs%npes ? 1 : 0);
  MAKECSR(i, npes, vtxdist);
  ASSERT(gnvtxs == vtxdist[npes]);

  nvtxs = graph->nvtxs = vtxdist[mype+1]-vtxdist[mype];

  //printf("[%03"PRIDX"] nvtxs: %"PRIDX", vtxdist: %"PRIDX" %"PRIDX" %"PRIDX" %"PRIDX"\n", 
  //    mype, nvtxs, vtxdist[0], vtxdist[1], vtxdist[2], vtxdist[3]);


  /* ======================================================= */
  /* read and distribute the edges */
  /* ======================================================= */
  {
    idx_t u, v, uu, vv, nlinesread, nchunks, chunk, chunksize, lnedges;
    idx_t *coo_buffers_cpos=NULL, *coo_chunks_len=NULL; 
    i2kv_t **coo_buffers=NULL, **coo_chunks=NULL, *lcoo=NULL;
    idx_t firstvtx, lastvtx;

    chunksize = CHUNKSIZE;
    nchunks = 1+gnedges/chunksize;
    nlinesread = 0;
    if (mype == 0) {
      sprintf(filename, "%s_edges.txt", fstem);
      fsize = 2*gk_getfsize(filename)/nchunks;  /* give it a 2x xtra space */
      fpin = gk_fopen(filename, "r", "DGLPart_ReadGraph: edges.txt");

      //printf("[edata0]nchunks: %"PRIDX", avg-chunksize: %"PRIDX"\n", nchunks, fsize);
  
      coo_buffers_cpos = imalloc(npes, "coo_buffers_cpos");
      coo_buffers      = (i2kv_t **)gk_malloc(npes*sizeof(i2kv_t *), "coo_buffers");
      for (pe=0; pe<npes; pe++)
        coo_buffers[pe]  = (i2kv_t *)gk_malloc(chunksize*sizeof(i2kv_t), "coo_buffers[pe]");
    }
  
    /* allocate memory for the chunks that will be collected by each PE */
    coo_chunks_len = imalloc(nchunks, "coo_chunks_len");
    coo_chunks     = (i2kv_t **)gk_malloc(nchunks*sizeof(i2kv_t *), "coo_chunks");
  
    /* start reading the edge file */
    for (chunk=0;;chunk++) {
      if (mype  == 0) {
        iset(npes, 0, coo_buffers_cpos);
        nlinesread = 0;
        while (gk_getline(&line, &lnlen, fpin) != -1) {
          rlen = strlen(gk_strtprune(line, "\n\r"));
          nlinesread++;

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
          ASSERT2(coo_buffers_cpos[pe] < chunksize);
          coo_buffers[pe][coo_buffers_cpos[pe]].key1 = u;
          coo_buffers[pe][coo_buffers_cpos[pe]].key2 = v;
          coo_buffers[pe][coo_buffers_cpos[pe]].val  = 1;
          coo_buffers_cpos[pe]++;

          /* record the edge in its oppositive direction */
          pe = vv%npes;
          ASSERT2(coo_buffers_cpos[pe] < chunksize);
          coo_buffers[pe][coo_buffers_cpos[pe]].key1 = v;
          coo_buffers[pe][coo_buffers_cpos[pe]].key2 = u;
          coo_buffers[pe][coo_buffers_cpos[pe]].val  = -1;
          coo_buffers_cpos[pe]++;

          /* the chunksize-1 is to account for the cases in which u and v are
           * assigned to the same pe */
          if (coo_buffers_cpos[uu%npes] >= chunksize-1 || 
              coo_buffers_cpos[vv%npes] >= chunksize-1) 
            break;
        }
      }
  
      /* distributed termination detection */
      if (GlobalSESumComm(comm, nlinesread) == 0)
        break;
  
      /* adjust memory if needed */
      if (chunk >= nchunks) {
        //printf("[%03"PRIDX"] Readjusting nchunks: %"PRIDX"\n", mype, nchunks);
        nchunks *= 1.2;
        coo_chunks_len = irealloc(coo_chunks_len, nchunks, "coo_chunks_len");
        coo_chunks     = (i2kv_t **)gk_realloc(coo_chunks, nchunks*sizeof(i2kv_t *), "coo_chunks");
      }
  
      if (mype == 0) {
        for (pe=1; pe<npes; pe++) {
          gkMPI_Send((void *)&(coo_buffers_cpos[pe]), 1, IDX_T, pe, 0, comm);
          gkMPI_Send((void *)coo_buffers[pe], coo_buffers_cpos[pe]*sizeof(i2kv_t), MPI_BYTE, pe, 0, comm);
        }
        coo_chunks_len[chunk] = coo_buffers_cpos[0];
        coo_chunks[chunk]     = (i2kv_t *)gk_malloc(coo_chunks_len[chunk]*sizeof(i2kv_t), "coo_chunks[chunk]");
        gk_ccopy(coo_buffers_cpos[0]*sizeof(i2kv_t), (char *)coo_buffers[0], (char *)coo_chunks[chunk]);
      }
      else {
        gkMPI_Recv((void *)&(coo_chunks_len[chunk]), 1, IDX_T, 0, 0, comm, &stat);
        coo_chunks[chunk] = (i2kv_t *)gk_malloc(coo_chunks_len[chunk]*sizeof(i2kv_t), "coo_chunks[chunk]");
        gkMPI_Recv((void *)coo_chunks[chunk], coo_chunks_len[chunk]*sizeof(i2kv_t), MPI_BYTE, 0, 0, comm, &stat);
      }
    }
    nchunks = chunk;

    /* done reading the edge file */
    if (mype == 0) {
      //printf("[edata1]nchunks: %"PRIDX", avg-chunksize: %"PRIDX"\n", nchunks, fsize);
      gk_fclose(fpin);
      
      for (pe=0; pe<npes; pe++) 
        gk_free((void **)&coo_buffers[pe], LTERM);
  
      gk_free((void **)&coo_buffers_cpos, &coo_buffers, LTERM);
    }
  
    /* consolidate the chunks into lcoo lnedges */
    lnedges = isum(nchunks, coo_chunks_len, 1);

    lcoo = (i2kv_t *)gk_malloc(sizeof(i2kv_t)*lnedges, "lcoo");
  
    lnedges = 0;
    for (chunk=0; chunk<nchunks; chunk++) {
      for (i=0; i<coo_chunks_len[chunk]; i++, lnedges++) 
        lcoo[lnedges] = coo_chunks[chunk][i];
  
      gk_free((void **)&coo_chunks[chunk], LTERM);
    }
    gk_free((void **)&coo_chunks_len, &coo_chunks, LTERM);

    //printf("[%03"PRIDX"] Done with consolidating the chunks into single arrays.\n", mype);
    gkMPI_Barrier(comm);
  

    /* sort and remove duplicates */
    i2kvsorti(lnedges, lcoo);
    for (j=0, i=1; i<lnedges; i++) {
      if (lcoo[i].key1 == lcoo[j].key1 && lcoo[i].key2 == lcoo[j].key2) {
        if (lcoo[i].val != -1)
          printf("[%03"PRIDX"]Duplicate edges with metadata: %"PRIDX"\n", mype, i);
      }
      else {
        lcoo[++j] = lcoo[i];
      }
    }
    lnedges = j+1;

    //printf("[%03"PRIDX"] Done with sorting and de-duplication.\n", mype);
    gkMPI_Barrier(comm);

    //printf("[%03"PRIDX"] Done with saving emdata into a file.\n", mype);
    gkMPI_Barrier(comm);
  

    /* convert the coo into the csr version */
    graph->nvtxs  = nvtxs;
    graph->nedges = lnedges;
    xadj   = graph->xadj   = ismalloc(nvtxs+1, 0, "DGLPart_ReadGraph: xadj");
    adjncy = graph->adjncy = imalloc(lnedges, "DGLPart_ReadGraph: adjncy");
  
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
  
    gk_free((void **)&lcoo, LTERM);
  }

  //printf("[%03"PRIDX"] Done with edges.\n", mype);
  gkMPI_Barrier(comm);



  /* ======================================================= */
  /* read and distribute the node weights */
  /* ======================================================= */
  {
    idx_t u, nlinesread, nchunks, chunk, chunksize;
    idx_t *con_buffers_cpos=NULL, **con_buffers=NULL;
    idx_t *con_chunks_len=NULL, **con_chunks=NULL;
    char *curstr, *newstr;
    idx_t *vmptr;
    char *vmdata;

    chunksize = CHUNKSIZE;
    nchunks = 1+gnvtxs/chunksize;
    nlinesread = 0;
    if (mype == 0) {
      sprintf(filename, "%s_nodes.txt", fstem);
      fsize = 2*gk_getfsize(filename)/nchunks;  /* give it a 2x xtra space */
      fpin = gk_fopen(filename, "r", "DGLPart_ReadGraph: nodes.txt");
  
      //printf("[vdata0]nchunks: %"PRIDX", avg-chunksize: %"PRIDX"\n", nchunks, fsize);

      con_buffers_cpos = imalloc(npes, "con_buffers_cpos");
      con_buffers      = (idx_t **)gk_malloc(npes*sizeof(idx_t *), "con_buffers");
      for (pe=0; pe<npes; pe++)
        con_buffers[pe]  = imalloc(ncon*chunksize, "con_buffers[pe]");

      u = 0;
    }
  
    /* allocate memory for the chunks that will be collected by each PE */
    con_chunks_len = imalloc(nchunks, "con_chunks_len");
    con_chunks     = (idx_t **)gk_malloc(nchunks*sizeof(idx_t *), "con_chunks");
  
    /* start reading the node file */
    for (chunk=0;;chunk++) {
      if (mype == 0) {
        iset(npes, 0, con_buffers_cpos);
        nlinesread = 0;
        while (gk_getline(&line, &lnlen, fpin) != -1) {
          rlen = strlen(gk_strtprune(line, "\n\r"));
          nlinesread++;
          pe = u%npes;
          u++;
  
          curstr = line;
          newstr = NULL;
          for (i=0; i<ncon; i++) {
            con_buffers[pe][ncon*con_buffers_cpos[pe]+i] = strtoidx(curstr, &newstr, 10);
            curstr = newstr;
          }
  
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
        nchunks *= 1.2;
        con_chunks_len = irealloc(con_chunks_len, nchunks, "con_chunks_len");
        con_chunks     = (idx_t **)gk_realloc(con_chunks, nchunks*sizeof(idx_t *), "con_chunks");
      }
  
      if (mype == 0) {
        for (pe=1; pe<npes; pe++) {
          gkMPI_Send((void *)&con_buffers_cpos[pe], 1, IDX_T, pe, 0, comm);
          gkMPI_Send((void *)con_buffers[pe], ncon*con_buffers_cpos[pe], IDX_T, pe, 0, comm);
        }
        con_chunks_len[chunk] = con_buffers_cpos[0];
        con_chunks[chunk]     = imalloc(ncon*con_chunks_len[chunk], "con_chunks[chunk]");
        icopy(ncon*con_chunks_len[chunk], con_buffers[0], con_chunks[chunk]);
      }
      else {
        gkMPI_Recv((void *)&con_chunks_len[chunk], 1, IDX_T, 0, 0, comm, &stat);
        con_chunks[chunk] = imalloc(ncon*con_chunks_len[chunk], "con_chunks[chunk]");
        gkMPI_Recv((void *)con_chunks[chunk], ncon*con_chunks_len[chunk], IDX_T, 0, 0, comm, &stat);
      }
    }
    nchunks = chunk;
  
    /* done reading the node file */
    if (mype == 0) {
      //printf("[vdata1]nchunks: %"PRIDX", avg-chunksize: %"PRIDX"\n", nchunks, fsize);
      gk_fclose(fpin);
      
      for (pe=0; pe<npes; pe++) 
        gk_free((void **)&con_buffers[pe], LTERM);
  
      gk_free((void **)&con_buffers_cpos, &con_buffers, LTERM);
    }
  
    /* populate vwgt and create (vmptr, vmdata) */
    ASSERT2(nvtxs == isum(nchunks, con_chunks_len, 1));

    vwgt  = graph->vwgt  = imalloc(nvtxs*(ncon-1), "DGLPart_ReadGraph: vwgt");
    vtype = graph->vtype = imalloc(nvtxs, "DGLPart_ReadGraph: vwgt");

    nvtxs = 0;
    for (chunk=0; chunk<nchunks; chunk++) {
      for (i=0; i<con_chunks_len[chunk]; i++, nvtxs++) {
        vtype[nvtxs] = con_chunks[chunk][ncon*i+0];
        for (j=1; j<ncon; j++) /* the 1st constraint is the vertex type */
          vwgt[(ncon-1)*nvtxs+j-1] = con_chunks[chunk][ncon*i+j];
      }

      gk_free((void **)&con_chunks[chunk], LTERM);
    }
    gk_free((void **)&con_chunks_len, LTERM);

    ASSERT2(nvtxs == vtxdist[mype+1]-vtxdist[mype]);
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

            printf("YYY%"PRIDX" [%"PRIDX" %"PRIDX"] [%"PRIDX" %"PRIDX"]\n", 
                mype,
                u, v, 
                DGLPart_mapFromCyclic(u, npes, vtxdist), 
                DGLPart_mapFromCyclic(v, npes, vtxdist),
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
void DGLPart_WritePartition(char *fstem, graph_t *graph, idx_t nparts, 
         idx_t *mypart, MPI_Comm comm)
{
  idx_t i, pe, nvtxs;
  idx_t *vtxdist, *where;
  char *filename;
  FILE *fpout;
  MPI_Status stat;

  idx_t npes, mype;
  gkMPI_Comm_size(comm, &npes);
  gkMPI_Comm_rank(comm, &mype);

  vtxdist = graph->vtxdist;
  nvtxs   = graph->nvtxs;

  if (mype == 0) {
    filename = gk_malloc(100+strlen(fstem), "DGLPart_WritePartition: filename");
    sprintf(filename, "%s_part.%"PRIDX, fstem, nparts);
    fpout = gk_fopen(filename, "w", "DGLPart_WritePartition: fpout");

    for (i=0; i<nvtxs; i++)
      fprintf(fpout, "%"PRIDX" %"PRIDX"\n", DGLPart_mapFromCyclic(i, npes, vtxdist), mypart[i]);

    for (pe=1; pe<npes; pe++) {
      nvtxs = vtxdist[pe+1]-vtxdist[pe];
      where = imalloc(nvtxs, "DGLPart_WritePartition: where");

      gkMPI_Recv((void *)where, nvtxs, IDX_T, pe, 0, comm, &stat);

      for (i=0; i<nvtxs; i++)
        fprintf(fpout, "%"PRIDX" %"PRIDX"\n", 
            DGLPart_mapFromCyclic(vtxdist[pe]+i, npes, vtxdist), where[i]);

      gk_free((void **)&where, LTERM);
    }
    gk_fclose(fpout);
  }
  else {
    gkMPI_Send((void *)mypart, nvtxs, IDX_T, 0, 0, comm);
  }
}



/*************************************************************************/
/*! To and From Cyclic node-ID mapping */
/*************************************************************************/
idx_t DGLPart_mapFromCyclic(idx_t u, idx_t npes, idx_t *vtxdist)
{
  idx_t i;

  for (i=0; i<npes; i++) {
    if (u < vtxdist[i+1])
      break;
  }
  return (u-vtxdist[i])*npes + i; 
}

idx_t DGLPart_mapToCyclic(idx_t u, idx_t npes, idx_t *vtxdist)
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


