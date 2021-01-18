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

#define CHUNKSIZE (1<<14)

int DistDGL_GPart(char *fstem, idx_t nparts_per_pe, MPI_Comm comm);
graph_t *DistDGL_ReadGraph(char *fstem, MPI_Comm comm);
graph_t *DistDGL_MoveGraph(graph_t *ograph, idx_t *part, MPI_Comm comm);
void i2kvsorti(size_t n, i2kv_t *base);

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
  idx_t *part, *cpart, *mpart;
  idx_t numflag=0, wgtflag=0, options[10], edgecut, ndims;
  real_t *tpwgts=NULL, *ubvec=NULL;

  gkMPI_Comm_size(comm, &npes);
  gkMPI_Comm_rank(comm, &mype);

  /* read and create the graph */
  graph = DistDGL_ReadGraph(fstem, comm);
  gkMPI_Barrier(comm);
  if (graph == NULL)
    return EXIT_FAILURE;

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

  return EXIT_SUCCESS;
}


/*************************************************************************/
/*! Reads, distributes, and pre-processes the DistDGL's input files to 
    create the graph used for partitioning */
/*************************************************************************/
graph_t *DistDGL_ReadGraph(char *fstem, MPI_Comm comm)
{
  idx_t i, j, pe;
  idx_t npes, mype, ier;
  graph_t *graph=NULL;
  idx_t gnvtxs, gnedges, nvtxs, lnvtxs, ncon; 
  idx_t *vtxdist, *xadj, *adjncy, *vwgt;
  MPI_Status stat;
  ssize_t rlen, fsize;
  size_t lnlen=0;
  char *filename=NULL, *line=NULL;
  FILE *fpin=NULL;


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
  graph->ncon = ncon;

  vtxdist = graph->vtxdist = ismalloc(npes+1, 0, "DistDGL_ReadGraph: vtxdist");
  lnvtxs = (gnvtxs+npes-1)/npes;
  for (pe=0; pe<npes; pe++)
    vtxdist[pe] = pe*lnvtxs;
  vtxdist[npes] = gnvtxs;

  nvtxs = graph->nvtxs = vtxdist[mype+1]-vtxdist[mype];

  //printf("[%03"PRIDX"] nvtxs: %"PRIDX", lnedges: %"PRIDX"\n", mype, nvtxs, lnvtxs);


  /* ======================================================= */
  /* read and distribute the edges and their metadata */
  /* ======================================================= */
  {
    idx_t u, v, nlinesread, nchunks, chunk, chunksize, lnedges, lnmeta;
    idx_t *coo_buffers_cpos=NULL, *coo_chunks_len=NULL; 
    i2kv_t **coo_buffers=NULL, **coo_chunks=NULL, *lcoo=NULL;
    idx_t *meta_buffers_cpos=NULL, *meta_buffers_len=NULL, *meta_chunks_len=NULL;
    char **meta_buffers=NULL, **meta_chunks=NULL, *lmeta;
    idx_t *emptr;
    char *emdata;
    idx_t firstvtx, lastvtx;

    chunksize = CHUNKSIZE;
    nlinesread = 0;
    if (mype == 0) {
      sprintf(filename, "%s_edges.txt", fstem);
      fsize = 3*gk_getfsize(filename)/(2*npes);  /* give it a 1.5x xtra space */
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
    nchunks = 2*gnedges/chunksize;
    coo_chunks_len  = imalloc(nchunks, "coo_chunks_len");
    meta_chunks_len = imalloc(nchunks, "meta_chunks_len");
    coo_chunks  = (i2kv_t **)gk_malloc(nchunks*sizeof(i2kv_t *), "coo_chunks");
    meta_chunks = (char **)gk_malloc(nchunks*sizeof(char *), "meta_chunks");
  
    /* start reading the edge file */
    for (chunk=0;;chunk++) {
      if (mype  == 0) {
        iset(npes, 0, coo_buffers_cpos);
        iset(npes, 0, meta_buffers_cpos);
        nlinesread = 0;
        while (gk_getline(&line, &lnlen, fpin) != -1) {
          rlen = strlen(gk_strtprune(line, "\n\r"));
          nlinesread++;
          sscanf(line, "%"SCIDX" %"SCIDX, &u, &v);
          u = vtxdist[u%npes] + u/npes;
          v = vtxdist[v%npes] + v/npes;
  
          //if (nlinesread == 1)
          //  printf("u: %"PRIDX" v:%"PRIDX" mdata: %s\n", u, v, line);

          /* see if you need to realloc the metadata buffer */
          if (meta_buffers_cpos[pe]+rlen+1 >= meta_buffers_len[pe]) {
            meta_buffers_len[pe] += meta_buffers_len[pe] + rlen + 1;
            meta_buffers[pe] = gk_crealloc(meta_buffers[pe], meta_buffers_len[pe], "meta_buffers[pe]");
          }

          /* record the edge in its input direction */
          pe = u/lnvtxs;
          coo_buffers[pe][coo_buffers_cpos[pe]].key1 = u;
          coo_buffers[pe][coo_buffers_cpos[pe]].key2 = v;
          coo_buffers[pe][coo_buffers_cpos[pe]].val  = meta_buffers_cpos[pe];
          coo_buffers_cpos[pe]++;

          gk_ccopy(rlen+1, line, meta_buffers[pe]+meta_buffers_cpos[pe]); 
          //if (nlinesread <= 2)
          //  printf("u: %"PRIDX" v:%"PRIDX" mdata: %s %s\n", u, v, line, meta_buffers[pe]+meta_buffers_cpos[pe]);
          meta_buffers_cpos[pe] += rlen + 1;
  
          /* record the edge in its oppositive direction */
          pe = v/lnvtxs;
          coo_buffers[pe][coo_buffers_cpos[pe]].key1 = v;
          coo_buffers[pe][coo_buffers_cpos[pe]].key2 = u;
          coo_buffers[pe][coo_buffers_cpos[pe]].val  = -1;
          coo_buffers_cpos[pe]++;

          if (coo_buffers_cpos[u/lnvtxs] >= chunksize || coo_buffers_cpos[v/lnvtxs] >= chunksize) 
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
    lnedges = j;
  
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
    emptr  = graph->emptr  = ismalloc(lnedges+1, 0, "DistDGL_ReadGraph: emptr");
    emdata = graph->emdata = gk_cmalloc(lnmeta, "DistDGL_ReadGraph: emdata");
    for (i=0; i<lnedges; i++) {
      if (lcoo[i].val == -1) {
        emptr[i+1] = emptr[i];
      }
      else { 
        j = strlen(lmeta+lcoo[i].val)+1;
        gk_ccopy(j, lmeta+lcoo[i].val, emdata+emptr[i]);
        emptr[i+1] = emptr[i]+j;
      }
    }
  
    gk_free((void **)&lcoo, &lmeta, LTERM);
  }

  //printf("[%03"PRIDX"] Done with edges.\n", mype);
  gkMPI_Barrier(comm);



  /* ======================================================= */
  /* read and distribute the node weights and their metadata */
  /* ======================================================= */
  {
    idx_t u, v, nlinesread, nchunks, chunk, chunksize, lnmeta;
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
  
    /* start reading the edge file */
    for (chunk=0;;chunk++) {
      if (mype == 0) {
        iset(npes, 0, con_buffers_cpos);
        iset(npes, 0, meta_buffers_cpos);
        nlinesread = 0;
        while (gk_getline(&line, &lnlen, fpin) != -1) {
          rlen = strlen(gk_strtprune(line, "\n\r"));
          nlinesread++;
          v = vtxdist[u%npes] + u/npes;
          u++;
          pe = v/lnvtxs;
  
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
  
    vmptr  = graph->vmptr  = ismalloc(nvtxs+1, 0, "DistDGL_ReadGraph: vmptr");
    vmdata = graph->vmdata = gk_cmalloc(lnmeta, "DistDGL_ReadGraph: vmdata");

    vwgt = graph->vwgt = imalloc(nvtxs*ncon, "DistDGL_ReadGraph: vwgt");

    nvtxs = lnmeta = 0;
    for (chunk=0; chunk<nchunks; chunk++) {
      for (i=0; i<con_chunks_len[chunk]; i++, nvtxs++) {
        for (j=0; j<ncon; j++)
          vwgt[ncon*nvtxs+j] = con_chunks[chunk][(ncon+1)*i+j];
        vmptr[nvtxs] = strlen(meta_chunks[chunk]+con_chunks[chunk][(ncon+1)*i+ncon])+1;
      }

      gk_ccopy(meta_chunks_len[chunk], meta_chunks[chunk], vmdata+lnmeta);
      lnmeta += meta_chunks_len[chunk];
  
      gk_free((void **)&con_chunks[chunk], &meta_chunks[chunk], LTERM);
    }
    ASSERT2(nvtxs == vtxdist[mype+1]-vtxdist[mype]);
    MAKECSR(i, nvtxs, vmptr);

    gk_free((void **)&con_chunks_len, &meta_chunks_len, LTERM);
  }


ERROR_EXIT:
  gk_free((void **)&filename, &line, LTERM);

  return graph;
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
  graph = SetupGraph(ctrl, 1, ograph->vtxdist, ograph->xadj, ograph->vwgt, ograph->vsize,
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
