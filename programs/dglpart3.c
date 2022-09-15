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

/*************************************************************************
* The following data structure stores key-key tuple 
**************************************************************************/
typedef struct edge_t {
  idx_t u, v;
} edge_t;

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

int DGLPart_GPart(char *fstem, idx_t nparts_per_pe, char *lstnfiles, char *lstefiles, MPI_Comm comm);
graph_t *DGLPart_ReadGraph(char *fstem, char *lstnfiles, char *lstefiles, MPI_Comm comm);
void edgesorti(size_t n, edge_t *base);
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

  if (argc != 5) {
    if (mype == 0)
      printf("Usage: %s <fstem> <nparts> <nodes_file> <edges_file>\n", argv[0]);

    MPI_Finalize();
    exit(0);
  }

  retval = DGLPart_GPart(argv[1], atoi(argv[2]), argv[3], argv[4], comm);

  gkMPI_Comm_free(&comm);

  MPI_Finalize();

  return retval;
}


/*************************************************************************/
/*! Partition, move, and save the local graphs */
/*************************************************************************/
int DGLPart_GPart(char *fstem, idx_t nparts, char *lstnfiles, char *lstefiles, MPI_Comm comm)
{
  idx_t i, npes, mype;
  graph_t *graph, *mgraph;
  idx_t *part;
  idx_t numflag=0, wgtflag=0, options[10], edgecut, ndims;
  real_t *tpwgts=NULL, *ubvec=NULL;

  gkMPI_Comm_size(comm, &npes);
  gkMPI_Comm_rank(comm, &mype);

  /* read and create the graph */
  graph = DGLPart_ReadGraph(fstem, lstnfiles, lstefiles, comm);
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
graph_t *DGLPart_ReadGraph(char *fstem, char *lstnfiles, char *lstefiles, MPI_Comm comm)
{
  idx_t i, j, pe, idxwidth;
  idx_t npes, mype, ier;
  graph_t *graph=NULL;
  idx_t gnvtxs, gnedges, nvtxs, ncon; 
  idx_t *vtxdist, *xadj, *adjncy, *vwgt, *vtype;
  MPI_Status stat;
  size_t lnlen=0;
  char *filename=NULL, *line=NULL;
  FILE *fpin=NULL;
  FILE *fpinaux=NULL;

  idxwidth = sizeof(idx_t);

  gkMPI_Comm_size(comm, &npes);
  gkMPI_Comm_rank(comm, &mype);

  /* check to see if the required files exist */ 
  ier = 0;
  filename = gk_malloc(100+strlen(fstem), "DGLPart_ReadGraph: filename");

  sprintf(filename, "%s_stats.txt", fstem);
  if (!gk_fexists(filename)) {
    printf("ERROR: File '%s' does not exists.\n", filename);
    ier++;
  }

  //sprintf(filename, "%s_edges_%02d.txt", fstem, mype);
  if (!gk_fexists(lstnfiles)) {
    printf("ERROR: File '%s' does not exists.\n", filename);
    ier++;
  }

  //sprintf(filename, "%s_nodes_%02d.txt", fstem, mype);
  if (!gk_fexists(lstefiles)) {
    printf("ERROR: File '%s' does not exists.\n", filename);
    ier++;
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
    idx_t edgecount;
    idx_t *coo_buffers_cpos=NULL, *coo_chunks_len=NULL;
    edge_t **coo_buffers=NULL, **coo_chunks=NULL, *lcoo=NULL;
    edge_t *snd_buffer=NULL, *rcv_buffer=NULL;
    idx_t firstvtx, lastvtx;
    idx_t mflinesread;
    idx_t *rcv_buffer_count=NULL;
    idx_t *snd_displs=NULL, *rcv_displs=NULL;
    char *newstr, *curstr;

    chunksize = CHUNKSIZE;
    nchunks = 1+gnedges/chunksize;
    edgecount = 1+gnedges/npes;
    nlinesread = 0;
    mflinesread = 0;
    printf("Reading edge file: %s\n", lstefiles);
    fpin = gk_fopen(lstefiles, "r", "DGLPart_ReadGraph: lstefile.txt");

    coo_buffers_cpos = imalloc(npes, "coo_buffers_cpos");
    rcv_buffer_count = imalloc(npes, "rcv_buffer_count");
    snd_displs = imalloc(npes, "snd_displacements");
    rcv_displs = imalloc(npes, "rcv_displacements");

    coo_buffers = (edge_t **)gk_malloc(npes*sizeof(edge_t *), "coo_buffers");
    for (pe=0; pe<npes; pe++){
      coo_buffers[pe]  = (edge_t *)gk_malloc(chunksize*sizeof(edge_t), "coo_buffers[pe]");
    }
    snd_buffer  = (edge_t *)gk_malloc(npes*chunksize*sizeof(edge_t), "snd_buffer");
    rcv_buffer  = (edge_t *)gk_malloc(npes*chunksize*sizeof(edge_t), "rcv_buffer");

    /* allocate memory for the chunks that will be collected by each PE */
    coo_chunks_len = imalloc(nchunks, "coo_chunks_len");
    coo_chunks     = (edge_t **)gk_malloc(nchunks*sizeof(edge_t *), "coo_chunks");

    /* allocate memory for storing the edges */
    lcoo = (edge_t *)gk_malloc(sizeof(edge_t)*edgecount, "lcoo");
    printf("Allocated mem. for edges: %d\n", edgecount);
  
    /* start reading the edge file */
    lnedges = 0;
    chunk = 0;
    mflinesread = -1;

    while(gk_getline(&line, &lnlen, fpin) != -1) {
      mflinesread ++;
      if((mflinesread % npes) != mype)
        continue;

      /* Token this line for file_name and global starting nid for this files nodes */
      curstr = line;
      newstr = NULL;

      idx_t ln = strlen(curstr) - 1; //newline character.
      if (*curstr && curstr[ln] == '\n')
        curstr[ln] = 0;
      fpinaux = gk_fopen(curstr, "r", "Reading files for edges");

      for (;;chunk++) {
        iset(npes, 0, coo_buffers_cpos);

        nlinesread = 0;
        while (gk_getline(&line, &lnlen, fpinaux) != -1) {
          nlinesread++;

          /* NOTE: The edges are read as (dest-id, source-id), as they 
            are assigned to the partitions based on the dest-id; 
            i.e., each partition stores the incoming edges */
          sscanf(line, "%"SCIDX" %"SCIDX, &vv, &uu);
          u = vtxdist[uu%npes] + uu/npes;
          v = vtxdist[vv%npes] + vv/npes;

          GKASSERT(u < gnvtxs);
          GKASSERT(v < gnvtxs);

          /* record the edge in its input direction */
          pe = uu%npes;
          ASSERT2(coo_buffers_cpos[pe] < chunksize);
          coo_buffers[pe][coo_buffers_cpos[pe]].u = u;
          coo_buffers[pe][coo_buffers_cpos[pe]].v = v;
          coo_buffers_cpos[pe]++;

          /* record the edge in its oppositive direction */
          pe = vv%npes;
          GKASSERT(coo_buffers_cpos[pe] < chunksize);
          coo_buffers[pe][coo_buffers_cpos[pe]].u = v;
          coo_buffers[pe][coo_buffers_cpos[pe]].v = u;
          coo_buffers_cpos[pe]++;

          /* the chunksize-1 is to account for the cases in which u and v are
          * assigned to the same pe */
          if (coo_buffers_cpos[uu%npes] >= chunksize-1 || 
              coo_buffers_cpos[vv%npes] >= chunksize-1) 
            break;

        } //end of while loop for reading individual edges file
  
        /* distributed termination detection */
        if (GlobalSESumComm(comm, nlinesread) == 0)
          break;

	      /* adjust memory if needed */
        if (chunk >= nchunks) {
          nchunks *= 1.2;
          coo_chunks_len = irealloc(coo_chunks_len, nchunks, "coo_chunks_len");
          coo_chunks     = (edge_t **)gk_realloc(coo_chunks, nchunks*sizeof(edge_t *), "coo_chunks");
        }

        //send/recv counts
        iset(npes, 0, rcv_buffer_count);
        gkMPI_Alltoall((void *)coo_buffers_cpos, 1, IDX_T, (void *)rcv_buffer_count, 1, IDX_T, comm);

        /* prepare message to send to others in the comm. world. */
        iset(npes, 0, rcv_displs);
        rcv_displs[0] = 0;
        for(pe=1; pe<npes; pe++)
          rcv_displs[pe] = rcv_displs[pe-1] + rcv_buffer_count[pe-1];

        i = 0; 
        for(pe=0; pe<npes; pe++){
          gk_ccopy(coo_buffers_cpos[pe]*sizeof(edge_t), (char *)coo_buffers[pe], ((char *)snd_buffer)+i*sizeof(edge_t));
          snd_displs[pe] = i*sizeof(edge_t);
          i += coo_buffers_cpos[pe];

          rcv_displs[pe] *= sizeof(edge_t);
          rcv_buffer_count[pe] *= sizeof(edge_t);
          coo_buffers_cpos[pe] *= sizeof(edge_t); 
        }

        gkMPI_Alltoallv(snd_buffer, coo_buffers_cpos, snd_displs, MPI_BYTE, 
	                rcv_buffer, rcv_buffer_count, rcv_displs, MPI_BYTE, comm);

        /* store the received edges in chunks */
        coo_chunks_len[chunk] = isum(npes, rcv_buffer_count, 1) / sizeof(edge_t);
        coo_chunks[chunk] = (edge_t *)gk_malloc(coo_chunks_len[chunk]*sizeof(edge_t), "coo_chunks[chunk]");
        gk_ccopy(coo_chunks_len[chunk]*sizeof(edge_t), (char *)rcv_buffer, (char *)coo_chunks[chunk]);

        for(i=0; i<coo_chunks_len[chunk]; i++){
          edge_t *temp = coo_chunks[chunk] + i;
          GKASSERT (temp->u < gnvtxs); GKASSERT( temp->u >= 0);
          GKASSERT (temp->v < gnvtxs); GKASSERT( temp->v >= 0);
        }
				
        if(mype == 0){
          printf("[Rank: %d] ChunkID: %d has edges: %d\n", mype, chunk, coo_chunks_len[chunk]);
        }
      } // end of chunks for loop

      /*
      * Note that each process has equal no. of nodes/edges files to read. 
      * And all the nodes/edges files will have same no. of lines, except the very last file which might have fewer lines
      * So this termination condition should work for all processes. 
      */
    }//end of main-while loop for the meta-edges file
    nchunks = chunk;

    /* done reading the edge file */
    gk_fclose(fpin);
    gk_fclose(fpinaux);
      
    for (pe=0; pe<npes; pe++) 
      gk_free((void **)&coo_buffers[pe], LTERM);
    gk_free((void **)&coo_buffers_cpos, &coo_buffers, LTERM);

    gk_free((void **)&snd_buffer, LTERM);
    gk_free((void **)&rcv_buffer, LTERM);
    gk_free((void **)&rcv_buffer_count, LTERM);

    lnedges = isum(nchunks, coo_chunks_len, 1);
    lcoo = (edge_t *)gk_malloc(sizeof(edge_t)*lnedges, "lcoo");
    lnedges = 0;
    for(chunk=0; chunk<nchunks; chunk++){
      for(i=0; i<coo_chunks_len[chunk]; i++, lnedges++)
        lcoo[lnedges] = coo_chunks[chunk][i];
      gk_free((void **)&coo_chunks[chunk], LTERM);
    }
    gk_free((void **)&coo_chunks_len, &coo_chunks, LTERM);
    printf("[Rank: %d] Copied edges to coo buffers: %d\n", mype, lnedges);

    /* sort and remove duplicates */
    edgesorti(lnedges, lcoo);
    for (j=0, i=1; i<lnedges; i++) {
      if (lcoo[i].u != lcoo[j].u || lcoo[i].v != lcoo[j].v) 
        lcoo[++j] = lcoo[i];
    }
    lnedges = j+1;
    printf("[Rank: %d] After sorting and removing duplicates: %d\n", mype, lnedges);

    /* convert the coo into the csr version */
    graph->nvtxs  = nvtxs;
    graph->nedges = lnedges;
    xadj   = graph->xadj   = ismalloc(nvtxs+1, 0, "DGLPart_ReadGraph: xadj");
    adjncy = graph->adjncy = imalloc(lnedges, "DGLPart_ReadGraph: adjncy");
    printf("[Rank: %d] Graph construction: Vertices: %d and edges: %d\n", mype, nvtxs+1, lnedges);
  
    firstvtx = vtxdist[mype];
    lastvtx  = vtxdist[mype+1];
    for (i=0; i<lnedges; i++) {
      GKASSERT(firstvtx <= lcoo[i].u && lcoo[i].u < lastvtx);
      xadj[lcoo[i].u-firstvtx]++;
    }
    printf("[Rank: %d] Completed coo entries\n", mype);
    MAKECSR(i, nvtxs, xadj);
    printf("[Rank: %d] Completed CSR conversion\n", mype);
  
    for (i=0; i<lnedges; i++) 
      adjncy[xadj[lcoo[i].u-firstvtx]++] = lcoo[i].v;
    SHIFTCSR(i, nvtxs, xadj);
  
    gk_free((void **)&lcoo, LTERM);
    printf("[Rank: %d] Done processing edges: %d\n", mype, lnedges);
  }


  /* ======================================================= */
  /* read and distribute the node weights */
  /* ======================================================= */
  {
    idx_t u, nlinesread, nchunks, chunk, chunksize;
    idx_t *con_buffers_cpos=NULL, **con_buffers=NULL;
    idx_t *con_chunks_len=NULL, **con_chunks=NULL, **con_id_chunks=NULL;
    char *curstr, *newstr;
    idx_t *vmptr;
    char *vmdata;
    idx_t lvtxs;
    idx_t *s_vwgt, *s_vtype;
    idx_t start_nid, end_nid;

    idx_t gnid;
    idx_t mlinesread;
    idx_t *snd_buffer_count=NULL, *snd_buffer=NULL;
    idx_t *rcv_buffer_count=NULL, *rcv_buffer=NULL;
    idx_t *snd_displs=NULL, *rcv_displs=NULL;
    idx_t **con_id_buffers=NULL;
    edge_t *idx_buffer=NULL;

    chunksize = CHUNKSIZE;
    nchunks = 1+gnvtxs/chunksize;
    lvtxs = 1+gnvtxs/npes;
    nlinesread = 0;
    mlinesread = 0;

    /*open the metadata file which has a list of node file names followed by starting global nid 
    * for the nodes present in any given file */
    printf("Opening the meta nodes file: %s\n", lstnfiles);
    fpin = gk_fopen(lstnfiles, "r", "DGLPart_ReadGraph: meta file for nodes");
  
    con_buffers_cpos = imalloc(npes, "con_buffers_cpos");
    con_buffers      = (idx_t **)gk_malloc(npes*sizeof(idx_t *), "con_buffers");

    /* to store global nids */
    con_id_buffers	 = (idx_t **)gk_malloc(npes*sizeof(idx_t *), "con_id_buffers");

    /* buffers to send/recv counts */
    snd_buffer_count = imalloc(npes, "snd_buffers_count");
    rcv_buffer_count = imalloc(npes, "rcv_buffers_count");
    snd_displs = imalloc(npes, "Snd_displs");
    rcv_displs = imalloc(npes, "Rcv_displs");

    /* buffers to send/recv data */
    snd_buffer			 = imalloc(npes*ncon*chunksize, "snd_buffer");
    rcv_buffer			 = imalloc(npes*ncon*chunksize, "rcv_buffer");

    for (pe=0; pe<npes; pe++){
      con_buffers[pe]  = imalloc(ncon*chunksize, "con_buffers[pe]");
      /*allocate buffers to store nids for each chunk*/
      /* this is later copied into con_id_chunks for later processing */
      con_id_buffers[pe] = imalloc(chunksize, "con_id_buffers[pe]");
    }

    /* allocate memory for the chunks that will be collected by each PE */
    con_chunks_len = imalloc(nchunks, "con_chunks_len");
    con_chunks     = (idx_t **)gk_malloc(nchunks*sizeof(idx_t *), "con_chunks");    
    con_id_chunks  = (idx_t **)gk_malloc(nchunks*sizeof(idx_t *), "con_id_chunks");

    u = 0;
    gnid = 0;
    lvtxs = 0;
    start_nid = end_nid = 0;
  
    /* start reading the node file */
    chunk = 0;
    mlinesread = -1;
    printf("Starting the main loop for edges: \n");
    while(gk_getline(&line, &lnlen, fpin) != -1){

      mlinesread ++;
      if((mlinesread % npes) != mype)
        continue;

      /* Token this line for file_name and global starting nid for this files nodes */
      curstr = line;
      curstr = strtok(curstr, " ");
      start_nid = atoi(strtok(NULL, " "));
      end_nid = atoi(strtok(NULL, " "));

      idx_t ln = strlen(curstr) - 1; //newline character.
      if (*curstr && curstr[ln] == '\n')
        curstr[ln] = 0;
      fpinaux = gk_fopen(curstr, "r", "Reading files for edges");

      gnid = start_nid; //should initialize this here... 

      for (;;chunk++) {
        iset(npes, 0, con_buffers_cpos);
      	nlinesread = 0;

      	while (gk_getline(&line, &lnlen, fpinaux) != -1) {
          nlinesread++;
          //pe = u%npes;
	  pe = gnid%npes;
  
          curstr = line;
          newstr = NULL;
          for (i=0; i<ncon; i++) {
            con_buffers[pe][ncon*con_buffers_cpos[pe]+i] = strtoidx(curstr, &newstr, 10);
	    if(i == 0){
	      if ((con_buffers[pe][ncon*con_buffers_cpos[pe]+i] >= 4) || (con_buffers[pe][ncon*con_buffers_cpos[pe]+i] < 0)){
	        printf("[Rank: %d] For vertex: %d node type read as : %d\n", mype, gnid, con_buffers[pe][ncon*con_buffers_cpos[pe]+i]);
	        GKASSERT( con_buffers[pe][ncon*con_buffers_cpos[pe]+i] >= 0);
	        GKASSERT( con_buffers[pe][ncon*con_buffers_cpos[pe]+i] < 4);
	      }
	    }
            curstr = newstr;
          }
 
          con_id_buffers[pe][con_buffers_cpos[pe]] = gnid; //u;
          con_buffers_cpos[pe]++;

	  gnid++;
          if (con_buffers_cpos[pe] >= chunksize) 
            break;
      	}
      	printf("[Rank: %d] Done reading nodes chnk: %d\n", mype, chunk);
  
      	/* distributed termination detection */
      	if (GlobalSESumComm(comm, nlinesread) == 0)
          break;

      	/* adjust memory if needed */
      	if (chunk >= nchunks) {
          nchunks *= 1.2;
          con_chunks_len = irealloc(con_chunks_len, nchunks, "con_chunks_len");
          con_chunks = (idx_t **)gk_realloc(con_chunks, nchunks*sizeof(idx_t *), "con_chunks");
	  con_id_chunks	= (idx_t **)gk_realloc(con_id_chunks, nchunks*sizeof(idx_t *), "con_id_chunks");
      	}

	/* prepare sending global nids to respective ranks here */
        iset(npes, 0, rcv_buffer_count);
        gkMPI_Alltoall((void *)con_buffers_cpos, 1, IDX_T, (void *)rcv_buffer_count, 1, IDX_T, comm);

        /* prepare send message */
        u = 0;
        for(pe=0; pe<npes; pe++){
          icopy(con_buffers_cpos[pe], (void *)con_id_buffers[pe], (void *)&snd_buffer[u]);
          u += con_buffers_cpos[pe];
        }

        snd_displs[0] = 0;
        rcv_displs[0] = 0;
        for(pe=1; pe<npes; pe++){
          snd_displs[pe] = snd_displs[pe-1] + con_buffers_cpos[pe-1];
          rcv_displs[pe] = rcv_displs[pe-1] + rcv_buffer_count[pe-1];
        }

        /* send out the global node_ids to all processes, and store the owned ones */
        gkMPI_Alltoallv((void *)snd_buffer, con_buffers_cpos, snd_displs, IDX_T, (void *)rcv_buffer, rcv_buffer_count, rcv_displs, IDX_T, comm);
        con_chunks_len[chunk] = isum(npes, rcv_buffer_count, 1);
        con_id_chunks[chunk] = imalloc(con_chunks_len[chunk], "con_chunks[chunk]");
        icopy(con_chunks_len[chunk], rcv_buffer, con_id_chunks[chunk]); /* con_chunk_nids has global-nids after cyclic exchange */
        printf("[Rank: %d] copied global_nids: %d \n", mype, con_chunks_len[chunk]);

        /* send out the node weights now to other ranks, as before */
        /* prepare snd_displs, rcv_displs, and message here */
        u = 0;
        for(pe=0; pe<npes; pe++){
          con_buffers_cpos[pe] = ncon*con_buffers_cpos[pe];
	  rcv_buffer_count[pe] = ncon*rcv_buffer_count[pe];

	  icopy(con_buffers_cpos[pe], (void *)con_buffers[pe], (void *)&snd_buffer[u]);	
	  u += con_buffers_cpos[pe];
	}

	snd_displs[0] = 0;
	rcv_displs[0] = 0;
	for(pe=1; pe<npes; pe++){
	  snd_displs[pe] = snd_displs[pe-1] + con_buffers_cpos[pe-1];
	  rcv_displs[pe] = rcv_displs[pe-1] + rcv_buffer_count[pe-1];
	}

        gkMPI_Alltoallv((void *)snd_buffer, con_buffers_cpos, snd_displs, IDX_T, (void *)rcv_buffer, rcv_buffer_count, rcv_displs, IDX_T, comm);

        /* now after alltoallv, rcv_buffer has the node weights and node type */
        con_chunks_len[chunk] = isum(npes, rcv_buffer_count, 1);
        con_chunks[chunk] = imalloc(con_chunks_len[chunk], "con_chunks[chunk]");
        icopy(con_chunks_len[chunk], rcv_buffer, con_chunks[chunk]);
        printf("[Rank: %d] Copied node weights: %d into chunk: %d\n", mype, con_chunks_len[chunk], chunk );
      }//end of processing ind. file.

      GKASSERT(gnid == end_nid);
      printf("[Rank: %d] Finished reading file... with ending chunk: %d\n", mype, chunk);

    } //end of file names
    nchunks = chunk;

    lvtxs = isum(nchunks, con_chunks_len, 1)/ncon;
    GKASSERT(nvtxs == lvtxs);

    vwgt = graph->vwgt = imalloc(lvtxs*(ncon-1), "DDGLPart_ReadGraph vwgt");
    s_vwgt = imalloc(lvtxs*(ncon-1), "DDGLPart_ReadGraph vwgt");

    vtype = graph->vtype = imalloc(lvtxs, "DGLPart_ReadGraph: vtype");
    s_vtype = imalloc(lvtxs, "DGLPart_ReadGraph: vtype");

    /* Temp buffers used to sort and store the nodes in the correct order to mimic the cyclic ordering */
    idx_buffer = (edge_t *)gk_malloc(sizeof(edge_t)*lvtxs, "idx_buffer_lvtsx");

    /* sort the global nids */
    printf("[Rank: %d] Arranging data in graph objects ...\n", mype);
    nvtxs = 0;
    for(chunk=0; chunk<nchunks; chunk++){
      for(i=0; i<con_chunks_len[chunk]/ncon; i++){
        idx_buffer[nvtxs].u = con_id_chunks[chunk][i]; //actual global nids

        GKASSERT( idx_buffer[nvtxs].u < graph->gnvtxs );

        idx_buffer[nvtxs].v = nvtxs; //contiguous ids
        s_vtype[nvtxs] = con_chunks[chunk][ncon*i+0];
        for(j=1; j<ncon; j++){
          s_vwgt[(ncon-1)*nvtxs+j-1] = con_chunks[chunk][ncon*i+j];
	}
	nvtxs ++;
      }
      gk_free((void **)&con_chunks[chunk], LTERM);
      gk_free((void **)&con_id_chunks[chunk], LTERM);
    }
    gk_free((void **)&con_chunks_len, LTERM);

    /*
    * sort is based on edge_t.u
    * Here u is the idx and v is the global nid
    * we sort global nids, and use idx to re-arrange the node weights
    */
    printf("[Rank: %d] initiating sorting of edges...\n", mype);
    edgesorti(nvtxs, idx_buffer);
    printf("[Rank: %d] Sorting edges done...\n", mype);

    //Assert here to check global nids are sorted.
    for(i=1; i<nvtxs; i++){
      if(idx_buffer[i].u < idx_buffer[i-1].u)
        break;
    }
    GKASSERT(i == nvtxs);

    /*
    * shuffle node types and node weights accordingly now
    * order is defined by u's in the idx buffer
    */
    for(i=0; i<nvtxs; i++){
      vtype[i] = s_vtype[idx_buffer[i].v];
      for(j=1; j<ncon; j++){
        vwgt[(ncon-1)*i + j-1] = s_vwgt[(ncon-1)*idx_buffer[i].v+j-1];
      }
    }
    printf( "Done reading nodes: %d\n", nvtxs);

    /* done reading the node file */
    gk_fclose(fpin);
    gk_fclose(fpinaux);

    /* free memory */
    gk_free((void **)&s_vwgt, LTERM);
    gk_free((void **)&s_vtype, LTERM);
    gk_free((void **)&snd_buffer_count, LTERM);
    gk_free((void **)&rcv_buffer_count, LTERM);
    gk_free((void **)&snd_buffer, LTERM);
    gk_free((void **)&rcv_buffer, LTERM);
    gk_free((void **)&idx_buffer, LTERM);
      
    for (pe=0; pe<npes; pe++) {
      gk_free((void **)&con_buffers[pe], LTERM);
      gk_free((void **)&con_id_buffers[pe], LTERM);
    }
  
    gk_free((void **)&con_buffers_cpos, &con_buffers, LTERM);
  
    /* populate vwgt and create (vmptr, vmdata) */
    //ASSERT2(nvtxs == isum(nchunks, con_chunks_len, 1));

    GKASSERT(nvtxs == vtxdist[mype+1]-vtxdist[mype]);
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
/*! Sorts based on increasing <u, v> */
/*************************************************************************/
void edgesorti(size_t n, edge_t *base)
{
#define ikeyval_lt(a, b) \
  ((a)->u < (b)->u || \
   ((a)->u == (b)->u && (a)->v < (b)->v))
  GK_MKQSORT(edge_t, base, n, ikeyval_lt);
#undef ikeyval_lt
}

