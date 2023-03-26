/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * gkmpi.c
 *
 * This function contains wrappers around MPI calls to allow
 * future de-coupling of sizes from 'int' datatypes.
 *
 * Started 5/30/11
 * George
 *
 * $Id: gkmpi.c 10022 2011-05-30 20:18:42Z karypis $
 */

#include <parmetislib.h>

int gkMPI_Comm_size(MPI_Comm comm, idx_t *size)
{
  int status, lsize;

  status = MPI_Comm_size(comm, &lsize);
  *size = lsize;

  return status;
}

int gkMPI_Comm_rank(MPI_Comm comm, idx_t *rank)
{
  int status, lrank;

  status = MPI_Comm_rank(comm, &lrank);
  *rank = lrank;

  return status;
}
    
int gkMPI_Get_count(MPI_Status *status, MPI_Datatype datatype,
        idx_t *count)
{
  int rstatus;

#if MPI_VERSION < 4
  int lcount;
  rstatus = MPI_Get_count(status, datatype, &lcount);
#else
  MPI_Count lcount;
  rstatus = MPI_Get_count_c(status, datatype, &lcount);
#endif

  *count = lcount;

  return rstatus;
}

int gkMPI_Send(void *buf, idx_t count, MPI_Datatype datatype, idx_t dest,
        idx_t tag, MPI_Comm comm)
{
#if MPI_VERSION < 4
  return MPI_Send(buf, count, datatype, dest, tag, comm);
#else
  return MPI_Send_c(buf, count, datatype, dest, tag, comm);
#endif
}

int gkMPI_Recv(void *buf, idx_t count, MPI_Datatype datatype,
        idx_t source, idx_t tag, MPI_Comm comm, MPI_Status *status)
{
#if MPI_VERSION < 4
  return MPI_Recv(buf, count, datatype, source, tag, comm, status);
#else
  return MPI_Recv_c(buf, count, datatype, source, tag, comm, status);
#endif
}

int gkMPI_Isend(void *buf, idx_t count, MPI_Datatype datatype, idx_t dest,
        idx_t tag, MPI_Comm comm, MPI_Request *request)
{
#if MPI_VERSION < 4
  return MPI_Isend(buf, count, datatype, dest, tag, comm, request);
#else
  return MPI_Isend_c(buf, count, datatype, dest, tag, comm, request);
#endif
}

int gkMPI_Irecv(void *buf, idx_t count, MPI_Datatype datatype,
        idx_t source, idx_t tag, MPI_Comm comm, MPI_Request *request)
{
#if MPI_VERSION < 4
  return MPI_Irecv(buf, count, datatype, source, tag, comm, request);
#else
  return MPI_Irecv_c(buf, count, datatype, source, tag, comm, request);
#endif
}

int gkMPI_Wait(MPI_Request *request, MPI_Status *status)
{
  return MPI_Wait(request, status);
}

int gkMPI_Waitall(idx_t count, MPI_Request *array_of_requests, 
        MPI_Status *array_of_statuses)
{
  return MPI_Waitall(count, array_of_requests, array_of_statuses);
}

int gkMPI_Barrier(MPI_Comm comm)
{
  return MPI_Barrier(comm);
}

int gkMPI_Bcast(void *buffer, idx_t count, MPI_Datatype datatype,
        idx_t root, MPI_Comm comm)
{
#if MPI_VERSION < 4
  return MPI_Bcast(buffer, count, datatype, root, comm);
#else
  return MPI_Bcast_c(buffer, count, datatype, root, comm);
#endif
}

int gkMPI_Reduce(void *sendbuf, void *recvbuf, idx_t count,
        MPI_Datatype datatype, MPI_Op op, idx_t root, MPI_Comm comm)
{
#if MPI_VERSION < 4
  return MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
#else
  return MPI_Reduce_c(sendbuf, recvbuf, count, datatype, op, root, comm);
#endif
}

int gkMPI_Allreduce(void *sendbuf, void *recvbuf, idx_t count,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
#if MPI_VERSION < 4
  return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
#else
  return MPI_Allreduce_c(sendbuf, recvbuf, count, datatype, op, comm);
#endif
}

int gkMPI_Scan(void *sendbuf, void *recvbuf, idx_t count,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
#if MPI_VERSION < 4
  return MPI_Scan(sendbuf, recvbuf, count, datatype, op, comm);
#else
  return MPI_Scan_c(sendbuf, recvbuf, count, datatype, op, comm);
#endif
}

int gkMPI_Allgather(void *sendbuf, idx_t sendcount,
        MPI_Datatype sendtype, void *recvbuf, idx_t recvcount,
        MPI_Datatype recvtype, MPI_Comm comm)
{
#if MPI_VERSION < 4
  return MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, 
             recvcount, recvtype, comm);
#else
  return MPI_Allgather_c(sendbuf, sendcount, sendtype, recvbuf, 
             recvcount, recvtype, comm);
#endif
}

int gkMPI_Alltoall(void *sendbuf, idx_t sendcount,
        MPI_Datatype sendtype, void *recvbuf, idx_t recvcount,
        MPI_Datatype recvtype, MPI_Comm comm)
{
#if MPI_VERSION < 4
  return MPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount,
             recvtype, comm);
#else
  return MPI_Alltoall_c(sendbuf, sendcount, sendtype, recvbuf, recvcount,
             recvtype, comm);
#endif
}

int gkMPI_Alltoallv(void *sendbuf, idx_t *sendcounts,
        idx_t *sdispls, MPI_Datatype sendtype, void *recvbuf, 
        idx_t *recvcounts, idx_t *rdispls, MPI_Datatype recvtype, 
        MPI_Comm comm)
{
#if MPI_VERSION < 4
#if IDXTYPEWIDTH == 32
  return MPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, 
               recvbuf, recvcounts, rdispls, recvtype, comm);
#else
  idx_t i; 
  int status, npes, *lsendcounts, *lsdispls, *lrecvcounts, *lrdispls;

  MPI_Comm_size(comm, &npes);

  /* bail-out if MPI 3.x cannot handle such large counts */
  for (i=0; i<npes; i++) { 
    if (sendcounts[i] >= INT_MAX || sdispls[i] >= INT_MAX || 
        recvcounts[i] >= INT_MAX || rdispls[i] >= INT_MAX)
      errexit("MPI_Gatherv message sizes goes over INT_MAX. Use MPI 4.x\n");
  }
    
  lsendcounts = gk_imalloc(npes, "lsendcounts");
  lsdispls    = gk_imalloc(npes, "lsdispls");
  lrecvcounts = gk_imalloc(npes, "lrecvcounts");
  lrdispls    = gk_imalloc(npes, "lrdispls");

  for (i=0; i<npes; i++) {
    lsendcounts[i] = sendcounts[i];
    lsdispls[i]    = sdispls[i];
    lrecvcounts[i] = recvcounts[i];
    lrdispls[i]    = rdispls[i];
  }

  status = MPI_Alltoallv(sendbuf, lsendcounts, lsdispls, sendtype, 
               recvbuf, lrecvcounts, lrdispls, recvtype, comm);

  for (i=0; i<npes; i++) {
    sendcounts[i] = lsendcounts[i];
    sdispls[i]    = lsdispls[i];
    recvcounts[i] = lrecvcounts[i];
    rdispls[i]    = lrdispls[i];
  }

  gk_free((void **)&lsendcounts, &lrecvcounts, &lsdispls, &lrdispls, LTERM);

  return status;
#endif
#else
#if IDXTYPEWIDTH == 32
  return MPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, 
               recvbuf, recvcounts, rdispls, recvtype, comm);
#else
  return MPI_Alltoallv_c(sendbuf, sendcounts, sdispls, sendtype, 
               recvbuf, recvcounts, rdispls, recvtype, comm);
#endif
#endif
}

int gkMPI_Allgatherv(void *sendbuf, idx_t sendcount, MPI_Datatype sendtype, 
        void *recvbuf, idx_t *recvcounts, idx_t *rdispls, 
        MPI_Datatype recvtype, MPI_Comm comm)
{
#if MPI_VERSION < 4
#if IDXTYPEWIDTH == 32
  return MPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, 
               recvcounts, rdispls, recvtype, comm);
#else
  idx_t i; 
  int status, npes, *lrecvcounts, *lrdispls;

  MPI_Comm_size(comm, &npes);

  /* bail-out if MPI 3.x cannot handle such large counts */
  for (i=0; i<npes; i++) { 
    if (sendcount >= INT_MAX || recvcounts[i] >= INT_MAX || rdispls[i] >= INT_MAX)
      errexit("MPI_Allgatherv message sizes goes over INT_MAX. Use MPI 4.x\n");
  }

  lrecvcounts = gk_imalloc(npes, "lrecvcounts");
  lrdispls    = gk_imalloc(npes, "lrdispls");

  for (i=0; i<npes; i++) {
    lrecvcounts[i] = recvcounts[i];
    lrdispls[i]    = rdispls[i];
  }

  status = MPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, 
               lrecvcounts, lrdispls, recvtype, comm);

  for (i=0; i<npes; i++) {
    recvcounts[i] = lrecvcounts[i];
    rdispls[i]    = lrdispls[i];
  }

  gk_free((void **)&lrecvcounts, &lrdispls, LTERM);

  return status;
#endif
#else
#if IDXTYPEWIDTH == 32
  return MPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, 
               recvcounts, rdispls, recvtype, comm);
#else
  return MPI_Allgatherv_c(sendbuf, sendcount, sendtype, recvbuf, 
               recvcounts, rdispls, recvtype, comm);
#endif
#endif
}

int gkMPI_Scatterv(void *sendbuf, idx_t *sendcounts, idx_t *sdispls,
        MPI_Datatype sendtype, void *recvbuf, idx_t recvcount,
        MPI_Datatype recvtype, idx_t root, MPI_Comm comm)
{
#if MPI_VERSION < 4
#if IDXTYPEWIDTH == 32
  return MPI_Scatterv(sendbuf, sendcounts, sdispls, sendtype, 
               recvbuf, recvcount, recvtype, root, comm);
#else
  idx_t i; 
  int status, npes, *lsendcounts, *lsdispls;

  MPI_Comm_size(comm, &npes);

  /* bail-out if MPI 3.x cannot handle such large counts */
  for (i=0; i<npes; i++) { 
    if (sendcounts[i] >= INT_MAX || recvcount >= INT_MAX || sdispls[i] >= INT_MAX)
      errexit("MPI_Scatterv message sizes goes over INT_MAX. Use MPI 4.x\n");
  }

  lsendcounts = gk_imalloc(npes, "lsendcounts");
  lsdispls    = gk_imalloc(npes, "lsdispls");

  for (i=0; i<npes; i++) {
    lsendcounts[i] = sendcounts[i];
    lsdispls[i]    = sdispls[i];
  }

  status = MPI_Scatterv(sendbuf, lsendcounts, lsdispls, sendtype, 
               recvbuf, recvcount, recvtype, root, comm);

  for (i=0; i<npes; i++) {
    sendcounts[i] = lsendcounts[i];
    sdispls[i]    = lsdispls[i];
  }

  gk_free((void **)&lsendcounts, &lsdispls, LTERM);

  return status;
#endif
#else
#if IDXTYPEWIDTH == 32
  return MPI_Scatterv(sendbuf, sendcounts, sdispls, sendtype, 
               recvbuf, recvcount, recvtype, root, comm);
#else
  return MPI_Scatterv_c(sendbuf, sendcounts, sdispls, sendtype, 
               recvbuf, recvcount, recvtype, root, comm);
#endif
#endif
}

int gkMPI_Gatherv(void *sendbuf, idx_t sendcount, MPI_Datatype sendtype,
        void *recvbuf, idx_t *recvcounts, idx_t *rdispls, MPI_Datatype recvtype,
        idx_t root, MPI_Comm comm)
{
#if MPI_VERSION < 4
#if IDXTYPEWIDTH == 32
  return MPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, 
               recvcounts, rdispls, recvtype, root, comm);
#else
  idx_t i; 
  int status, npes, *lrecvcounts, *lrdispls;

  MPI_Comm_size(comm, &npes);

  /* bail-out if MPI 3.x cannot handle such large counts */
  for (i=0; i<npes; i++) { 
    if (sendcount >= INT_MAX || recvcounts[i] >= INT_MAX || rdispls[i] >= INT_MAX)
      errexit("MPI_Gatherv message sizes goes over INT_MAX. Use MPI 4.x\n");
  }


  lrecvcounts = gk_imalloc(npes, "lrecvcounts");
  lrdispls    = gk_imalloc(npes, "lrdispls");

  for (i=0; i<npes; i++) {
    lrecvcounts[i] = recvcounts[i];
    lrdispls[i]    = rdispls[i];
  }

  status = MPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, 
               lrecvcounts, lrdispls, recvtype, root, comm);

  for (i=0; i<npes; i++) {
    recvcounts[i] = lrecvcounts[i];
    rdispls[i]    = lrdispls[i];
  }

  gk_free((void **)&lrecvcounts, &lrdispls, LTERM);

  return status;
#endif
#else
#if IDXTYPEWIDTH == 32
  return MPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, 
               recvcounts, rdispls, recvtype, root, comm);
#else
  return MPI_Gatherv_c(sendbuf, sendcount, sendtype, recvbuf, 
               recvcounts, rdispls, recvtype, root, comm);
#endif
#endif
}

int gkMPI_Comm_split(MPI_Comm comm, idx_t color, idx_t key,
        MPI_Comm *newcomm)
{
  return MPI_Comm_split(comm, color, key, newcomm);
}

int gkMPI_Comm_free(MPI_Comm *comm)
{
  return MPI_Comm_free(comm);
}

int gkMPI_Finalize()
{
  return MPI_Finalize();
}
