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
  int rstatus, lcount;

  rstatus = MPI_Get_count(status, datatype, &lcount);
  *count = lcount;

  return rstatus;
}

int gkMPI_Send(void *buf, idx_t count, MPI_Datatype datatype, idx_t dest,
        idx_t tag, MPI_Comm comm)
{
  return MPI_Send(buf, count, datatype, dest, tag, comm);
}

int gkMPI_Recv(void *buf, idx_t count, MPI_Datatype datatype,
        idx_t source, idx_t tag, MPI_Comm comm, MPI_Status *status)
{
  return  MPI_Recv(buf, count, datatype, source, tag, comm, status);
}

int gkMPI_Isend(void *buf, idx_t count, MPI_Datatype datatype, idx_t dest,
        idx_t tag, MPI_Comm comm, MPI_Request *request)
{
  return MPI_Isend(buf, count, datatype, dest, tag, comm, request);
}

int gkMPI_Irecv(void *buf, idx_t count, MPI_Datatype datatype,
        idx_t source, idx_t tag, MPI_Comm comm, MPI_Request *request)
{
  return MPI_Irecv(buf, count, datatype, source, tag, comm, request);
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
  return MPI_Bcast(buffer, count, datatype, root, comm);
}

int gkMPI_Reduce(void *sendbuf, void *recvbuf, idx_t count,
        MPI_Datatype datatype, MPI_Op op, idx_t root, MPI_Comm comm)
{
  return MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
}

int gkMPI_Allreduce(void *sendbuf, void *recvbuf, idx_t count,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
}

int gkMPI_Scan(void *sendbuf, void *recvbuf, idx_t count,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  return MPI_Scan(sendbuf, recvbuf, count, datatype, op, comm);
}

int gkMPI_Allgather(void *sendbuf, idx_t sendcount,
        MPI_Datatype sendtype, void *recvbuf, idx_t recvcount,
        MPI_Datatype recvtype, MPI_Comm comm)
{
  return MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, 
             recvcount, recvtype, comm);
}

int gkMPI_Alltoall(void *sendbuf, idx_t sendcount,
        MPI_Datatype sendtype, void *recvbuf, idx_t recvcount,
        MPI_Datatype recvtype, MPI_Comm comm)
{
  return MPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount,
             recvtype, comm);
}

int gkMPI_Alltoallv(void *sendbuf, idx_t *sendcounts,
        idx_t *sdispls, MPI_Datatype sendtype, void *recvbuf, 
        idx_t *recvcounts, idx_t *rdispls, MPI_Datatype recvtype, 
        MPI_Comm comm)
{
#if IDXTYPEWIDTH == 32
  return MPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, 
               recvbuf, recvcounts, rdispls, recvtype, comm);
#else
  idx_t i; 
  int status, npes, *lsendcounts, *lsdispls, *lrecvcounts, *lrdispls;

  MPI_Comm_size(comm, &npes);

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
}

int gkMPI_Allgatherv(void *sendbuf, idx_t sendcount, MPI_Datatype sendtype, 
        void *recvbuf, idx_t *recvcounts, idx_t *rdispls, 
        MPI_Datatype recvtype, MPI_Comm comm)
{
#if IDXTYPEWIDTH == 32
  return MPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, 
               recvcounts, rdispls, recvtype, comm);
#else
  idx_t i; 
  int status, npes, *lrecvcounts, *lrdispls;
  int exceed_max_int;
  idx_t *child_recvcounts, *child_rdispls;
  idx_t *child_recvbuf;
  idx_t child_sendcount;
  int mype;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);
 
  exceed_max_int = 0; 
  for (i=0; i<npes; i++) {
    if (recvcounts[i] >= INT_MAX || rdispls[i] >= INT_MAX) {
        exceed_max_int = 1;
        break;
    }
  }

  if (exceed_max_int || sendcount >= INT_MAX) {
    // Length of message is larger than int, we use divide and conquer
    // Take care of idx_t is 64bit long
 
    child_recvcounts = gk_i64malloc(npes, "child_recvcounts");
    child_rdispls    = gk_i64malloc(npes+1, "child_displs");
    child_recvbuf = gk_i64malloc((recvcounts[npes-1]+rdispls[npes-1])/2+2, "child_recvbuf");
    
    // Allgather first half
    for (i=0; i<npes; i++) {
      child_recvcounts[i] = recvcounts[i]/2;
    }
    for (child_rdispls[0]=0, i=1; i<npes+1; i++) {
      child_rdispls[i] = child_rdispls[i-1] + child_recvcounts[i-1];
    }
    child_sendcount = sendcount/2;
    
    status = gkMPI_Allgatherv(sendbuf, child_sendcount,
      sendtype, child_recvbuf, child_recvcounts, child_rdispls, recvtype, comm);
 
    for (i=0; i<npes; i++) {
      memcpy((idx_t*)recvbuf+rdispls[i],
        (idx_t*)child_recvbuf+child_rdispls[i], sizeof(idx_t)*child_recvcounts[i]);
    }
 
    // Allgather last half
    for (i=0; i<npes; i++) {
      child_recvcounts[i] += recvcounts[i]%2;
    }
    for (child_rdispls[0]=0,i=1; i<npes+1; i++) {
      child_rdispls[i] = child_rdispls[i-1]+child_recvcounts[i-1];
    }
    child_sendcount += sendcount%2;
 
    status = gkMPI_Allgatherv((void*)((idx_t*)sendbuf+sendcount/2),
      child_sendcount, sendtype, child_recvbuf, child_recvcounts, child_rdispls, recvtype, comm);
 
    for (i=0; i<npes; i++) {
      memcpy((idx_t*)recvbuf+rdispls[i]+recvcounts[i]/2,
        (idx_t*)child_recvbuf+child_rdispls[i], sizeof(idx_t)*child_recvcounts[i]);
    }
 
    gk_free((void **)&child_recvcounts, &child_rdispls, &child_recvbuf, LTERM);
    return status;
  }

  // Length of message is shorter than int max, run allgatherv directly
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
}

int gkMPI_Scatterv(void *sendbuf, idx_t *sendcounts, idx_t *sdispls,
        MPI_Datatype sendtype, void *recvbuf, idx_t recvcount,
        MPI_Datatype recvtype, idx_t root, MPI_Comm comm)
{
#if IDXTYPEWIDTH == 32
  return MPI_Scatterv(sendbuf, sendcounts, sdispls, sendtype, 
               recvbuf, recvcount, recvtype, root, comm);
#else
  idx_t i; 
  int status, npes, *lsendcounts, *lsdispls;

  MPI_Comm_size(comm, &npes);

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
}

int gkMPI_Gatherv(void *sendbuf, idx_t sendcount, MPI_Datatype sendtype,
        void *recvbuf, idx_t *recvcounts, idx_t *rdispls, MPI_Datatype recvtype,
        idx_t root, MPI_Comm comm)
{
#if IDXTYPEWIDTH == 32
  return MPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, 
               recvcounts, rdispls, recvtype, root, comm);
#else
  idx_t i; 
  int status, npes, *lrecvcounts, *lrdispls;

  MPI_Comm_size(comm, &npes);

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
