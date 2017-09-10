/*
 * Copyright (c) 2017, NVIDIA CORPORATION. All rights reserved.
 *  
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif
#include <mpi.h>

#ifdef USE_CUDA
#define CHECK_CUDART(x) do { \
  cudaError_t res = (x); \
  if(res != cudaSuccess) { \
    fprintf(stderr, "CUDART: %s = %d (%s) at (%s:%d)\n",#x, res, cudaGetErrorString(res),__FILE__,__LINE__); \
    exit(1); \
  } \
} while(0)
#endif

#define MAXBUF (4)
#define MAXPEER (16)

#ifdef USE_CUDA
static cudaStream_t streams[MAXPEER];
static cudaEvent_t events[MAXPEER];
#endif
static int first_time = 1;
static double* buff_rem[MAXBUF][MAXPEER];
static double* buff_base[MAXBUF];
static int can_access_peer[MAXPEER];
static int gprocs;
static int grank;

void init_ipc_( double* buffer, int *buff_id_p, MPI_Fint *Fcomm, int* ipc_peers)
{
#ifdef USE_CUDA
  int buff_id = *buff_id_p;
  int rank, nprocs, i;
  cudaIpcMemHandle_t loc_handle;
  cudaIpcMemHandle_t rem_handle[MAXPEER];
  cudaError_t istat;

  MPI_Comm comm = MPI_Comm_f2c(*Fcomm);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);

  if(rank==0) printf("initializing IPC buffer: %d \n",buff_id);

  if(first_time){
    for(i=0; i<MAXPEER; i++) CHECK_CUDART( cudaStreamCreate( &streams[i] ) );
    for(i=0; i<MAXPEER; i++) CHECK_CUDART( cudaEventCreateWithFlags( &events[i], cudaEventDisableTiming ) );
    first_time = 0;
  }
  
  buff_base[ buff_id ] = buffer;
  CHECK_CUDART(cudaIpcGetMemHandle(&loc_handle, buffer));

  MPI_Allgather( &loc_handle, sizeof(cudaIpcMemHandle_t), MPI_CHAR,
                  rem_handle, sizeof(cudaIpcMemHandle_t), MPI_CHAR,
                  comm);

  for(i=0; i<nprocs; i++){
    if(i!=rank){
      istat = cudaIpcOpenMemHandle((void **)(&buff_rem[buff_id][i]), rem_handle[i], cudaIpcMemLazyEnablePeerAccess);
      if(istat == cudaSuccess){
        can_access_peer[i] = 1;
      }else if(istat == 64){
        can_access_peer[i] = 0;
      }else{
        printf("Error in cudaIpcOpenMemHandle: %d %s\n",istat, cudaGetErrorString(istat));
      }
    }
  }

  for(i=0; i<nprocs; i++){
    if(i!=rank){
      ipc_peers[i] = can_access_peer[i];
    }
  }
  gprocs = nprocs;
  grank = rank;
#endif
  return;
}

void get_ipc_peers_( int* ipc_peers)
{
  int i;
  for(i=0; i<gprocs; i++){
    if(i!=grank){
      ipc_peers[i] = can_access_peer[i];
    }
  }
  return;
}

void sync_ipc_sends_( MPI_Fint *Fcomm )
{
#ifdef USE_CUDA
  int rank, nprocs, i;
  MPI_Comm comm = MPI_Comm_f2c(*Fcomm);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);

  for(i=0; i<nprocs; i++){
    if(rank != i){
        if(can_access_peer[i] == 1)  CHECK_CUDART( cudaEventSynchronize( events[i] ) );
    }
  }
#endif
  return;
}

void ipc_send_( double *sendbuff, int *elems_p, double *recvbuff, int *buff_id_p, int *peer_id_p, MPI_Fint *Fcomm, int *ierr_p )
{
#ifdef USE_CUDA
  int buff_id = *buff_id_p;
  int elems = 2*(*elems_p);
  int dest = *peer_id_p;
  MPI_Comm comm = MPI_Comm_f2c(*Fcomm);
  int ierr;

  size_t offset = recvbuff - buff_base[buff_id];
  if(can_access_peer[dest] == 1){
    CHECK_CUDART( cudaMemcpyAsync( buff_rem[buff_id][dest] + offset, sendbuff, elems*sizeof(double), cudaMemcpyDefault, streams[dest] ) );
    CHECK_CUDART( cudaEventRecord( events[dest], streams[dest] ) );
  }else{
    printf("CANNOT ACCESS PEER!!!! \n");
    exit(1);
  }
#endif
  return;
}

#if 0
void my_a2a_( double *sendbuff, int *elems_p, double *recvbuff, int *buff_id_p, MPI_Fint *Fcomm, int *ierr_p )
{
  //MPI_Request rqst[16];
  //MPI_Status stts[16];
  int buff_id = *buff_id_p;
  int elems = *elems_p;
  MPI_Comm comm = MPI_Comm_f2c(*Fcomm);
  int ierr;
  int rank, nprocs;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);

  int i;
  //int rqst_count=0;

#if 0
//  ierr = MPI_Alltoall( sendbuff_p, elems, MPI_DOUBLE, recvbuff_p, elems, MPI_DOUBLE, MPI_COMM_WORLD );
  for(i=0; i<nprocs; i++){
    int partner = rank ^ i;
    //printf("iter %d rank %d <--> %d \n",i,rank,partner);
    if(rank != partner){
      MPI_Irecv( recvbuff_p + elems*partner, elems, MPI_DOUBLE, partner, 0, comm, &rqst[rqst_count++] );
      MPI_Isend( sendbuff_p + elems*partner, elems, MPI_DOUBLE, partner, 0, comm, &rqst[rqst_count++] ); 
    }
  }

  MPI_Waitall( rqst_count,  rqst, stts );
#else

    size_t offset = recvbuff - buff_base[buff_id]; 
/*
if(rank==0){
    if(smA<smB)
    printf("smA = %llu < smB = %llu \n",smA,smB);
    if(smB<smA)
    printf("smA = %llu > smB = %llu \n",smA,smB);
}
*/
//  CHECK_CUDART( cudaDeviceSynchronize() );

// need barrier to prevent writing to peer before peer is ready (no longer needed with aux2_d
// MPI_Barrier( comm );

  for(i=0; i<nprocs; i++){
//TODO: for non power-of-2 partner=(rank+i)%nprocs
//TODO: use mpi for non-ipc peers

    int partner = rank ^ i;
    //printf("iter %d rank %d <--> %d \n",i,rank,partner);
    if(rank != partner){
/*
      MPI_Irecv( recvbuff_p + elems*partner, elems, MPI_DOUBLE, partner, 0, comm, &rqst[rqst_count++] );
      MPI_Isend( sendbuff_p + elems*partner, elems, MPI_DOUBLE, partner, 0, comm, &rqst[rqst_count++] ); 
*/
/*
      if(smA<smB){
        CHECK_CUDART( cudaMemcpyAsync( B_rem[partner] + smA + elems*rank, sendbuff_p + elems*partner, elems*sizeof(double), cudaMemcpyDefault, streams[partner] ) );
      }else{
        CHECK_CUDART( cudaMemcpyAsync( A_rem[partner] + smB + elems*rank, sendbuff_p + elems*partner, elems*sizeof(double), cudaMemcpyDefault, streams[partner] ) );
      }
*/
        CHECK_CUDART( cudaMemcpyAsync( buff_rem[buff_id][partner] + offset + elems*rank, sendbuff + elems*partner, elems*sizeof(double), cudaMemcpyDefault, streams[partner] ) );
        CHECK_CUDART( cudaEventRecord( events[i], streams[partner] ) );
    }
  }

//  MPI_Waitall( rqst_count,  rqst, stts );
/*
  for(i=0; i<16; i++) CHECK_CUDART( cudaStreamDestroy( streams[i] ) );
  for(i=0; i<nprocs; i++){
    if(i!=rank){
      CHECK_CUDART( cudaIpcCloseMemHandle(rem_rcv_buf[i]) );
      CHECK_CUDART( cudaIpcCloseMemHandle(rem_snd_buf[i]) );
    }
  }
*/
//  CHECK_CUDART( cudaDeviceSynchronize() );
//  MPI_Barrier( comm );

#endif
  //ierr_p[0] = ierr;
  return;
}
#endif


/*

int main(int argc, const char *argv[])
{
  int rank, nprocs,name_len;
  char host_name[MPI_MAX_PROCESSOR_NAME];
  char (*host_names)[MPI_MAX_PROCESSOR_NAME];
  MPI_Comm nodeComm;
  MPI_Request rqst[16];
  MPI_Status stts[16];

//  profiler_on=1;

  MPI_Init(&argc, (char***)&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Get_processor_name(host_name,&name_len);
//  printf("%d %s HERE\n",rank,host_name);

  double *send_buffer_d, *recv_buffer_d;
  double *send_buffer_h, *recv_buffer_h;

  int bytes_per = 4*1024*1024;
  int elems_per = bytes_per/(sizeof(double));

  int elems = nprocs*elems_per;
  int bytes = nprocs*bytes_per;

  CHECK_CUDART( cudaMalloc( (void**)&send_buffer_d, bytes ) );
  CHECK_CUDART( cudaMalloc( (void**)&recv_buffer_d, bytes ) );
  CHECK_CUDART( cudaMallocHost( (void**)&send_buffer_h, bytes ) );
  CHECK_CUDART( cudaMallocHost( (void**)&recv_buffer_h, bytes ) );

  int i;
  for(i=0; i<elems; i++) send_buffer_h[i] = rank;
  for(i=0; i<elems; i++) recv_buffer_h[i] = 99;

  CHECK_CUDART( cudaMemcpy( send_buffer_d, send_buffer_h, bytes, cudaMemcpyHostToDevice ) );
  CHECK_CUDART( cudaMemcpy( recv_buffer_d, recv_buffer_h, bytes, cudaMemcpyHostToDevice ) );

  cudaStream_t streams[16];
  for(i=0; i<16; i++) CHECK_CUDART( cudaStreamCreate( &streams[i] ) );

  MPI_Barrier( MPI_COMM_WORLD );

  cudaIpcMemHandle_t rcv_loc;
  cudaIpcMemHandle_t snd_loc;
  cudaIpcMemHandle_t rcv_hdl[16];
  cudaIpcMemHandle_t snd_hdl[16];
  double* rem_rcv_buf[16];
  double* rem_snd_buf[16];

  CHECK_CUDART(cudaIpcGetMemHandle(&rcv_loc, recv_buffer_d));
  CHECK_CUDART(cudaIpcGetMemHandle(&snd_loc, send_buffer_d));

  MPI_Allgather( &rcv_loc, sizeof(cudaIpcMemHandle_t), MPI_CHAR,
                  rcv_hdl, sizeof(cudaIpcMemHandle_t), MPI_CHAR,
                  MPI_COMM_WORLD);

  MPI_Allgather( &snd_loc, sizeof(cudaIpcMemHandle_t), MPI_CHAR,
                  snd_hdl, sizeof(cudaIpcMemHandle_t), MPI_CHAR,
                  MPI_COMM_WORLD);

  for(i=0; i<nprocs; i++){
    if(i!=rank){
      CHECK_CUDART( cudaIpcOpenMemHandle((void **)(&rem_rcv_buf[i]), rcv_hdl[i], cudaIpcMemLazyEnablePeerAccess) );
      CHECK_CUDART( cudaIpcOpenMemHandle((void **)(&rem_snd_buf[i]), snd_hdl[i], cudaIpcMemLazyEnablePeerAccess) );
    }
  }


  double timer;
  int j;

for(j=0; j<ITERS; j++){
  MPI_Barrier( MPI_COMM_WORLD );
START_RANGE("IPC PUSH",j);
  timer = MPI_Wtime();
  for(i=0; i<nprocs; i++){
    int partner = rank ^ i;
    //printf("iter %d rank %d <--> %d \n",i,rank,partner);
    if(rank != partner){
      CHECK_CUDART( cudaMemcpyAsync( rem_rcv_buf[partner] + elems_per*i, send_buffer_d + elems_per*i, bytes_per, cudaMemcpyDefault, streams[partner] ) );
    }
  }
  CHECK_CUDART( cudaDeviceSynchronize() );
  MPI_Barrier( MPI_COMM_WORLD );
END_RANGE
  timer = MPI_Wtime() - timer;

  if(rank==0) printf("IPC PUSH: time = %3.3f ms, %3.1f GB/s \n",timer*1000.0,2*(nprocs-1)*bytes_per*1e-9/timer);
}


for(j=0; j<ITERS; j++){
  MPI_Barrier( MPI_COMM_WORLD );
START_RANGE("IPC PULL",j);
  timer = MPI_Wtime();
  for(i=0; i<nprocs; i++){
    int partner = rank ^ i;
    //printf("iter %d rank %d <--> %d \n",i,rank,partner);
    if(rank != partner){
      CHECK_CUDART( cudaMemcpyAsync( recv_buffer_d + elems_per*i, rem_snd_buf[partner] + elems_per*i, bytes_per, cudaMemcpyDefault, streams[partner] ) );
    }
  }
  CHECK_CUDART( cudaDeviceSynchronize() );
  MPI_Barrier( MPI_COMM_WORLD );
END_RANGE
  timer = MPI_Wtime() - timer;

  if(rank==0) printf("IPC PULL: time = %3.3f ms, %3.1f GB/s \n",timer*1000.0,2*(nprocs-1)*bytes_per*1e-9/timer);
}

for(j=0; j<ITERS; j++){
  MPI_Barrier( MPI_COMM_WORLD );
START_RANGE("MPI_COMM",j);
  timer = MPI_Wtime();
  int rqst_count=0;
  for(i=0; i<nprocs; i++){
    int partner = rank ^ i;
    //printf("iter %d rank %d <--> %d \n",i,rank,partner);
    if(rank != partner){
      MPI_Irecv( recv_buffer_d + elems_per*i, elems_per, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD, &rqst[rqst_count++] );
      MPI_Isend( send_buffer_d + elems_per*i, elems_per, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD, &rqst[rqst_count++] ); 
    }
  }

  MPI_Waitall( rqst_count,  rqst, stts );
END_RANGE
  timer = MPI_Wtime() - timer;
  
  if(rank==0) printf("MPI: time = %3.3f ms, %3.1f GB/s \n",timer*1000.0,rqst_count*bytes_per*1e-9/timer);
}



  MPI_Finalize();
  CHECK_CUDART( cudaDeviceReset() );
  return( 0 );
}

*/
