#ifndef OPERATORS_H
#define OPERATORS_H

#include "bc.h"
#ifdef GPAW_CUDA
#include "cuda/gpaw-cuda.h"
#endif

typedef struct
{
  PyObject_HEAD
  bmgsstencil stencil;
  boundary_conditions* bc;
  MPI_Request recvreq[2];
  MPI_Request sendreq[2];
#ifdef GPAW_CUDA
  int cuda;
  bmgsstencil_gpu stencil_gpu;
  int alloc_blocks;
  double* buf_gpu;
  double *sendbuf;
  double *recvbuf;
  double *sendbuf_gpu;
  double *recvbuf_gpu;
#endif
} OperatorObject;

#endif
