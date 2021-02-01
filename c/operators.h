#ifndef __OPERATORS_H
#define __OPERATORS_H

/*  Copyright (C) 2009-2012  CSC - IT Center for Science Ltd.
 *  Please see the accompanying LICENSE file for further information. */

#include "bc.h"

#ifdef GPAW_CUDA
#include "cuda/gpaw-cuda.h"
#endif

#ifdef __OPERATORS_C
typedef struct
{
  PyObject_HEAD
  bmgsstencil stencil;
  boundary_conditions* bc;
  MPI_Request recvreq[2];
  MPI_Request sendreq[2];
  int nthreads;
#ifdef GPAW_CUDA
  int cuda;
  bmgsstencil_gpu stencil_gpu;
#endif
} OperatorObject;
#else
// Provide opaque type for routines outside operators.c
struct _OperatorObject;
typedef struct _OperatorObject OperatorObject;
#endif

#ifdef GPAW_CUDA
void operator_init_cuda(OperatorObject *self);
void operator_dealloc_cuda(int force);
#endif

void apply_worker(OperatorObject *self, int chunksize, int start,
		  int end, int thread_id, int nthreads,
		  const double* in, double* out,
		  bool real, const double_complex* ph);

#endif
