#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL GPAW_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>
#include <stdlib.h>
#include <pthread.h>

#include "../extensions.h"
#define __OPERATORS_C
#include "../operators.h"
#undef __OPERATORS_C
#include "bmgs.h"
#include "gpu.h"

extern int gpaw_cuda_debug;

#define OPERATOR_NSTREAMS (2)

static gpuStream_t operator_stream[OPERATOR_NSTREAMS];
static gpuEvent_t operator_event[2];
static int operator_streams = 0;

static double *operator_buf_gpu = NULL;
static int operator_buf_size = 0;
static int operator_init_count = 0;

static int debug_size_arr = 0;
static int debug_size_buf = 0;
static double *debug_sendbuf;
static double *debug_recvbuf;
static double *debug_buf_cpu;
static double *debug_buf_gpu;
static double *debug_out_cpu;
static double *debug_out_gpu;
static double *debug_in_cpu;

/*
 * Increment reference count to register a new operator object
 * and copy the stencil to the GPU.
 */
void operator_init_cuda(OperatorObject *self)
{
    self->stencil_gpu = bmgs_stencil_to_gpu(&(self->stencil));
    operator_init_count++;
}

/*
 * Ensure buffer is allocated and is big enough. Reallocate only if
 * size has increased.
 *
 * Create also CUDA streams and events if not already created.
 */
void operator_alloc_buffers(OperatorObject *self, int blocks)
{
    const boundary_conditions* bc = self->bc;
    const int* size2 = bc->size2;
    int ng2 = (bc->ndouble * size2[0] * size2[1] * size2[2]) * blocks;

    if (ng2 > operator_buf_size) {
        gpuFree(operator_buf_gpu);
        gpuCheckLastError();
        gpuMalloc(&operator_buf_gpu, sizeof(double) * ng2);
        operator_buf_size = ng2;
    }
    if (!operator_streams) {
        for (int i=0; i < OPERATOR_NSTREAMS; i++) {
            gpuStreamCreate(&(operator_stream[i]));
        }
        for (int i=0; i < 2; i++) {
            gpuEventCreateWithFlags(
                    &operator_event[i],
                    gpuEventDefault|gpuEventDisableTiming);
        }
        operator_streams = OPERATOR_NSTREAMS;
    }
}

/*
 * Reset reference count and unset buffer.
 */
void operator_init_buffers_cuda()
{
    operator_buf_gpu = NULL;
    operator_buf_size = 0;
    operator_init_count = 0;
    operator_streams = 0;
}

/*
 * Deallocate buffer and destroy CUDA streams and events,
 * or decrease reference count
 *
 * arguments:
 *   (int) force -- if true, force deallocation etc.
 */
void operator_dealloc_cuda(int force)
{
    if (force) {
        operator_init_count = 1;
    }
    if (operator_init_count == 1) {
        gpuFree(operator_buf_gpu);
        if (operator_streams) {
            for (int i=0; i < OPERATOR_NSTREAMS; i++) {
                gpuStreamSynchronize(operator_stream[i]);
                gpuStreamDestroy(operator_stream[i]);
            }
            for (int i=0; i < 2; i++) {
                gpuEventDestroy(operator_event[i]);
            }
        }
        operator_init_buffers_cuda();
        return;
    }
    if (operator_init_count > 0) {
        operator_init_count--;
    }
}

/*
 * Allocate debug buffers and precalculate sizes.
 */
void debug_operator_allocate(OperatorObject* self, int nin, int blocks)
{
    boundary_conditions* bc = self->bc;
    const int *size1 = bc->size1;
    const int *size2 = bc->size2;
    int ng = bc->ndouble * size1[0] * size1[1] * size1[2];
    int ng2 = bc->ndouble * size2[0] * size2[1] * size2[2];

    debug_size_arr = ng * nin;
    debug_size_buf = ng2 * blocks;

    debug_sendbuf = GPAW_MALLOC(double, bc->maxsend * blocks);
    debug_recvbuf = GPAW_MALLOC(double, bc->maxrecv * blocks);
    debug_buf_cpu = GPAW_MALLOC(double, debug_size_buf);
    debug_buf_gpu = GPAW_MALLOC(double, debug_size_buf);
    debug_out_cpu = GPAW_MALLOC(double, debug_size_arr);
    debug_out_gpu = GPAW_MALLOC(double, debug_size_arr);
    debug_in_cpu = GPAW_MALLOC(double, debug_size_arr);
}

/*
 * Deallocate debug buffers and set sizes to zero.
 */
void debug_operator_deallocate()
{
    free(debug_sendbuf);
    free(debug_recvbuf);
    free(debug_buf_cpu);
    free(debug_buf_gpu);
    free(debug_out_cpu);
    free(debug_out_gpu);
    free(debug_in_cpu);
    debug_size_arr = 0;
    debug_size_buf = 0;
}

/*
 * Copy initial GPU arrays to debug buffers on the CPU.
 */
void debug_operator_memcpy_pre(const double *in, double *out)
{
    gpuMemcpy(debug_in_cpu, in, sizeof(double) * debug_size_arr,
              gpuMemcpyDeviceToHost);
    gpuMemcpy(debug_out_cpu, out, sizeof(double) * debug_size_arr,
              gpuMemcpyDeviceToHost);
}

/*
 * Copy final GPU arrays to debug buffers on the CPU.
 */
void debug_operator_memcpy_post(double *out, double *buf)
{
    gpuMemcpy(debug_out_gpu, out, sizeof(double) * debug_size_arr,
              gpuMemcpyDeviceToHost);
    gpuMemcpy(debug_buf_gpu, buf, sizeof(double) * debug_size_buf,
              gpuMemcpyDeviceToHost);
}

/*
 * Run the relax algorithm (see Operator_relax() in ../operators.c)
 * on the CPU and compare to results from the GPU.
 */
void debug_operator_relax(OperatorObject* self, int relax_method, int nrelax,
                          double w)
{
    MPI_Request recvreq[2];
    MPI_Request sendreq[2];
    int i;

    boundary_conditions* bc = self->bc;

    const double_complex *ph;
    ph = 0;

    for (int n=0; n < nrelax; n++ ) {
        for (i=0; i < 3; i++) {
            bc_unpack1(bc, debug_out_cpu, debug_buf_cpu, i, recvreq, sendreq,
                       debug_recvbuf, debug_sendbuf, ph + 2 * i, 0, 1);
            bc_unpack2(bc, debug_buf_cpu, i, recvreq, sendreq,
                       debug_recvbuf, 1);
        }
        bmgs_relax(relax_method, &self->stencil, debug_buf_cpu, debug_out_cpu,
                   debug_in_cpu, w);
    }

    double buf_err = 0;
    for (i=0; i < debug_size_buf; i++) {
        buf_err = MAX(buf_err, fabs(debug_buf_cpu[i] - debug_buf_gpu[i]));
    }
    double fun_err = 0;
    for (i=0; i < debug_size_arr; i++) {
        fun_err = MAX(fun_err, fabs(debug_out_cpu[i] - debug_out_gpu[i]));
    }
    int rank = 0;
    if (bc->comm != MPI_COMM_NULL)
        MPI_Comm_rank(bc->comm, &rank);
    if (buf_err > GPU_ERROR_ABS_TOL) {
        fprintf(stderr,
                "[%d] Debug CUDA operator relax (buf): error %g\n",
                rank, buf_err);
    }
    if (fun_err > GPU_ERROR_ABS_TOL) {
        fprintf(stderr,
                "[%d] Debug CUDA operator relax (fun): error %g\n",
                rank, fun_err);
    }
}

/*
 * Run the relax algorithm (see Operator_relax() in ../operators.c)
 * on the GPU.
 */
static void _operator_relax_cuda_gpu(OperatorObject* self, int relax_method,
                                     double *fun, const double *src,
                                     int nrelax, double w)
{
    boundary_conditions* bc = self->bc;
    const int *size2 = bc->size2;
    const int *size1 = bc->size1;
    int ng = bc->ndouble * size1[0] * size1[1] * size1[2];
    int ng2 = bc->ndouble * size2[0] * size2[1] * size2[2];

    MPI_Request recvreq[3][2];
    MPI_Request sendreq[3][2];

    const double_complex *ph;
    ph = 0;

    int blocks = 1;
    operator_alloc_buffers(self, blocks);

    int boundary = 0;
    if (bc->sendproc[0][0] != DO_NOTHING)
        boundary |= GPAW_BOUNDARY_X0;
    if (bc->sendproc[0][1] != DO_NOTHING)
        boundary |= GPAW_BOUNDARY_X1;
    if (bc->sendproc[1][0] != DO_NOTHING)
        boundary |= GPAW_BOUNDARY_Y0;
    if (bc->sendproc[1][1] != DO_NOTHING)
        boundary |= GPAW_BOUNDARY_Y1;
    if (bc->sendproc[2][0] != DO_NOTHING)
        boundary |= GPAW_BOUNDARY_Z0;
    if (bc->sendproc[2][1] != DO_NOTHING)
        boundary |= GPAW_BOUNDARY_Z1;

    int cuda_overlap = bmgs_fd_boundary_test(&self->stencil_gpu, boundary,
                                             bc->ndouble);
    int nsendrecvs = 0;
    for (int i=0; i < 3; i++) {
        for (int j=0; j < 2; j++) {
            nsendrecvs += MAX(bc->nsend[i][j], bc->nrecv[i][j])
                        * blocks * sizeof(double);
        }
    }
    cuda_overlap &= (nsendrecvs > GPU_OVERLAP_SIZE);
    if (cuda_overlap)
        gpuEventRecord(operator_event[1], 0);

    for (int n=0; n < nrelax; n++ ) {
        if (cuda_overlap) {
            gpuStreamWaitEvent(operator_stream[0], operator_event[1], 0);
            bc_unpack_paste_cuda_gpu(bc, fun, operator_buf_gpu, recvreq,
                                     operator_stream[0], 1);
            gpuEventRecord(operator_event[0], operator_stream[0]);

            bmgs_relax_cuda_gpu(relax_method, &self->stencil_gpu,
                                operator_buf_gpu, fun, src, w,
                                boundary|GPAW_BOUNDARY_SKIP,
                                operator_stream[0]);
            gpuStreamWaitEvent(operator_stream[1], operator_event[0], 0);
            for (int i=0; i < 3; i++) {
                bc_unpack_cuda_gpu_async(bc, fun, operator_buf_gpu, i,
                                         recvreq, sendreq[i], ph + 2 * i,
                                         operator_stream[1], 1);
            }
            bmgs_relax_cuda_gpu(relax_method, &self->stencil_gpu,
                                operator_buf_gpu, fun, src, w,
                                boundary|GPAW_BOUNDARY_ONLY,
                                operator_stream[1]);
            gpuEventRecord(operator_event[1], operator_stream[1]);
        } else {
            bc_unpack_paste_cuda_gpu(bc, fun, operator_buf_gpu, recvreq,
                                     0, 1);
            for (int i=0; i < 3; i++) {
                bc_unpack_cuda_gpu(bc, fun, operator_buf_gpu, i,
                                   recvreq, sendreq[i], ph + 2 * i, 0, 1);
            }
            bmgs_relax_cuda_gpu(relax_method, &self->stencil_gpu,
                                operator_buf_gpu, fun, src, w,
                                GPAW_BOUNDARY_NORMAL, 0);
        }
    }

    if (cuda_overlap) {
        gpuStreamWaitEvent(0, operator_event[1], 0);
        gpuStreamSynchronize(operator_stream[0]);
    }
}

/*
 * Python interface for the GPU version of the relax algorithm
 * (similar to Operator_relax() for CPUs).
 *
 * arguments:
 *   relax_method -- relaxation method (int)
 *   func_gpu     -- pointer to device memory (GPUArray.gpudata)
 *   source_gpu   -- pointer to device memory (GPUArray.gpudata)
 *   nrelax       -- number of iterations (int)
 *   w            -- weight (float)
 */
PyObject* Operator_relax_cuda_gpu(OperatorObject* self, PyObject* args)
{
    int relax_method;
    void *func_gpu;
    void *source_gpu;
    double w = 1.0;
    int nrelax;

    if (!PyArg_ParseTuple(args, "inni|d", &relax_method, &func_gpu,
                          &source_gpu, &nrelax, &w))
        return NULL;

    double *fun = (double*) func_gpu;
    const double *src = (double*) source_gpu;

    if (gpaw_cuda_debug) {
        debug_operator_allocate(self, 1, 1);
        debug_operator_memcpy_pre(src, fun);
    }

    _operator_relax_cuda_gpu(self, relax_method, fun, src, nrelax, w);

    if (gpaw_cuda_debug) {
        gpuDeviceSynchronize();
        debug_operator_memcpy_post(fun, operator_buf_gpu);
        debug_operator_relax(self, relax_method, nrelax, w);
        debug_operator_deallocate();
    }

    if (PyErr_Occurred())
        return NULL;
    else
        Py_RETURN_NONE;
}

/*
 * Run the FD algorithm (see apply_worker() in ../operators.c)
 * on the CPU and compare to results from the GPU.
 */
void debug_operator_apply(OperatorObject* self, int nin, int blocks,
                          bool real, const double_complex *ph)
{
    MPI_Request recvreq[2];
    MPI_Request sendreq[2];
    int i;
    double err;

    boundary_conditions* bc = self->bc;
    const int *size1 = bc->size1;
    const int *size2 = bc->size2;
    int ng = bc->ndouble * size1[0] * size1[1] * size1[2];
    int ng2 = bc->ndouble * size2[0] * size2[1] * size2[2];

    for (int n=0; n < nin; n += blocks) {
        const double *in = debug_in_cpu + n * ng;
        double *out = debug_out_cpu + n * ng;
        int myblocks = MIN(blocks, nin - n);
        for (i=0; i < 3; i++) {
            bc_unpack1(bc, in, debug_buf_cpu, i, recvreq, sendreq,
                       debug_recvbuf, debug_sendbuf, ph + 2 * i, 0, myblocks);
            bc_unpack2(bc, debug_buf_cpu, i, recvreq, sendreq,
                       debug_recvbuf, myblocks);
        }
        for (int m=0; m < myblocks; m++) {
            if (real)
                bmgs_fd(&self->stencil, debug_buf_cpu + m * ng2,
                        out + m * ng);
            else
                bmgs_fdz(&self->stencil,
                         (const double_complex*) (debug_buf_cpu + m * ng2),
                         (double_complex*) (out + m * ng));
        }
    }

    double buf_err = 0;
    int buf_err_n = 0;
    for (i = 0; i < debug_size_buf; i++) {
        err = fabs(debug_buf_cpu[i] - debug_buf_gpu[i]);
        if (err > GPU_ERROR_ABS_TOL)
            buf_err_n++;
        buf_err = MAX(buf_err, err);
    }
    double out_err = 0;
    int out_err_n = 0;
    for (i = 0; i < debug_size_arr; i++) {
        err = fabs(debug_out_cpu[i] - debug_out_gpu[i]);
        if (err > GPU_ERROR_ABS_TOL)
            out_err_n++;
        out_err = MAX(out_err, err);
    }
    int rank = 0;
    if (bc->comm != MPI_COMM_NULL)
        MPI_Comm_rank(bc->comm, &rank);
    if (buf_err > GPU_ERROR_ABS_TOL) {
        fprintf(stderr,
                "[%d] Debug CUDA operator apply (buf): error %g (count %d/%d)\n",
                rank, buf_err, buf_err_n, debug_size_buf);
    }
    if (out_err > GPU_ERROR_ABS_TOL) {
        fprintf(stderr,
                "[%d] Debug CUDA operator apply (out): error %g (count %d/%d)\n",
                rank, out_err, out_err_n, debug_size_arr);
    }
}

/*
 * Run the FD algorithm (see apply_worker() in ../operators.c)
 * on the GPU.
 */
static void _operator_apply_cuda_gpu(OperatorObject *self,
                                     const double *in, double *out,
                                     int nin, int blocks, bool real,
                                     const double_complex *ph)
{
    boundary_conditions* bc = self->bc;
    const int *size1 = bc->size1;
    const int *size2 = bc->size2;
    int ng = bc->ndouble * size1[0] * size1[1] * size1[2];
    int ng2 = bc->ndouble * size2[0] * size2[1] * size2[2];

    MPI_Request recvreq[3][2];
    MPI_Request sendreq[3][2];

    operator_alloc_buffers(self, blocks);

    int boundary = 0;
    if (bc->sendproc[0][0] != DO_NOTHING)
        boundary |= GPAW_BOUNDARY_X0;
    if (bc->sendproc[0][1] != DO_NOTHING)
        boundary |= GPAW_BOUNDARY_X1;
    if (bc->sendproc[1][0] != DO_NOTHING)
        boundary |= GPAW_BOUNDARY_Y0;
    if (bc->sendproc[1][1] != DO_NOTHING)
        boundary |= GPAW_BOUNDARY_Y1;
    if (bc->sendproc[2][0] != DO_NOTHING)
        boundary |= GPAW_BOUNDARY_Z0;
    if (bc->sendproc[2][1] != DO_NOTHING)
        boundary |= GPAW_BOUNDARY_Z1;

    int cuda_overlap = bmgs_fd_boundary_test(&self->stencil_gpu, boundary,
                                             bc->ndouble);
    int nsendrecvs = 0;
    for (int i=0; i < 3; i++) {
        for (int j=0; j < 2; j++) {
            nsendrecvs += MAX(bc->nsend[i][j], bc->nrecv[i][j])
                        * blocks * sizeof(double);
        }
    }
    cuda_overlap &= (nsendrecvs > GPU_OVERLAP_SIZE);
    if  (cuda_overlap)
        gpuEventRecord(operator_event[1], 0);

    for (int n=0; n < nin; n += blocks) {
        const double *in2 = in + n * ng;
        double *out2 = out + n * ng;
        int myblocks = MIN(blocks, nin - n);
        if (cuda_overlap) {
            gpuStreamWaitEvent(operator_stream[0], operator_event[1], 0);
            bc_unpack_paste_cuda_gpu(bc, in2, operator_buf_gpu, recvreq,
                                     operator_stream[0], myblocks);
            gpuEventRecord(operator_event[0], operator_stream[0]);

            if (real) {
                bmgs_fd_cuda_gpu(&self->stencil_gpu, operator_buf_gpu, out2,
                                 boundary|GPAW_BOUNDARY_SKIP, myblocks,
                                 operator_stream[0]);
            } else {
                bmgs_fd_cuda_gpuz(&self->stencil_gpu,
                                  (const gpuDoubleComplex*) operator_buf_gpu,
                                  (gpuDoubleComplex*) out2,
                                  boundary|GPAW_BOUNDARY_SKIP, myblocks,
                                  operator_stream[0]);
            }
            gpuStreamWaitEvent(operator_stream[1], operator_event[0], 0);
            for (int i=0; i < 3; i++) {
                bc_unpack_cuda_gpu_async(bc, in2, operator_buf_gpu, i,
                                         recvreq, sendreq[i], ph + 2 * i,
                                         operator_stream[1], myblocks);
            }
            if (real) {
                bmgs_fd_cuda_gpu(&self->stencil_gpu, operator_buf_gpu, out2,
                                 boundary|GPAW_BOUNDARY_ONLY, myblocks,
                                 operator_stream[1]);
            } else {
                bmgs_fd_cuda_gpuz(&self->stencil_gpu,
                                  (const gpuDoubleComplex*) operator_buf_gpu,
                                  (gpuDoubleComplex*) out2,
                                  boundary|GPAW_BOUNDARY_ONLY, myblocks,
                                  operator_stream[1]);
            }
            gpuEventRecord(operator_event[1], operator_stream[1]);
        } else {
            bc_unpack_paste_cuda_gpu(bc, in2, operator_buf_gpu, recvreq,
                                     0, myblocks);
            for (int i=0; i < 3; i++) {
                bc_unpack_cuda_gpu(bc, in2, operator_buf_gpu, i,
                                   recvreq, sendreq[i], ph + 2 * i,
                                   0, myblocks);
            }
            if (real) {
                bmgs_fd_cuda_gpu(&self->stencil_gpu, operator_buf_gpu, out2,
                                 GPAW_BOUNDARY_NORMAL, myblocks, 0);
            } else {
                bmgs_fd_cuda_gpuz(&self->stencil_gpu,
                                  (const gpuDoubleComplex*) (operator_buf_gpu),
                                  (gpuDoubleComplex*) out2,
                                  GPAW_BOUNDARY_NORMAL, myblocks, 0);
            }
        }
    }

    if (cuda_overlap) {
        gpuStreamWaitEvent(0, operator_event[1], 0);
        gpuStreamSynchronize(operator_stream[0]);
    }
}

/*
 * Python interface for the GPU version of the FD algorithm
 * (similar to Operator_apply() for CPUs).
 *
 * arguments:
 *   input_gpu  -- pointer to device memory (GPUArray.gpudata)
 *   output_gpu -- pointer to device memory (GPUArray.gpudata)
 *   shape      -- shape of the array (tuple)
 *   type       -- datatype of array elements
 *   phases     -- phase (complex) (ignored if type is NPY_DOUBLE)
 */
PyObject * Operator_apply_cuda_gpu(OperatorObject* self, PyObject* args)
{
    PyArrayObject* phases = 0;
    void *input_gpu;
    void *output_gpu;
    PyObject *shape;
    PyArray_Descr *type;

    if (!PyArg_ParseTuple(args, "nnOO|O", &input_gpu, &output_gpu, &shape,
                          &type, &phases))
        return NULL;

    int nin = 1;
    if (PyTuple_Size(shape) == 4)
        nin = (int) PyLong_AsLong(PyTuple_GetItem(shape, 0));

    const double *in = (double*) input_gpu;
    double *out = (double*) output_gpu;

    bool real = (type->type_num == NPY_DOUBLE);

    const double_complex *ph;
    if (real)
        ph = 0;
    else
        ph = COMPLEXP(phases);

    boundary_conditions* bc = self->bc;
    int mpi_size = 1;
    if ((bc->maxsend || bc->maxrecv) && bc->comm != MPI_COMM_NULL) {
        MPI_Comm_size(bc->comm, &mpi_size);
    }
    int blocks = MAX(1, MIN(nin, MIN((GPU_BLOCKS_MIN) * mpi_size,
                                     (GPU_BLOCKS_MAX) / bc->ndouble)));

    if (gpaw_cuda_debug) {
        debug_operator_allocate(self, nin, blocks);
        debug_operator_memcpy_pre(in, out);
    }

    _operator_apply_cuda_gpu(self, in, out, nin, blocks, real, ph);

    if (gpaw_cuda_debug) {
        gpuDeviceSynchronize();
        debug_operator_memcpy_post(out, operator_buf_gpu);
        debug_operator_apply(self, nin, blocks, real, ph);
        debug_operator_deallocate();
    }

    if (PyErr_Occurred())
        return NULL;
    else
        Py_RETURN_NONE;
}
