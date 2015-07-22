#include "extensions.h"
//#include <stdlib.h>

PyObject *pc_potential(PyObject *self, PyObject *args)
{
    PyArrayObject* beg_v_obj;
    PyArrayObject* h_v_obj;
    PyArrayObject* q_p_obj;
    PyArrayObject* R_pv_obj;
    double rc;
    PyArrayObject* vext_G_obj;
    PyArrayObject* rhot_G_obj = 0;
    PyArrayObject* F_pv_obj = 0;
    if (!PyArg_ParseTuple(args, "OOOOdO|OO", &beg_v_obj, &h_v_obj, &q_p_obj,
                          &R_pv_obj, &rc, &vext_G_obj,
                          &rhot_G_obj, &F_pv_obj))
    return NULL;

    const long *beg_v = PyArray_DATA(beg_v_obj);
    const double *h_v = PyArray_DATA(h_v_obj);
    const double *q_p = PyArray_DATA(q_p_obj);
    const double *R_pv = PyArray_DATA(R_pv_obj);
    double *vext_G = PyArray_DATA(vext_G_obj);

    int np = PyArray_DIM(R_pv_obj, 0);
    npy_intp* n = PyArray_DIMS(vext_G_obj);

    const double* rhot_G = 0;
    double* F_pv = 0;
    double dV = 0.0;
    if (F_pv_obj != 0) {
        rhot_G = PyArray_DATA(rhot_G_obj);
        F_pv = PyArray_DATA(F_pv_obj);
        dV = h_v[0] * h_v[1] * h_v[2];
        }
    
    for (int i = 0; i < n[0]; i++) {
        double x = (beg_v[0] + i) * h_v[0];
        for (int j = 0; j < n[1]; j++) {
            double y = (beg_v[1] + j) * h_v[1];
            int ij = (i * n[1] + j) * n[2];
            for (int k = 0; k < n[2]; k++) {
                double z = (beg_v[2] + k) * h_v[2];
                for (int p = 0; p < np; p++) {
                    const double* R_v = R_pv + 3 * p;
                    double dx = R_v[0] - x;
                    double dy = R_v[1] - y;
                    double dz = R_v[2] - z;
                    double d = sqrt(dx * dx + dy * dy + dz * dz);
                    int G = ij + k;
                    if (F_pv == 0) {
                        double v;
                        if (d > rc)
                            v = q_p[p] / d;
                        else {
                            double x = d / rc;
                            double x2 = x * x;
                            v = q_p[p] * (3.28125 +
                                          x2 * (-5.46875 +
                                                x2 * (4.59375 +
                                                      x2 * -1.40625))) / rc;
                            }
                        vext_G[G] += v;
                        }
                    else {
                        double w;  // -(dv/dr)/r
                        if (d > rc)
                            w = 1 / (d * d * d);
                        else {
                            double x = d / rc;
                            double x2 = x * x;
                            w = (-2 * (-5.46875 +
                                       x2 * (2 * 4.59375 +
                                             x2 * 3 * -1.40625)) /
                                 (rc * rc * rc));
                            }
                        w *= q_p[p] * rhot_G[G] * dV;
                        double* F_v = F_pv + 3 * p;
                        F_v[0] += w * dx;
                        F_v[1] += w * dy;
                        F_v[2] += w * dz;
                        }
                    }
                }
            }
        }
    Py_RETURN_NONE;
}
