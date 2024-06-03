#ifndef MPA_H
#define MPA_H

#include <Python.h>

#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>
PyObject* evaluate_mpa_poly(PyObject *self, PyObject *args);

#include "extensions.h"

#endif
