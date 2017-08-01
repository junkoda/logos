#ifndef CORR_H
#define CORR_H 1

#include "Python.h"

PyMODINIT_FUNC
corr_module_init();

PyObject* py_corr_compute_ucorr(PyObject* self, PyObject* args);

#endif
