#include "Python.h"

#include "corr.h"

static PyMethodDef methods[] = {
  {"_corr_compute_ucorr", py_corr_compute_ucorr, METH_VARARGS,
   "compute velocity correlation function"},
  {"_corr_compute_umoments", py_corr_compute_umoments, METH_VARARGS,
   "compute relative velocity moments"},

    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef module = {
  PyModuleDef_HEAD_INIT,
  "_logos", // name of this module
  "A package for 1D redshift distortion", // Doc String
  -1,
  methods
};

PyMODINIT_FUNC
PyInit__logos(void) {
  corr_module_init();    
  
  return PyModule_Create(&module);
}

