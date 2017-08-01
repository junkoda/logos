#include "corr.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

PyMODINIT_FUNC
corr_module_init()
{
  import_array();

  return NULL;
}

void compute_ucorr(double const * const u, const int n, double* const corr)
{
  for(int r=0; r<n; ++r) {
    double uu= 0.0;
    for(int i=0; i<n; ++i) {
      int j= (i + r) % n;
      uu += u[i]*u[j];
    }
    corr[r]= uu/n;
  }
}

PyObject* py_corr_compute_ucorr(PyObject* self, PyObject* args)
{
  PyObject *py_u, *py_corr;

  if(!PyArg_ParseTuple(args, "OO", &py_u, &py_corr)) {
    return NULL;
  }

  Py_buffer buf_u;
  if(PyObject_GetBuffer(py_u, &buf_u,
			PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
    return NULL;

  Py_buffer buf_corr;
  
  if(PyObject_GetBuffer(py_u, &buf_corr,
			PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1)
    return NULL;

  if(buf_u.ndim != 1 || buf_corr.ndim != 1) {
    PyErr_SetString(PyExc_TypeError, "Expected a 1-dimensional array");
    PyBuffer_Release(&buf_u);
    PyBuffer_Release(&buf_corr);
    return NULL;
  }

  if(strcmp(buf_u.format,"d") != 0 || strcmp(buf_corr.format,"d") != 0) {
    PyErr_SetString(PyExc_TypeError, "Expected an array of doubles");
    PyBuffer_Release(&buf_u);
    PyBuffer_Release(&buf_corr);
    return NULL;
  }

  if(buf_u.shape[0] != buf_corr.shape[0]) {
    PyErr_SetString(PyExc_TypeError, "Expected arrays of same length");
    PyBuffer_Release(&buf_u);
    PyBuffer_Release(&buf_corr);
    return NULL;
  }

  compute_ucorr((double*) buf_u.buf, buf_u.shape[0],
		(double*) buf_corr.buf);

  PyBuffer_Release(&buf_u);
  PyBuffer_Release(&buf_corr);

  Py_RETURN_NONE;
}

