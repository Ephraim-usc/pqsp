#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#include "numpy/arrayobject.h"

typedef struct {
    PyObject_HEAD
    int n_particles;
    int n_sites;
    //int *compartment;
    //int *states;
    int *bindings;
} LigandObject;

static void
Ligand_dealloc(LigandObject *self)
{
    //Py_XDECREF(self->first);
    free(self->bindings);
    Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *
Ligand_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    LigandObject *self;
    self = (LigandObject *) type->tp_alloc(type, 0);
    if (self != NULL) {
        self->n_particles = 0;
        self->n_sites = 0;
        self->bindings = NULL;
    }
    return (PyObject *) self;
}

static int
Ligand_init(LigandObject *self, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"n_particles", "n_sites", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|ii", &self->n_particles, &self->n_sites))
        return -1;
    self->bindings = calloc(self->n_particles * self->n_sites, sizeof(int));
    return 0;
}

static PyMemberDef Ligand_members[] = {
    {"n_particles", Py_T_INT, offsetof(LigandObject, first), Py_READONLY, "first name"},
    {"n_sites", Py_T_INT, offsetof(LigandObject, last), Py_READONLY, "last name"},
    {NULL}  /* Sentinel */
};

static PyTypeObject LigandType = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "ligand.Ligand",
    .tp_doc = PyDoc_STR("Ligand objects"),
    .tp_basicsize = sizeof(LigandObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = Ligand_new,
    .tp_member = Ligand_member,
    .tp_dealloc = (destructor) Ligand_dealloc,
};

static PyModuleDef ligandmodule = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "ligand",
    .m_doc = "Example module that creates an extension type.",
    .m_size = -1,
};

PyMODINIT_FUNC
PyInit_ligand(void)
{
    PyObject *m;
    if (PyType_Ready(&LigandType) < 0)
        return NULL;

    m = PyModule_Create(&ligandmodule);
    if (m == NULL)
        return NULL;

    Py_INCREF(&LigandType);
    if (PyModule_AddObject(m, "Ligand", (PyObject *) &LigandType) < 0) {
        Py_DECREF(&LigandType);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}












/*

typedef double DTYPE;
typedef unsigned long ITYPE;

typedef struct matrix {
  ITYPE n;
  DTYPE *data;
} matrix;

matrix* new_matrix(ITYPE n)
{
  matrix* mat;
  mat = (matrix *)malloc(sizeof(matrix));
  mat->n = n;
  mat->data = calloc(n*n,sizeof(DTYPE));
  return mat;
}

void destroy_matrix(matrix* mat)
{
  free(mat->data);
  free(mat);
}

void print_matrix(matrix* mat)
{
  ITYPE n = mat->n;
  printf("matrix of dimension %lu\n", n);
  
  ITYPE i, j;
  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++)
      printf("%lf ", mat->data[i*n + j]);
    printf("\n");
  }
}

void add_square(matrix* mat, ITYPE* idx, ITYPE len_idx, DTYPE q)
{
  ITYPE* x = idx;
  ITYPE* y = idx;
  ITYPE* end = idx + len_idx;
  DTYPE* data = mat->data;
  ITYPE n = mat->n;
  
  for (; x < end; x++)
    for (y=idx; y <= x; y++)
    {
      data[*x * n + *y] += q;
    }
}

void set_square(matrix* mat, ITYPE* idx, ITYPE len_idx, DTYPE q)
{
  ITYPE* x = idx;
  ITYPE* y = idx;
  ITYPE* end = idx + len_idx;
  DTYPE* data = mat->data;
  ITYPE n = mat->n;
  
  for (; x < end; x++)
    for (y=idx; y <= x; y++)
    {
      data[*x * n + *y] = q;
    }
}

void add(matrix* mat_1, matrix* mat_2)
{
  ITYPE n = mat_1->n;
  DTYPE* p = mat_1->data;
  DTYPE* q = mat_2->data;
  ITYPE i;
  
  for (i = 0; i < n*n; i++)
  {
    *(p+i) += *(q+i);
  }
}

void set_zeros(matrix* mat_1)
{
  ITYPE n = mat_1->n;
  DTYPE* data = mat_1->data;
  
  memset(data, 0, n*n * sizeof(DTYPE));
}

static void parse_py_int_seq(PyObject *py_int_seq, ITYPE** pr, ITYPE* len)
{
  *len = (ITYPE)PySequence_Fast_GET_SIZE(py_int_seq);
  *pr = (ITYPE *) malloc(sizeof(ITYPE) * *len);
  ITYPE i;
  for (i = 0; i < *len; i++) 
  {
    PyObject *item = PySequence_Fast_GET_ITEM(py_int_seq, i);
    (*pr)[i] = (ITYPE)PyLong_AsLong(item);
  }
}

static PyObject* py_new_matrix(PyObject* self, PyObject* args)
{ 
  ITYPE n;
  PyArg_ParseTuple(args, "k", &n);
  
  matrix* mat = new_matrix((ITYPE)n);
  PyObject* py_mat = PyCapsule_New((void *)mat, "matrix._matrix_C_API", NULL);
  
  return py_mat;
}

static PyObject* py_destroy_matrix(PyObject* self, PyObject* args)
{
  PyObject* py_mat;
  PyArg_UnpackTuple(args, NULL, 1, 1, &py_mat);
  
  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");
  destroy_matrix(mat);
  
  Py_RETURN_NONE;
}

static PyObject* py_add_square(PyObject* self, PyObject* args)
{
  PyObject* py_mat;
  PyObject* py_idx;
  PyObject* py_q;
  PyArg_UnpackTuple(args, NULL, 3, 3, &py_mat, &py_idx, &py_q);
  
  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");
  DTYPE q = (DTYPE)PyFloat_AS_DOUBLE(py_q);
  
  PyObject* py_int_seq;
  py_int_seq = PySequence_Fast(py_idx, NULL);
  
  ITYPE* idx; ITYPE len_idx;
  parse_py_int_seq(py_int_seq, &idx, &len_idx);
  add_square(mat, idx, len_idx, q);
  free(idx);
  
  Py_DECREF(py_int_seq);
  Py_RETURN_NONE;
}

static PyObject* py_set_square(PyObject* self, PyObject* args)
{
  PyObject* py_mat;
  PyObject* py_idx;
  PyObject* py_q;
  PyArg_UnpackTuple(args, NULL, 3, 3, &py_mat, &py_idx, &py_q);
  
  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");
  DTYPE q = (DTYPE)PyFloat_AS_DOUBLE(py_q);
  
  PyObject* py_int_seq;
  py_int_seq = PySequence_Fast(py_idx, NULL);
  
  ITYPE* idx; ITYPE len_idx;
  parse_py_int_seq(py_int_seq, &idx, &len_idx);
  set_square(mat, idx, len_idx, q);
  free(idx);
  
  Py_DECREF(py_int_seq);
  Py_RETURN_NONE;
}

static PyObject* py_add(PyObject* self, PyObject* args)
{
  PyObject* py_mat_1;
  PyObject* py_mat_2;
  PyArg_UnpackTuple(args, NULL, 2, 2, &py_mat_1, &py_mat_2);
  
  matrix* mat_1 = (matrix *)PyCapsule_GetPointer(py_mat_1, "matrix._matrix_C_API");
  matrix* mat_2 = (matrix *)PyCapsule_GetPointer(py_mat_2, "matrix._matrix_C_API");
  
  add(mat_1, mat_2);
  
  Py_RETURN_NONE;
}

static PyObject* py_set_zeros(PyObject* self, PyObject* args)
{
  PyObject* py_mat;
  PyArg_UnpackTuple(args, NULL, 1, 1, &py_mat);
  
  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");
  set_zeros(mat);
  
  Py_RETURN_NONE;
}

static PyObject* py_print_matrix(PyObject* self, PyObject* args)
{
  PyObject* py_mat; 
  PyArg_UnpackTuple(args, NULL, 1, 1, &py_mat);
  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");
  print_matrix(mat);
  
  Py_RETURN_NONE;
}

static PyObject* py_export_list(PyObject* self, PyObject* args) // export to python list, elements are copied. To release memory, you have to destroy the mat_C.
{
  PyObject* py_mat; 
  PyArg_UnpackTuple(args, NULL, 1, 1, &py_mat);
  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");
  DTYPE* data = mat->data; ITYPE n = mat->n;
  
  PyObject *py_list = Py_BuildValue("[]");
  ITYPE i = 0;
  for (; i < n*n; i++)
  {
    PyList_Append(py_list, Py_BuildValue("d", data[i]));
  }
  
  return py_list;
}

static PyObject* py_export_ndarray(PyObject* self, PyObject* args) // export to numpy ndarray, elements are used in place. Do not destroy the mat_C.
{
  PyObject* py_mat; 
  PyArg_UnpackTuple(args, NULL, 1, 1, &py_mat);
  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");
  DTYPE* data = mat->data; ITYPE n = mat->n;
  
  import_array();
  ITYPE dims[2];
  dims[0] = dims[1] = n;
  PyObject *py_ndarray = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, data);
  
  return py_ndarray;
}

static PyMethodDef myMethods[] = 
{
  {"new_matrix", py_new_matrix, METH_VARARGS, "new matrix"},
  {"print_matrix", py_print_matrix, METH_VARARGS, "print matrix"},
  {"add_square", py_add_square, METH_VARARGS, "add_square"},
  {"set_square", py_set_square, METH_VARARGS, "set_square"},
  {"add", py_add, METH_VARARGS, "add"},
  {"set_zeros", py_set_zeros, METH_VARARGS, "set_zeros"},
  {"destroy_matrix", py_destroy_matrix, METH_VARARGS, "destroy matrix"},
  {"export_list", py_export_list, METH_VARARGS, "export as list"},
  {"export_ndarray", py_export_ndarray, METH_VARARGS, "export as ndarray"},
  {NULL, NULL, 0, NULL},
};

static struct PyModuleDef ligandModule =
{
  PyModuleDef_HEAD_INIT,
  "ligandModule",
  "ligand Module",
  -1,
  myMethods
};

PyMODINIT_FUNC PyInit_ligand(void)
{
  return PyModule_Create(&ligandModule);
}

*/
