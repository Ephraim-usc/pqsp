#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "structmember.h"
#include "expm.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>


/******************************************************************************
                                the Transition type
******************************************************************************/

typedef struct {
    PyObject_HEAD
    int n;
    int *targets;
    double *Q;
    double *P;
} TransitionObject;

static void
Transition_dealloc(TransitionObject *self)
{
    free(self->targets);
    free(self->Q);
    free(self->P);
    Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *
Transition_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    TransitionObject *self;
    self = (TransitionObject *) type->tp_alloc(type, 0);
    if (self != NULL)
    {
        self->n = 0;
        self->targets = NULL;
        self->P = NULL;
        self->Q = NULL;
    }
    return (PyObject *) self;
}

static int
Transition_init(TransitionObject *self, PyObject *args, PyObject *kwds)
{
    PyObject * targetsObj;
    int i;
    
    static char *kwlist[] = {"targets", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O!", kwlist, &PyList_Type, &targetsObj))
        return -1;
    
    self->n = PyList_Size(targetsObj);
    for (i = 0; i < self->n; i++)
      self->targets[i] = i //PyFloat_AsDouble(PyList_GetItem(targetsObj, i));
    self->P = calloc(self->n, sizeof(*double));
    self->Q = r8mat_expm1(self->n, self->P);
    return 0;
}

static PyMemberDef Transition_members[] = {
    {"n", T_INT, offsetof(TransitionObject, n), READONLY, "number of states"},
    {NULL}  /* Sentinel */
};

static PyObject *
Transition_print(TransitionObject *self, PyObject *Py_UNUSED(ignored))
{
    int i;
    for (i = 0; i < self->n; i++) 
        printf("%d ", self->targets[i]);
    printf("\n");
    for (i = 0; i < self->n * self->n; i++)
        printf("%f ", self->P[i]);
    printf("\n");
    for (i = 0; i < self->n * self->n; i++)
        printf("%f ", self->Q[i]);
    printf("\n");
    
    Py_RETURN_NONE;
}

static PyMethodDef Transition_methods[] = {
    {"print", (PyCFunction) Transition_print, METH_NOARGS, "print the transition matrix"},
    {NULL}  /* Sentinel */
};

static PyTypeObject TransitionType = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "ligand.Transition",
    .tp_doc = PyDoc_STR("Transition objects"),
    .tp_basicsize = sizeof(TransitionObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = Transition_new,
    .tp_init = (initproc) Transition_init,
    .tp_members = Transition_members,
    .tp_methods = Transition_methods,
    .tp_dealloc = (destructor) Transition_dealloc,
};



/******************************************************************************
                                the Ligand type
******************************************************************************/

typedef struct {
    PyObject_HEAD
    long n_particles;
    int n_sites;
    //int *compartment;
    //int *states;
    int *bindings;
} LigandObject;

static void
Ligand_dealloc(LigandObject *self)
{
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
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|li", kwlist, &self->n_particles, &self->n_sites))
        return -1;
    self->bindings = calloc(self->n_particles * self->n_sites, sizeof(int));
    return 0;
}

static PyMemberDef Ligand_members[] = {
    {"n_particles", T_INT, offsetof(LigandObject, n_particles), READONLY, "number of particles"},
    {"n_sites", T_INT, offsetof(LigandObject, n_sites), READONLY, "number of binding sites"},
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
    .tp_init = (initproc) Ligand_init,
    .tp_members = Ligand_members,
    .tp_dealloc = (destructor) Ligand_dealloc,
};


/******************************************************************************
                                the module
******************************************************************************/

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
    if (PyType_Ready(&TransitionType) < 0)
        return NULL;
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
