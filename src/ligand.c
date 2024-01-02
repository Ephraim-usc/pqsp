#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "structmember.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

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
    {"n_particles", T_INT, offsetof(LigandObject, n_particles), READONLY, "first name"},
    {"n_sites", T_INT, offsetof(LigandObject, n_sites), READONLY, "last name"},
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
