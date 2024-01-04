#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "structmember.h"
#include "expm.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>


/******************************************************************************
                                the Site type
******************************************************************************/

typedef struct {
    PyObject_HEAD
    int n;
    int __max_states__;
    int *targets;
    double **Qs;
    double **Ps;
} SiteObject;

static void
Site_dealloc(SiteObject *self)
{
    int i;
    free(self->targets);
    for (i = 0; i < self->__max_states__; i++) if (self->Qs[i]) free(self->Qs[i]);
    free(self->Qs);
    for (i = 0; i < self->__max_states__; i++) if (self->Ps[i]) free(self->Ps[i]);
    free(self->Ps);
    Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *
Site_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    SiteObject *self;
    self = (SiteObject *) type->tp_alloc(type, 0);
    if (self != NULL)
    {
        self->n = 0;
        self->__max_states__ = 1024;
        self->targets = NULL;
        self->Qs = NULL;
        self->Ps = NULL;
    }
    return (PyObject *) self;
}

static int
Site_init(SiteObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *targetsObj, *onsObj, *offsObj;
    int i;
    
    static char *kwlist[] = {"targets", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O!O!O!", kwlist, &PyList_Type, &targetsObj, &PyList_Type, &onsObj, &PyList_Type, &offsObj))
        return -1;
    
    self->n = (int) PyList_Size(targetsObj) + 1;
    self->targets = calloc(self->n, sizeof(int));
    for (i = 0; i < self->n; i++)
      self->targets[i] = (int) PyLong_AsLong(PyList_GetItem(targetsObj, i - 1));
    
    self->Qs = calloc(self->__max_states__, sizeof(double *));
    self->Qs[0] = calloc(self->n * self->n, sizeof(double));
    for (i = 1; i < self->n; i++)
      self->Qs[0][i] = (int) PyFloat_AsDouble(PyList_GetItem(onsObj, i - 1));
    for (i = 1; i < self->n; i++)
      self->Qs[0][self->n * i] = (int) PyFloat_AsDouble(PyList_GetItem(offsObj, i - 1));
    
    self->Ps = calloc(self->__max_states__, sizeof(double *));
    self->Ps[0] = r8mat_expm1(self->n, self->Qs[0]);
    return 0;
}

static PyMemberDef Site_members[] = {
    {"n", T_INT, offsetof(SiteObject, n), READONLY, "number of states"},
    {NULL}  /* Sentinel */
};

static PyObject *
Site_print(SiteObject *self, PyObject *Py_UNUSED(ignored))
{
    int i;
    for (i = 0; i < self->n; i++)
        printf("%d ", self->targets[i]);
    printf("\n");
    
    if (self->Ps[0] != NULL){
    for (i = 0; i < self->n * self->n; i++)
        printf("%f ", self->Ps[0][i]);
    printf("\n");
    }
    
    if (self->Qs[0] != NULL){
    for (i = 0; i < self->n * self->n; i++)
        printf("%f ", self->Qs[0][i]);
    printf("\n");
    }
    
    Py_RETURN_NONE;
}

static PyMethodDef Site_methods[] = {
    {"print", (PyCFunction) Site_print, METH_NOARGS, "print the Site matrix"},
    {NULL}  /* Sentinel */
};

static PyTypeObject SiteType = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "ligand.Site",
    .tp_doc = PyDoc_STR("Site objects"),
    .tp_basicsize = sizeof(SiteObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = Site_new,
    .tp_init = (initproc) Site_init,
    .tp_members = Site_members,
    .tp_methods = Site_methods,
    .tp_dealloc = (destructor) Site_dealloc,
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
    if (PyType_Ready(&SiteType) < 0)
        return NULL;
    if (PyType_Ready(&LigandType) < 0)
        return NULL;

    m = PyModule_Create(&ligandmodule);
    if (m == NULL)
        return NULL;
    
    Py_INCREF(&SiteType);
    if (PyModule_AddObject(m, "Site", (PyObject *) &SiteType) < 0) {
        Py_DECREF(&SiteType);
        Py_DECREF(m);
        return NULL;
    }
    
    Py_INCREF(&LigandType);
    if (PyModule_AddObject(m, "Ligand", (PyObject *) &LigandType) < 0) {
        Py_DECREF(&LigandType);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}
