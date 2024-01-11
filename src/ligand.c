#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "structmember.h"
#include "expm.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>


static int *
PyList2Array_INT(PyObject *listObj)
{
    int len, i;
    len = (int) PyList_Size(listObj);
    int *array = calloc(len, sizeof(int));
    for (i = 0; i < len; i++)
      array[i] = (int) PyLong_AsLong(PyList_GetItem(listObj, i));
    return array;
}

static PyObject *
Array2PyList_INT(int *array, int len)
{
    int i;
    PyObject *listObj;
    listObj = PyList_New(len);
    for (i = 0; i < len; i++)
        PyList_SetItem(listObj, i, (PyObject *) Py_BuildValue("i", array[i]));
    return listObj;
}



/******************************************************************************
                                the Transition type
******************************************************************************/

typedef struct {
    PyObject_HEAD
    
    /* fixed attributes */
    int n_compartments; // number of compartments
    int n_states; // number of analytes
    int n_targets; // number of targets
    int *targets;
    
    /* variable attributes */
    double ***Pses; // list of P matrices, for each state, for each compartment
} TransitionObject;

static void
Transition_dealloc(TransitionObject *self)
{
    int c, s;
    free(self->targets);
    for (c = 0; c < self->n_compartments; c++)
    {
        for (s = 0; s < self->n_states; s++)
            if (self->Pses[c][s])
                free(self->Pses[c][s]);
        free(self->Pses[c]);
    }
    Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *
Transition_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    TransitionObject *self;
    self = (TransitionObject *) type->tp_alloc(type, 0);
    if (self != NULL)
    {
        self->n_compartments = 0;
        self->n_states = 0;
        self->n_targets = 0;
        self->targets = NULL;
        self->Pses = NULL;
    }
    return (PyObject *) self;
}

static int
Transition_init(TransitionObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *targetsObj;
    int c, i;
    
    static char *kwlist[] = {"n_compartments", "n_states", "targets", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|iiO!", kwlist, &self->n_compartments, &self->n_states, &PyList_Type, &targetsObj))
        return -1;
    
    self->n_targets = (int) PyList_Size(targetsObj);
    self->targets = calloc(self->n_targets, sizeof(int));
    for (i = 0; i < self->n_targets; i++)
        self->targets[i] = (int) PyLong_AsLong(PyList_GetItem(targetsObj, i));
    
    self->Pses = calloc(self->n_compartments, sizeof(double **));
    for (c = 0; c < self->n_compartments; c++)
        self->Pses[c] = calloc(self->n_states, sizeof(double *));
    
    return 0;
}

static PyMemberDef Transition_members[] = {
    {"n_compartments", T_INT, offsetof(TransitionObject, n_compartments), READONLY, "number of compartments"},
    {"n_states", T_INT, offsetof(TransitionObject, n_states), READONLY, "number of states"},
    {"n_targets", T_INT, offsetof(TransitionObject, n_targets), READONLY, "number of targets"},
    {NULL}  /* Sentinel */
};


static PyObject *
Transition_print(TransitionObject *self, PyObject *Py_UNUSED(ignored))
{
    int c, s, i;
    
    printf("targets:\n");
    for (i = 0; i < self->n_targets; i++)
        printf("%d ", self->targets[i]);
    printf("\n\n");
    
    printf("P matrices:\n");
    for (c = 0; c < self->n_compartments; c++)
    {
        for (s = 0; s < self->n_states; s++)
        if (self->Pses[c][s])
        {
            printf("[compartment %d, state %d] ", c, s);
            for (i = 0; i < (self->n_targets + 1) * (self->n_targets + 1); i++)
                printf("%f ", self->Pses[c][s][i]);
            printf("\n");
        }
    }
    printf("\n");
    
    Py_RETURN_NONE;
}

static PyObject *
Transition_set_P(TransitionObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *PObj;
    int c, s, i;
    
    static char *kwlist[] = {"compartment", "state", "P", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|iiO!", kwlist, &c, &s, &PyList_Type, &PObj))
        Py_RETURN_NONE;
    
    if (!self->Pses[c][s])
        self->Pses[c][s] = calloc(self->n_targets + 1, sizeof(double));
    for (i = 0; i < (self->n_targets + 1) * (self->n_targets + 1); i++)
        self->Pses[c][s][i] = (double) PyFloat_AsDouble(PyList_GetItem(PObj, i));
    
    Py_RETURN_NONE;
}

static int 
_Transition_apply(TransitionObject *self, int n_particles, int *compartments, int *states, int *values, int *deltas)
{
    int p, x, x_;
    double *P;
    double tmp;
    
    for (p = 0; p < n_particles; p++)
    {
        x = values[p];
        if (self->Pses[compartments[p]][states[p]])
            P = self->Pses[compartments[p]][states[p]] + x * (self->n_targets + 1); // transition matrix + shift for starting state x = transition vector for x
        else
            P = self->Pses[compartments[p]][0] + x * (self->n_targets + 1);
        
        tmp = drand48();
        x_ = 0;
        for (x_ = 0; tmp -= P[x_], tmp > 0; x_++);

        deltas[x] -= 0;
        deltas[x_] += 0;
        values[p] = x_;
    }
    
    return 0;
}

static PyObject *
Transition_apply(TransitionObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *compartmentsObj, *statesObj, *valuesObj, *deltasObj;
    int n_particles;
    int *compartments, *states, *values, *deltas;
    
    static char *kwlist[] = {"compartments", "states", "values", "deltas", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O!O!O!O!", kwlist, &PyList_Type, &compartmentsObj, &PyList_Type, &statesObj, &PyList_Type, &valuesObj, &PyList_Type, &deltasObj))
        Py_RETURN_NONE;
    
    n_particles = (int) PyList_Size(valuesObj);
    
    compartments = PyList2Array_INT(compartmentsObj);
    states = PyList2Array_INT(statesObj);
    values = PyList2Array_INT(valuesObj);
    deltas = PyList2Array_INT(deltasObj);
    
    _Transition_apply(self, n_particles, compartments, states, values, deltas);
    
    valuesObj = Array2PyList_INT(values, n_particles);
    deltasObj = Array2PyList_INT(deltas, n_particles);
    
    free(compartments);
    free(states);
    free(values);
    free(deltas);
    
    return Py_NewRef(Py_BuildValue("(OO)", valuesObj, deltasObj));
}

static PyMethodDef Transition_methods[] = {
    {"print", (PyCFunction) Transition_print, METH_NOARGS, "print"},
    {"set_P", (PyCFunction) Transition_set_P, METH_VARARGS | METH_KEYWORDS, "set P matrix for a given compartment for a given state"},
    {"apply", (PyCFunction) Transition_apply, METH_VARARGS | METH_KEYWORDS, "apply Transition to particles"},
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
                                the Dict type
******************************************************************************/


typedef struct Node {
    struct Node **children;
    int value;
} Node;

typedef struct {
    PyObject_HEAD
    
    /* fixed attributes */
    int n; // number of possible tokens
    Node *root;
    
    /* variable attributes */
} DictObject;



/******************************************************************************
                                the Site type
******************************************************************************/

typedef struct {
    PyObject_HEAD

    /* fixed attributes */
    int n; // number of targets
    int __max_states__;
    int *targets; // list of targets
    double *offs; // list of off rates
    double **onses; // list of lists of on rates, for each state
    
    /* variable attributes */
    double **Ps; // list of P matrices (each a (n+1)x(n+1) matrix stored in linear form), for each state, for each compartment
} SiteObject;

static void
Site_dealloc(SiteObject *self)
{
    int i;
    free(self->targets);
    free(self->offs);
    for (i = 0; i < self->__max_states__; i++) if (self->onses[i]) free(self->onses[i]);
    free(self->onses);
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
        self->offs = NULL;
        self->onses = NULL;
        self->Ps = NULL;
    }
    return (PyObject *) self;
}

static int
Site_init(SiteObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *targetsObj, *onsObj, *offsObj;
    int i;
    
    static char *kwlist[] = {"targets", "ons", "offs", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O!O!O!", kwlist, &PyList_Type, &targetsObj, &PyList_Type, &onsObj, &PyList_Type, &offsObj))
        return -1;
    
    self->n = (int) PyList_Size(targetsObj);
    
    self->targets = calloc(self->n, sizeof(int));
    for (i = 0; i < self->n; i++)
        self->targets[i] = (int) PyLong_AsLong(PyList_GetItem(targetsObj, i));

    self->onses = calloc(self->__max_states__, sizeof(double *));
    self->onses[0] = calloc(self->n, sizeof(double));
    for (i = 0; i < self->n; i++)
        self->onses[0][i] = (double) PyFloat_AsDouble(PyList_GetItem(onsObj, i));
    
    self->offs = calloc(self->n, sizeof(double));
    for (i = 0; i < self->n; i++)
        self->offs[i] = (double) PyFloat_AsDouble(PyList_GetItem(offsObj, i));
    
    self->Ps = calloc(self->__max_states__, sizeof(double *));
    //self->Ps[0] = r8mat_expm1(self->n, self->Qs[0]);
    return 0;
}

static PyMemberDef Site_members[] = {
    {"n", T_INT, offsetof(SiteObject, n), READONLY, "number of states"},
    {NULL}  /* Sentinel */
};

static PyObject *
Site_print(SiteObject *self, PyObject *Py_UNUSED(ignored))
{
    int state, i;
    
    printf("targets:\n");
    for (i = 0; i < self->n; i++)
        printf("%d ", self->targets[i]);
    printf("\n\n");
    
    printf("off rates:\n");
    for (i = 0; i < self->n; i++)
        printf("%f ", self->offs[i]);
    printf("\n\n");
    
    printf("on rates:\n");
    for (state = 0; state < self->__max_states__; state++)
        if (self->onses[state] != NULL)
        {
            printf("[state %d] ", state);
            for (i = 0; i < self->n; i++)
                printf("%f ", self->onses[state][i]);
            printf("\n");
        }
    printf("\n");
    
    printf("Ps:\n");
    for (state = 0; state < self->__max_states__; state++)
        if (self->Ps[state] != NULL)
        {
            printf("[state %d] ", state);
            for (i = 0; i < (self->n + 1) * (self->n + 1); i++)
                printf("%f ", self->Ps[state][i]);
            printf("\n");
        }
    printf("\n");
    
    Py_RETURN_NONE;
}

static PyObject *
Site_add_state(SiteObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *onsObj;
    int state, i;
    
    static char *kwlist[] = {"state", "ons", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|iO!", kwlist, &state, &PyList_Type, &onsObj))
        Py_RETURN_NONE;
    
    if (!self->onses[state]) // if already added, then replace
      self->onses[state] = calloc(self->n, sizeof(double));
    for (i = 0; i < self->n; i++)
        self->onses[state][i] = (double) PyFloat_AsDouble(PyList_GetItem(onsObj, i));
    
    Py_RETURN_NONE;
}


static PyObject *
Site_compute_Ps(SiteObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *xsObj;
    int state, i;
    double t;
    double *xs = calloc(self->n, sizeof(double));
    double *Q = calloc((self->n + 1) * (self->n + 1), sizeof(double));
    
    static char *kwlist[] = {"t", "xs", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dO!", kwlist, &t, &PyList_Type, &xsObj))
        Py_RETURN_NONE;
    
    for (i = 0; i < self->n; i++)
        xs[i] = (double) PyFloat_AsDouble(PyList_GetItem(xsObj, i));
    
    for (state = 0; state < self->__max_states__; state++)
        if (self->onses[state] != NULL)
        {
            Q[0] = 0.0;
            for (i = 0; i < self->n; i++)
            {
                Q[0] -= self->onses[state][i] * xs[i] * t;
                Q[i + 1] = self->onses[state][i] * xs[i] * t;
            }
            for (i = 0; i < self->n; i++)
            {
                Q[(self->n + 1) * (i + 1)] = self->offs[i] * t;
                Q[(self->n + 2) * (i + 1)] = - self->offs[i] * t;
            }
            /* printf("[state %d] ", state);
            for (i = 0; i < (self->n + 1) * (self->n + 1); i++)
                printf("%f ", Q[i]);
            printf("\n"); */
            self->Ps[state] = r8mat_expm1(self->n + 1, Q);
        }
    
    Py_RETURN_NONE;
}

// mpp: mole per particle
// bounds: list of bound targets for each particle
// xs: list of concentrations for each target
// xs_: list of updated concentrations for each target
static int
Site_bind(SiteObject *self, int n_particles, double mpp, double *bounds, double *xs, double *xs_) 
{
    ;
    return 0;
}

static PyMethodDef Site_methods[] = {
    {"print", (PyCFunction) Site_print, METH_NOARGS, "print the Site"},
    {"add_state", (PyCFunction) Site_add_state, METH_VARARGS | METH_KEYWORDS, "add a state to the Site"},
    {"compute_Ps", (PyCFunction) Site_compute_Ps,  METH_VARARGS | METH_KEYWORDS, "compute P matrices"},
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
    SiteObject **sites;
    int n_sites;
    long n_particles;
    double mpp; // mole per particle
    int *states; // list of states for each particle
    int **boundses; // list of lists of bound targets for each binding site
} LigandObject;

static void
Ligand_dealloc(LigandObject *self)
{
    int s;
    for (s = 0; s < self->n_sites; s++)
      Py_XDECREF(self->sites[s]);
      free(self->boundses[s]);
    free(self->states);
    free(self->boundses);
    Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *
Ligand_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    LigandObject *self;
    self = (LigandObject *) type->tp_alloc(type, 0);
    if (self != NULL) {
        self->n_sites = 0;
        self->sites = NULL;
        self->n_particles = 0;
        self->mpp = 0.0;
        self->states = NULL;
        self->boundses = NULL;
    }
    return (PyObject *) self;
}

static int
Ligand_init(LigandObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *sitesObj;
    int s;
    
    static char *kwlist[] = {"sites", "n_particles", "mpp", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O!ld", kwlist, &PyList_Type, &sitesObj, &self->n_particles, &self->mpp))
        return -1;
    
    self->n_sites = (int) PyList_Size(sitesObj);
    self->sites = calloc(self->n_sites, sizeof(SiteObject *));
    for(s = 0; s < self->n_sites; s++)
      self->sites[s] = (SiteObject *) PyList_GetItem(sitesObj, s);
    self->states = calloc(self->n_particles, sizeof(int));
    self->boundses = calloc(self->n_sites, sizeof(int *));
    for(s = 0; s < self->n_sites; s++)
      self->boundses[s] = calloc(self->n_particles, sizeof(int));
    return 0;
}

static PyMemberDef Ligand_members[] = {
    {"n_sites", T_INT, offsetof(LigandObject, n_sites), READONLY, "number of binding sites"},
    {"n_particles", T_INT, offsetof(LigandObject, n_particles), READONLY, "number of particles"},
    {"mpp", T_DOUBLE, offsetof(LigandObject, mpp), READONLY, "mole per particle"},
    {NULL}  /* Sentinel */
};

static PyObject *
Ligand_getsites(LigandObject *self, void *closure)
{
    int s;
    PyObject *sitesObj = PyList_New(self->n_sites);
    for (s = 0; s < self->n_sites; s++)
        PyList_SetItem(sitesObj, s, (PyObject *) self->sites[s]);
    return Py_NewRef(sitesObj);
}

static PyObject *
Ligand_getstates(LigandObject *self, void *closure)
{
    int i;
    PyObject *statesObj = PyList_New(self->n_particles);
    for (i = 0; i < self->n_particles; i++)
        PyList_SetItem(statesObj, i, (PyObject *) Py_BuildValue("i", self->states[i]));
    return Py_NewRef(statesObj);
}

static PyObject *
Ligand_getboundses(LigandObject *self, void *closure)
{
    int s, i;
    PyObject *boundsesObj = PyList_New(self->n_sites);
    PyObject *boundsObj;
    for (s = 0; s < self->n_sites; s++)
    {
        boundsObj = PyList_New(self->n_particles);
        for (i = 0; i < self->n_particles; i++)
            PyList_SetItem(boundsObj, i, (PyObject *) Py_BuildValue("i", self->boundses[s][i]));
        PyList_SetItem(boundsesObj, s, (PyObject *) boundsObj);
    }
    return Py_NewRef(boundsesObj);
}

static PyGetSetDef Ligand_getsetters[] = {
    {"sites", (getter) Ligand_getsites, NULL, "binding sites", NULL},
    {"states", (getter) Ligand_getstates, NULL, "states of each particle", NULL},
    {"boundses", (getter) Ligand_getboundses, NULL, "bound targets on each site for each particle", NULL},
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
    .tp_getset = Ligand_getsetters,
};


/******************************************************************************
                                the System type
******************************************************************************/

typedef struct {
    PyObject_HEAD

    /* fixed attributes */
    int n_compartments; // number of compartments
    int n_analytes; // number of analytes
    int __max_ligands__;
    
    /* variable attributes */
    double **xses; // list of lists of analyte concentrations, for each compartment
    int n_ligands;
    LigandObject **ligands;
} SystemObject;

static void
System_dealloc(SystemObject *self)
{
    int c;
    free(self->ligands);
    for (c = 0; c < self->n_compartments; c++) free(self->xses[c]);
    Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *
System_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    SystemObject *self;
    self = (SystemObject *) type->tp_alloc(type, 0);
    if (self != NULL)
    {
        self->n_compartments = 0;
        self->n_analytes = 0;
        self->xses = NULL;
        self->__max_ligands__ = 1024;
        self->n_ligands = 0;
        self->ligands = NULL;
    }
    return (PyObject *) self;
}

static int
System_init(SystemObject *self, PyObject *args, PyObject *kwds)
{
    int c;
    
    static char *kwlist[] = {"n_compartments", "n_analytes", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|ii", kwlist, &self->n_compartments, &self->n_analytes))
        return -1;
    
    self->xses = calloc(self->n_compartments, sizeof(double *));
    for (c = 0; c < self->n_compartments; c++)
        self->xses[c] = calloc(self->n_analytes, sizeof(double));
    
    self->n_ligands = 0;
    self->ligands = calloc(self->__max_ligands__, sizeof(LigandObject *));
    return 0;
}

static PyMemberDef System_members[] = {
    {"n_compartments", T_INT, offsetof(SystemObject, n_compartments), READONLY, "number of compartments"},
    {"n_analytes", T_INT, offsetof(SystemObject, n_analytes), READONLY, "number of analytes"},
    {"n_ligands", T_INT, offsetof(SystemObject, n_ligands), READONLY, "number of ligands"},
    {NULL}  /* Sentinel */
};


static PyObject *
System_getxses(SystemObject *self, void *closure)
{
    int c, a;
    PyObject *xsesObj = PyList_New(self->n_compartments);
    PyObject *xsObj;
    for (c = 0; c < self->n_compartments; c++)
    {
        xsObj = PyList_New(self->n_analytes);
        for (a = 0; a < self->n_analytes; a++)
            PyList_SetItem(xsObj, a, (PyObject *) Py_BuildValue("d", self->xses[c][a]));
        PyList_SetItem(xsesObj, c, (PyObject *) xsObj);
    }
    return Py_NewRef(xsesObj);
}

static PyObject *
System_getligands(SystemObject *self, void *closure)
{
    int l;
    PyObject *ligandsObj = PyList_New(self->n_ligands);
    for (l = 0; l < self->n_ligands; l++)
        PyList_SetItem(ligandsObj, l, (PyObject *) self->ligands[l]);
    return Py_NewRef(ligandsObj);
}

static PyGetSetDef System_getsetters[] = {
    {"xses", (getter) System_getxses, NULL, "analyte concentrations", NULL},
    {"ligands", (getter) System_getligands, NULL, "ligands", NULL},
    {NULL}  /* Sentinel */
};


static PyObject *
System_add_x(SystemObject *self, PyObject *args, PyObject *kwds)
{
    int c, a;
    double x;
    
    static char *kwlist[] = {"compartment", "analyte", "value", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|iid", kwlist, &c, &a, &x))
        Py_RETURN_NONE;
    
    self->xses[c][a] += x;
    Py_RETURN_NONE;
}

static PyObject *
System_add_ligand(SystemObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *ligandObj;
    
    static char *kwlist[] = {"ligand", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O!", kwlist, &LigandType, &ligandObj))
        Py_RETURN_NONE;
    
    self->ligands[self->n_ligands] = (LigandObject *) ligandObj;
    self->n_ligands += 1;
    Py_RETURN_NONE;
}


static PyMethodDef System_methods[] = {
    {"add_x", (PyCFunction) System_add_x, METH_VARARGS | METH_KEYWORDS, "add analyte"},
    {"add_ligand", (PyCFunction) System_add_ligand, METH_VARARGS | METH_KEYWORDS, "add ligand"},
    {NULL}  /* Sentinel */
};

static PyTypeObject SystemType = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "ligand.System",
    .tp_doc = PyDoc_STR("System objects"),
    .tp_basicsize = sizeof(SystemObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = System_new,
    .tp_init = (initproc) System_init,
    .tp_members = System_members,
    .tp_methods = System_methods,
    .tp_dealloc = (destructor) System_dealloc,
    .tp_getset = System_getsetters,
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
    if (PyType_Ready(&SiteType) < 0)
        return NULL;
    if (PyType_Ready(&LigandType) < 0)
        return NULL;
    if (PyType_Ready(&SystemType) < 0)
        return NULL;

    m = PyModule_Create(&ligandmodule);
    if (m == NULL)
        return NULL;

    Py_INCREF(&TransitionType);
    if (PyModule_AddObject(m, "Transition", (PyObject *) &TransitionType) < 0) {
        Py_DECREF(&TransitionType);
        Py_DECREF(m);
        return NULL;
    }
    
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

    Py_INCREF(&SystemType);
    if (PyModule_AddObject(m, "System", (PyObject *) &SystemType) < 0) {
        Py_DECREF(&SystemType);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}
