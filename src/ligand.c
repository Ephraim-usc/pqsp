#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "structmember.h"
#include "expm.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>


/******************************************************************************
                          utility functions
******************************************************************************/

//  Increment the counters, lexicographic (dictionary/odometer) style.
static int increment_counters(int n_dims, int *dims, int *counters)
{
    int i;
    for (i = n_dims - 1; i > 0; i--)
    {
        if (++counters[i] < dims[i])
            return 1;
        else
            counters[i] = 0;
    }
    
    if (++counters[0] < dims[0])
        return 1;
    else
        return 0;
}

static int *
PyList2Array_INT(PyObject *listObj) // have to switch from int to long for security!!!
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

static double *
PyList2Array_DOUBLE(PyObject *listObj)
{
    int len, i;
    len = (int) PyList_Size(listObj);
    double *array = calloc(len, sizeof(double));
    for (i = 0; i < len; i++)
      array[i] = (double) PyFloat_AsDouble(PyList_GetItem(listObj, i));
    return array;
}

static PyObject *
Array2PyList_DOUBLE(double *array, int len)
{
    int i;
    PyObject *listObj;
    listObj = PyList_New(len);
    for (i = 0; i < len; i++)
        PyList_SetItem(listObj, i, (PyObject *) Py_BuildValue("d", array[i]));
    return listObj;
}


/******************************************************************************
                     type declarations
******************************************************************************/

typedef struct {
    /* fixed attributes */
    int __max_states__;
    int n_compartments; // number of compartments
    int n_targets; // number of targets
    
    /* variable attributes */
    double ***Qses; // list of P matrices, for each state, for each compartment
    double ***Pses; // list of P matrices, for each state, for each compartment
} Transition;

typedef struct {
    PyObject_HEAD

    /* fixed attributes */
    int __max_states__;
    int n_targets; // number of targets
    int *targets; // list of targets
    double **onses; // list of lists of on rates, for each state
    double **offses; // list of lists of off rates, for each state
    
    /* variable attributes */
    int n_particles;
    int *values; // list of indices of bound targets, instead of targets themselves, for faster transition application
} SiteObject;

typedef struct {
    PyObject_HEAD
    
    /* fixed attributes */
    int n_sites;
    SiteObject **sites;
    int n_forms;
    
    /* state map */
    int *statemap;
    
    /* variable attributes */
    double mpp; // mole per particle
    long n_particles;
    int *compartments; // list of compartments for each particle
    int *forms; // list of forms for each particle
    int *states; // list of states for each particle
} LigandObject;

typedef struct {
    PyObject_HEAD

    /* fixed attributes */
    int n_compartments; // number of compartments
    double *volumes; // list of volumes for each compartment
    int n_analytes; // number of analytes
    int __max_ligands__;
    
    /* variable attributes */
    double **xses; // list of lists of analyte concentrations, for each compartment
    int n_ligands;
    LigandObject **ligands;
} SystemObject;


/******************************************************************************
                           Transition functions
******************************************************************************/

static Transition *
Transition_create(SystemObject *systemObj, SiteObject *siteObj, double t)
{
    Transition *transition = (Transition *)malloc(sizeof(Transition));
    transition->__max_states__ = 1024;
    transition->n_compartments = systemObj->n_compartments;
    transition->n_targets = siteObj->n_targets;
    
    int c, s, i;
    double *ons, *offs, *xs, *Q;
    
    transition->Qses = calloc(systemObj->n_compartments, sizeof(double **));
    transition->Pses = calloc(systemObj->n_compartments, sizeof(double **));
    
    for (c = 0; c < systemObj->n_compartments; c++)
    {
        transition->Qses[c] = calloc(transition->__max_states__, sizeof(double *));
        transition->Pses[c] = calloc(transition->__max_states__, sizeof(double *));
        
        xs = systemObj->xses[c];
        for(s = 0; s < siteObj->__max_states__; s++)
            if(siteObj->onses[s])
            {
                transition->Qses[c][s] = calloc((siteObj->n_targets + 1) * (siteObj->n_targets + 1), sizeof(double));
                
                ons = siteObj->onses[s];
                offs = siteObj->offses[s];
                Q = transition->Qses[c][s];
                
                Q[0] = 0.0;
                for (i = 0; i < siteObj->n_targets; i++)
                {
                    Q[0] -= ons[i] * xs[siteObj->targets[i]] * t;
                    Q[i + 1] = ons[i] * xs[siteObj->targets[i]] * t;
                }
                for (i = 0; i < siteObj->n_targets; i++)
                {
                    Q[(siteObj->n_targets + 1) * (i + 1)] = offs[i] * t;
                    Q[(siteObj->n_targets + 2) * (i + 1)] = - offs[i] * t;
                }
                
                transition->Pses[c][s] = r8mat_expm1(siteObj->n_targets + 1, Q);
            }
    }
    return transition;
}

static int
Transition_free(Transition *transition)
{
    int c, s;
    for (c = 0; c < transition->n_compartments; c++)
    {
        for(s = 0; s < transition->__max_states__; s++)
            if (transition->Qses[c][s])
            {
              free(transition->Qses[c][s]);
              free(transition->Pses[c][s]);
            }
        free(transition->Qses[c]);
        free(transition->Pses[c]);
    }
    free(transition->Qses);
    free(transition->Pses);
    free(transition);
    return 0;
}

static int
Transition_print(Transition *transition)
{   
    int c, s, i;
    double *Q, *P;
    
    for (c = 0; c < transition->n_compartments; c++)
    {
        for(s = 0; s < transition->__max_states__; s++)
            if(transition->Qses[c][s])
            {
                Q = transition->Qses[c][s];
                P = transition->Pses[c][s];
                
                printf("[compartment %d, state %d]\n", c, s);
                printf("[Q] ");
                for (i = 0; i < (transition->n_targets + 1) * (transition->n_targets + 1); i++)
                    printf("%f ", Q[i]);
                printf("\n");
                printf("[P] ");
                for (i = 0; i < (transition->n_targets + 1) * (transition->n_targets + 1); i++)
                    printf("%f ", P[i]);
                printf("\n");
            }
    }
    return 0;
}

static int
Transition_apply(Transition *transition, int compartment, int state, int value)
{
    int value_;
    double tmp;
    
    double *P;
    if (transition->Pses[compartment][state])
        P = transition->Pses[compartment][state] + value * (transition->n_targets + 1); // transition matrix + shift for starting state x = transition vector for x
    else
        P = transition->Pses[compartment][0] + value * (transition->n_targets + 1);
    
    tmp = drand48();
    for (value_ = 0; tmp -= P[value_], tmp > 0; value_++);
    
    return value_;
}


/******************************************************************************
                                the Site type
******************************************************************************/

static void
Site_dealloc(SiteObject *self)
{
    int i;
    free(self->targets);
    for (i = 0; i < self->__max_states__; i++) if (self->onses[i]) free(self->onses[i]);
    free(self->onses);
    for (i = 0; i < self->__max_states__; i++) if (self->offses[i]) free(self->offses[i]);
    free(self->offses);
    free(self->values);
    
    Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *
Site_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    SiteObject *self;
    self = (SiteObject *) type->tp_alloc(type, 0);
    if (self != NULL)
    {
        self->__max_states__ = 1024;
        
        self->n_targets = 0;
        self->targets = NULL;
        self->onses = NULL;
        self->offses = NULL;
        
        self->n_particles = 0;
        self->values = NULL;
    }
    return (PyObject *) self;
}

static int
Site_init(SiteObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *targetsObj, *onsObj, *offsObj;
    
    static char *kwlist[] = {"targets", "ons", "offs", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O!O!O!", kwlist, &PyList_Type, &targetsObj, &PyList_Type, &onsObj, &PyList_Type, &offsObj))
        return -1;
    
    self->n_targets = (int) PyList_Size(targetsObj);
    
    self->targets = PyList2Array_INT(targetsObj);
    
    self->onses = calloc(self->__max_states__, sizeof(double *));
    self->onses[0] = PyList2Array_DOUBLE(onsObj);
    
    self->offses = calloc(self->__max_states__, sizeof(double *));
    self->offses[0] = PyList2Array_DOUBLE(offsObj);
    
    //self->transition = PyObject_Call((PyObject *) &TransitionType, PyTuple_New(0), Py_BuildValue("{s:i, s:i, s:O}", "n_compartments", self));
    return 0;
}

static PyMemberDef Site_members[] = {
    {"n_targets", T_INT, offsetof(SiteObject, n_targets), READONLY, "number of targets"},
    {NULL}  /* Sentinel */
};

static PyObject *
Site_getvalues(SiteObject *self, void *closure)
{
    PyObject *valuesObj = Array2PyList_INT(self->values, self->n_particles);
    return Py_NewRef(valuesObj);
}

static PyGetSetDef Site_getsetters[] = {
    {"values", (getter) Site_getvalues, NULL, "list of values for each particle", NULL},
    {NULL}  /* Sentinel */
};


static PyObject *
Site_print(SiteObject *self, PyObject *Py_UNUSED(ignored))
{
    int state, i;
    
    printf("targets:\n");
    for (i = 0; i < self->n_targets; i++)
        printf("%d ", self->targets[i]);
    printf("\n\n");
    
    printf("on rates:\n");
    for (state = 0; state < self->__max_states__; state++)
        if (self->onses[state] != NULL)
        {
            printf("[state %d] ", state);
            for (i = 0; i < self->n_targets; i++)
                printf("%f ", self->onses[state][i]);
            printf("\n");
        }
    printf("\n");

    printf("off rates:\n");
    for (state = 0; state < self->__max_states__; state++)
        if (self->offses[state] != NULL)
        {
            printf("[state %d] ", state);
            for (i = 0; i < self->n_targets; i++)
                printf("%f ", self->offses[state][i]);
            printf("\n");
        }
    printf("\n");
    
    Py_RETURN_NONE;
}

static PyObject *
Site_set_state(SiteObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *onsObj, *offsObj;
    int state, i;
    
    static char *kwlist[] = {"state", "ons", "offs", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|iO!O!", kwlist, &state, &PyList_Type, &onsObj, &PyList_Type, &offsObj))
        Py_RETURN_NONE;
    
    if (!self->onses[state]) // if already added, then replace
        self->onses[state] = calloc(self->n_targets, sizeof(double));
    for (i = 0; i < self->n_targets; i++)
        self->onses[state][i] = (double) PyFloat_AsDouble(PyList_GetItem(onsObj, i));
    
    if (!self->offses[state]) // if already added, then replace
        self->offses[state] = calloc(self->n_targets, sizeof(double));
    for (i = 0; i < self->n_targets; i++)
        self->offses[state][i] = (double) PyFloat_AsDouble(PyList_GetItem(offsObj, i));
    
    Py_RETURN_NONE;
}

static int
Site_add_particles(SiteObject *self, int n)
{
    int n_particles_old;
    int p;
    int *values_old;
    
    n_particles_old = self->n_particles;
    values_old = self->values;
    
    self->n_particles += n;
    self->values = calloc(self->n_particles, sizeof(int));
    for(p = 0; p < n_particles_old; p++)
        self->values[p] = values_old[p];
    free(values_old);
    
    return 0;
}

static PyMethodDef Site_methods[] = {
    {"print", (PyCFunction) Site_print, METH_NOARGS, "print the Site"},
    {"set_state", (PyCFunction) Site_set_state, METH_VARARGS | METH_KEYWORDS, "set a state (i.e., on and off rates) of the Site"},
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
    .tp_getset = Site_getsetters,
    .tp_methods = Site_methods,
    .tp_dealloc = (destructor) Site_dealloc,
};



/******************************************************************************
                                the Ligand type
******************************************************************************/

static void
Ligand_dealloc(LigandObject *self)
{
    free(self->compartments);
    free(self->states);
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
        
        self->statemap = NULL;
        
        self->n_particles = 0;
        self->mpp = 0.0;
        self->compartments = NULL;
        self->forms = NULL;
        self->states = NULL;
    }
    return (PyObject *) self;
}

static int
Ligand_init(LigandObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *sitesObj;
    int st;
    int len_statemap;
    
    static char *kwlist[] = {"sites", "n_forms", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O!i", kwlist, &PyList_Type, &sitesObj, &self->n_forms))
        return -1;
    
    self->n_sites = (int) PyList_Size(sitesObj);
    self->sites = calloc(self->n_sites, sizeof(SiteObject *));
    for(st = 0; st < self->n_sites; st++)
        self->sites[st] = (SiteObject *) PyList_GetItem(sitesObj, st);
    
    len_statemap = self->n_forms;
    for (st = 0; st < self->n_sites; st++)
        len_statemap *= self->sites[st]->n_targets + 1;
    self->statemap = calloc(len_statemap, sizeof(int));
    
    return 0;
}

static PyMemberDef Ligand_members[] = {
    {"n_sites", T_INT, offsetof(LigandObject, n_sites), READONLY, "number of binding sites"},
    {"n_forms", T_INT, offsetof(LigandObject, n_forms), READONLY, "number of forms"},
    {"mpp", T_DOUBLE, offsetof(LigandObject, mpp), READONLY, "mole per particle"},
    {"n_particles", T_LONG, offsetof(LigandObject, n_particles), READONLY, "number of particles"},
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
Ligand_getstatemap(LigandObject *self, void *closure)
{
    int st;
    int len_statemap = self->n_forms;
    for (st = 0; st < self->n_sites; st++)
        len_statemap *= self->sites[st]->n_targets + 1;
    
    PyObject *statemapObj = Array2PyList_INT(self->statemap, len_statemap);
    return Py_NewRef(statemapObj);
}

static PyObject *
Ligand_getcompartments(LigandObject *self, void *closure)
{
    PyObject *compartmentsObj = Array2PyList_INT(self->compartments, self->n_particles);
    return Py_NewRef(compartmentsObj);
}

static PyObject *
Ligand_getforms(LigandObject *self, void *closure)
{
    PyObject *formsObj = Array2PyList_INT(self->forms, self->n_particles);
    return Py_NewRef(formsObj);
}

static PyObject *
Ligand_getstates(LigandObject *self, void *closure)
{
    PyObject *statesObj = Array2PyList_INT(self->states, self->n_particles);
    return Py_NewRef(statesObj);
}

static PyGetSetDef Ligand_getsetters[] = {
    {"sites", (getter) Ligand_getsites, NULL, "binding sites", NULL},
    {"statemap", (getter) Ligand_getstatemap, NULL, "state map", NULL},
    {"compartments", (getter) Ligand_getcompartments, NULL, "list of compartments for each particle", NULL},
    {"forms", (getter) Ligand_getforms, NULL, "list of forms for each particle", NULL},
    {"states", (getter) Ligand_getstates, NULL, "list of states, for each particle", NULL},
    {NULL}  /* Sentinel */
};

/*
static PyObject *
Ligand_update_states(LigandObject *self, PyObject *Py_UNUSED(ignored))
{
    long p;
    for(p = 0; p < self->n_particles; p++)
    {
        
    }
    
    Py_RETURN_NONE;
}
*/

static PyObject *
Ligand_set_mpp(LigandObject *self, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"mpp", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|d", kwlist, &self->mpp))
        Py_RETURN_NONE;
    
    Py_RETURN_NONE;
}

static PyObject *
Ligand_define_state(LigandObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *formsObj;
    PyObject *valuesesObj, *valuesObj;
    int st, s, i, idx;
    int **valueses, *dims, *counters;
    int n_dims;
    
    static char *kwlist[] = {"forms", "valueses", "state", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O!O!i", kwlist, &PyList_Type, &formsObj, &PyList_Type, &valuesesObj, &s))
        Py_RETURN_NONE;
    
    n_dims = 1 + self->n_sites;
    
    dims = calloc(1 + self->n_sites, sizeof(int));
    dims[0] = self->n_forms;
    for (st = 0; st < self->n_sites; st++)
        dims[st + 1] = self->sites[st]->n_targets + 1;
    
    valueses = (int **)calloc(1 + self->n_sites, sizeof(int *));
    valueses[0] = PyList2Array_INT(formsObj);
    for (st = 0; st < self->n_sites; st++)
        valueses[st + 1] = PyList2Array_INT(PyList_GetItem(valuesesObj, st));
    
    counters = (int *)calloc(n_dims, sizeof(int));
    do
    {
        for (i = 1; i < n_dims; i++)
        {
            printf("%d ", valueses[i][counters[i]]);
        }
        printf("\n");
        
        idx = valueses[0][counters[0]];
        for (i = 1; i < n_dims; i++)
        {
            idx *= dims[i];
            idx += valueses[i][counters[i]];
        }
        self->statemap[idx] = s;
    }
    while (increment_counters(n_dims, dims, counters));
    
    Py_RETURN_NONE;
}

static PyObject *
Ligand_add_particles(LigandObject *self, PyObject *args, PyObject *kwds)
{
    long n;
    int compartment, form;
    long p; int st;
    long n_particles_old;
    int *compartments_old;
    int *forms_old;
    int *states_old;
    
    static char *kwlist[] = {"n", "compartment", "form", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|lii", kwlist, &n, &compartment, &form))
        Py_RETURN_NONE;
    
    n_particles_old = self->n_particles;
    compartments_old = self->compartments;
    forms_old = self->forms;
    states_old = self->states;
    
    self->n_particles += n;
    
    self->compartments = calloc(self->n_particles, sizeof(int));
    for(p = 0; p < n_particles_old; p++)
        self->compartments[p] = compartments_old[p];
    for(p = n_particles_old; p < self->n_particles; p++)
        self->compartments[p] = compartment;
    free(compartments_old);
    
    self->forms = calloc(self->n_particles, sizeof(int));
    for(p = 0; p < n_particles_old; p++)
        self->forms[p] = forms_old[p];
    for(p = n_particles_old; p < self->n_particles; p++)
        self->forms[p] = form;
    free(forms_old);
    
    self->states = calloc(self->n_particles, sizeof(int));
    for(p = 0; p < n_particles_old; p++)
        self->states[p] = states_old[p];
    free(states_old);
    
    for(st = 0; st < self->n_sites; st++)
        Site_add_particles(self->sites[st], n);
    
    Py_RETURN_NONE;
}

static PyMethodDef Ligand_methods[] = {
    {"set_mpp", (PyCFunction) Ligand_set_mpp, METH_VARARGS | METH_KEYWORDS, "set the mpp value of the ligand"},
    {"add_particles", (PyCFunction) Ligand_add_particles, METH_VARARGS | METH_KEYWORDS, "add particles of the ligand"},
    {"define_state", (PyCFunction) Ligand_define_state, METH_VARARGS | METH_KEYWORDS, "define state"},
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
    .tp_methods = Ligand_methods,
    .tp_dealloc = (destructor) Ligand_dealloc,
    .tp_getset = Ligand_getsetters,
};


/******************************************************************************
                                the System type
******************************************************************************/

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
    PyObject *volumesObj;
    int c;
    
    static char *kwlist[] = {"volumes", "n_analytes", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O!i", kwlist, &PyList_Type, &volumesObj, &self->n_analytes))
        return -1;
    
    self->n_compartments = (int) PyList_Size(volumesObj);
    self->volumes = PyList2Array_DOUBLE(volumesObj);
    
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
System_getvolumes(SystemObject *self, void *closure)
{
    PyObject *volumesObj = Array2PyList_DOUBLE(self->volumes, self->n_compartments);
    return Py_NewRef(volumesObj);
}

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
    {"volumes", (getter) System_getvolumes, NULL, "list volumes of each compartment", NULL},
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

static PyObject *
System_interact(SystemObject *self, PyObject *args, PyObject *kwds)
{
    LigandObject *ligandObj;
    double t;
    
    static char *kwlist[] = {"ligand", "t", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O!d", kwlist, &LigandType, &ligandObj, &t))
        Py_RETURN_NONE;
    
    int *compartments = ligandObj->compartments;
    int *states = ligandObj->states;
    
    SiteObject *siteObj;
    Transition *transition;
    int c, a, s, st, p, value, value_; // index of compartment, analyte, state, site, particle.
    int *targets, *values;
    
    int **deltas = calloc(self->n_compartments, sizeof(int *));
    for (c = 0; c < self->n_compartments; c++)
        deltas[c] = calloc(self->n_analytes, sizeof(int));
    
    for (st = 0; st < ligandObj->n_sites; st++)
    {
        siteObj = ligandObj->sites[st];
        transition = Transition_create(self, siteObj, t);
        targets = siteObj->targets - 1; // for faster mapping from value to target
        values = siteObj->values;
        
        for (p = 0; p < siteObj->n_particles; p++)
        {
            c = compartments[p]; s = states[p]; value = values[p];
            value_ = Transition_apply(transition, c, s, value);
            
            values[p] = value_;
            if (value) deltas[c][targets[value]] -= 1;
            if (value_) deltas[c][targets[value_]] += 1;
        }
        Transition_free(transition);
    }

    /*
    for (c = 0; c < self->n_compartments; c++)
    {
        for (a = 0; a < self->n_analytes + 1; a++)
            printf("%f ", (double) deltas[c][a]);
        printf("\n");
    }
    */
    
    for (c = 0; c < self->n_compartments; c++)
        for (a = 0; a < self->n_analytes; a++)
            self->xses[c][a] = fmax(0.0, self->xses[c][a] - (double) deltas[c][a] * ligandObj->mpp / self->volumes[c]);
    
    for (c = 0; c < self->n_compartments; c++)
        free(deltas[c]);
    free(deltas);
    
    Py_RETURN_NONE;
}

static PyMethodDef System_methods[] = {
    {"add_x", (PyCFunction) System_add_x, METH_VARARGS | METH_KEYWORDS, "add analyte"},
    {"add_ligand", (PyCFunction) System_add_ligand, METH_VARARGS | METH_KEYWORDS, "add ligand"},
    {"interact", (PyCFunction) System_interact, METH_VARARGS | METH_KEYWORDS, "interact with a ligand"},
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
    if (PyType_Ready(&SiteType) < 0)
        return NULL;
    if (PyType_Ready(&LigandType) < 0)
        return NULL;
    if (PyType_Ready(&SystemType) < 0)
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

    Py_INCREF(&SystemType);
    if (PyModule_AddObject(m, "System", (PyObject *) &SystemType) < 0) {
        Py_DECREF(&SystemType);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}
