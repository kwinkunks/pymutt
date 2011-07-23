#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdio.h>
#include "NumPy_macros.h"        /* useful macros */

#define MT_STATIC
#include "mtbase.c"

//
// mtft(
//   array data,
//   double dt = 1.0,
//   int kind = 2,           1 for hires, 2 for adaptive
//   int nwin = 5,
//   double npi = 3.0,
//   int paddenlen = 0,      effective if > len(data)
//   int dodof = 0,          produces dof and F keys
//   array lines = 0,        produces lamp, reshaped, and linevar keys.
//                           forces dodof.
// )
//

static PyObject* pymutt_mtft(PyObject* self,
                             PyObject* args,
                             PyObject* keywds
                             )
{
    PyObject *inseries;

    double dt = 1.0;
    int kind = 2;
    int nwin = 5;
    double npi = 3.0;
    int paddedlen = 0;
    int dodof = 0;
    int doweights = 0;
    PyArrayObject* lines = 0;

    int truelen;
    int finallen;

    PyArrayObject* series;
    double* datap;
    double* mtpower;

    int nlines = 0;
    double *linedatap = 0;
    double* flines = 0;

    double W = 0.0;
    double* dof = 0;
    double* Fvalues = 0;
    double complex* lamp = 0;
    double* reshaped = 0;
    double* linevar = 0;

    double* tweights = 0;
    double complex* tspectra = 0;

    PyObject* rd = PyDict_New();

    PyArrayObject* retary;
    PyObject* rv;
    int outlen;
    int woutlen;

    static char *kwlist[] = {"series", "dt", "kind", "nwin", "npi", 
                             "paddedlen", "dodof", "doweights",
                             "lines", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|diidiiiO", kwlist,
                                     &inseries,
                                     &dt, &kind, &nwin, &npi,
                                     &paddedlen, &dodof, &doweights,
                                     &lines))
        return NULL;

    if(lines) {
        if(lines->nd == 0 || lines->dimensions[0] == 0) {
            linedatap = 0;
            nlines = 0;
        } else {
            linedatap = (double *) lines->data;
            nlines = lines->dimensions[0];
        }
    }

    series = (PyArrayObject *) PyArray_FROM_OTF(inseries,
                                                NPY_DOUBLE,
                                                NPY_IN_ARRAY);

    NDIM_CHECK(series, 1);

    truelen = series->dimensions[0];
    datap = (double *) series->data;

    W = npi / (truelen * dt);
    finallen = (truelen >= paddedlen) ? truelen : paddedlen;
    outlen = (finallen / 2) + 1;
    woutlen = nwin * outlen;
    mtpower = (double *) malloc(outlen * sizeof(double)); 

    if(dodof || linedatap) {
        dof = (double *) malloc(outlen * sizeof(double));
        Fvalues = (double *) malloc(outlen * sizeof(double));
        if(linedatap) {
            lamp = (double complex *) malloc(nlines
                                             * sizeof(double complex));
            flines = (double *) malloc(nlines * sizeof(flines[0]));
            reshaped = (double *) malloc(outlen * sizeof(double));
            linevar = (double *) malloc(nlines * sizeof(double));

        }
    }

    if(doweights)
        tweights = (double *) malloc(woutlen * sizeof(double));

    do_mtap_spec(datap, truelen, dt, kind, nwin, npi, finallen,
                 mtpower,
                 dof, Fvalues, nlines, linedatap, lamp, flines,
                 reshaped, linevar,
                 tweights, tspectra);

    retary = (PyArrayObject *) PyArray_SimpleNew(1, &outlen, NPY_DOUBLE);
    memcpy(retary->data, mtpower, outlen * sizeof(double));
    free(mtpower);
    PyDict_SetItemString(rd, "power", (PyObject *) retary);
    Py_DECREF(retary);

    rv = PyFloat_FromDouble(dt);
    PyDict_SetItemString(rd, "dt", rv);
    Py_DECREF(rv);

    rv = PyFloat_FromDouble(1.0 / (finallen * dt));
    PyDict_SetItemString(rd, "df", rv);
    Py_DECREF(rv);

    rv = PyInt_FromLong(finallen);
    PyDict_SetItemString(rd, "n", rv);
    Py_DECREF(rv);

    rv = PyInt_FromLong(outlen);
    PyDict_SetItemString(rd, "nspec", rv);
    Py_DECREF(rv);

    rv = PyInt_FromLong(truelen);
    PyDict_SetItemString(rd, "nin", rv);
    Py_DECREF(rv);

    rv = PyFloat_FromDouble(npi);
    PyDict_SetItemString(rd, "npi", rv);
    Py_DECREF(rv);

    rv = PyInt_FromLong(nwin);
    PyDict_SetItemString(rd, "nwin", rv);
    Py_DECREF(rv);

    rv = PyFloat_FromDouble(W);
    PyDict_SetItemString(rd, "W", rv);
    Py_DECREF(rv);

    rv = PyInt_FromLong(kind);
    PyDict_SetItemString(rd, "kind", rv);
    Py_DECREF(rv);

    if(tweights) {

        int tdims[2];
        tdims[1] = outlen;
        tdims[0] = nwin;
        retary = (PyArrayObject *) PyArray_SimpleNew(2, tdims, NPY_DOUBLE);
        memcpy(retary->data, tweights, woutlen * sizeof(double));
        free(tweights);
        PyDict_SetItemString(rd, "weights", (PyObject *) retary);
        Py_DECREF(retary);
    }

    if(dodof || linedatap) {

        retary = (PyArrayObject *) PyArray_SimpleNew(1, &outlen, NPY_DOUBLE);
        memcpy(retary->data, dof, outlen * sizeof(double));
        free(dof);
        PyDict_SetItemString(rd, "dof", (PyObject *) retary);
        Py_DECREF(retary);
        retary = (PyArrayObject *) PyArray_SimpleNew(1, &outlen, NPY_DOUBLE);
        memcpy(retary->data, Fvalues, outlen * sizeof(double));
        free(Fvalues);
        PyDict_SetItemString(rd, "F", (PyObject *) retary);
        Py_DECREF(retary);

    }

    if(linedatap) {

            // this is the only place (so far) where we assume that
            // the C99 double complex and NPY_CDOUBLE have the same binary
            // format

        retary = (PyArrayObject *) PyArray_SimpleNew(1, &nlines, NPY_CDOUBLE);
        memcpy(retary->data, lamp, nlines * sizeof(double complex));
        free(lamp);
        PyDict_SetItemString(rd, "linea", (PyObject *) retary);
        Py_DECREF(retary);

        retary = (PyArrayObject *) PyArray_SimpleNew(1, &nlines, NPY_DOUBLE);
        memcpy(retary->data, flines, nlines * sizeof(double));
        free(flines);
        PyDict_SetItemString(rd, "linef", (PyObject *) retary);
        Py_DECREF(retary);

        retary = (PyArrayObject *) PyArray_SimpleNew(1, &outlen, NPY_DOUBLE);
        memcpy(retary->data, reshaped, outlen * sizeof(double));
        free(reshaped);
        PyDict_SetItemString(rd, "reshaped", (PyObject *) retary);
        Py_DECREF(retary);

        retary = (PyArrayObject *) PyArray_SimpleNew(1, &nlines, NPY_DOUBLE);
        memcpy(retary->data, linevar, nlines * sizeof(double));
        free(linevar);
        PyDict_SetItemString(rd, "linevar", (PyObject *) retary);
        Py_DECREF(retary);

    }

    return rd;

}

/* Doc strings: */

static char pymutt_mtft_doc[] =
    "Estimate the spectral power (and auxiliary quantities, see below)\n"
    "of a one-dimenstional time series using Thomson's multi-taper algorithm.\n"
    "\n"
    "result_dict = mtft(\n"
    "  array data,      a numpy array\n"
    "  dt = 1.0,\n"
    "  kind = 2,        1 for hires, 2 for adaptive\n"
    "  nwin = 5,        \n"
    "  npi = 3.0,       \n"
    "  paddenlen = 0,   if > length(data), zero-pad to this length\n"
    "  dodof = 0,       produce dof and F arrays\n"
    "  doweights = 0,   returns all of the taper weights\n"
    "  lines = [],      desired frequencies of spectral lines at which \n"
    "                   to estimate complex spectral amplitudes.\n"
    "                   Produces lamp and reshaped spectrum outputs.\n"
    "                   Implies dodof = 1 as well.  Must be a numpy array.\n"
    "  )\n"
    "\n"
    "mtft returns a dictionary containing both fixed and optional items.\n"
    "\n"
    "  Copies of input values:\n"
    "    'dt'        time series sample interval\n"
    "    'nin'       number of points in original time series.\n"
    "    'npi'       dimensionless width of the spectral averaging window\n"
    "    'nwin'      maximum number of tapers to combine\n"
    "    'kind'      spectral algorithm: 1 for hires, 2 for adaptive\n"
    "    'linef'     frequencies at which to estimate spectral line amplitudes\n"
    "                by default this numpy array is empty\n"
    "\n"
    "  Fixed output values:\n"
    "    'df'        frequency-domain spacing of spectral estimates\n"
    "    'n'         number of points in the padded input array\n"
    "    'nspec'     number of points in the spectrum\n"
    "    'W'         dimensionless averaging width in terms of the\n"
    "                possibly-padded input series\n"
    "    'power'     estimated power (not amplitude) spectrum\n"
    "\n"
    "  Optional output values:\n"
    "    'dof'       degrees-of-freedom vs frequency\n"
    "    'F'         F-test for spectral line vs frequency\n"
    "    'weights'   nwin-by-nf array of spectral weights\n"
    "    'linea'     estimated line amplitudes (at frequencies in linef)\n"
    "    'linevar'   estimated line amplitude variance\n"
    "    'reshaped'  spectrum with lines removed and reshaped around holes\n"
    ;

static char pymutt_module_doc[] =
    "Multi-taper 1D digital fourier transform: Code is adapted from\n"
    "a package made freely available by J.M. Lees and J. Park.\n"
    "See http://www.unc.edu/~leesj/mtm/index.html for more.  Also see\n"
    "the file lees_and_park.pdf in the docs directory with this module.\n"
    ;


static PyMethodDef pymutt_methods[] = {

    {"mtft",
     (PyCFunction) pymutt_mtft,
     METH_VARARGS | METH_KEYWORDS,
     pymutt_mtft_doc},

    {NULL, NULL, 0, NULL}
};

void initpymutt(void)
{
    Py_InitModule3("pymutt", pymutt_methods, pymutt_module_doc);
    import_array();   /* required NumPy initialization */
}

