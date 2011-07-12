//@+leo-ver=4
//@+node:@file mtbase.h
//@@language c
#ifndef MTBASE_H
#include <complex.h>

void do_mtap_spec(double* data, int npoints, double dt, int kind,
                  int nwin, double npi, int nextd,
                  double* ospec,
                  double* dof,
                  double* Fvalues,
                  int nlines,
                  double* lines,
                  double complex* lamp,
                  double* flamp,
                  double* sigsq,
                  double* linevar,
                  double* tweights,
                  double complex* tspectra
                  );

#define MTBASE_H
#endif
//@-node:@file mtbase.h
//@-leo
