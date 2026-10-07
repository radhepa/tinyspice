/*
 * analysis.h - DC operating point and transient analysis.
 */
#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <stdio.h>
#include "circuit.h"

typedef struct {
    double reltol;   /* relative tolerance for Newton convergence        */
    double vntol;    /* absolute tolerance on node voltages, V           */
    double abstol;   /* absolute tolerance on branch currents, A         */
    double gmin;     /* tiny conductance from every node to ground, S    */
    int maxiter;     /* Newton iterations allowed per time point         */
    int trap;        /* 1 = trapezoidal rule (default), 0 = backward Euler */
} Options;

typedef struct {
    long steps;      /* accepted time steps            */
    long rejected;   /* steps retried with a smaller h */
    long newton;     /* total Newton iterations        */
} Stats;

Options default_options(void);

/* Solve for the DC operating point. x must hold circuit_size(c) doubles.
 * Returns 0 on success. */
int dc_operating_point(Circuit *c, const Options *o, double *x, Stats *st, char *err, size_t errlen);

/* Run .tran and write a CSV (time, node voltages, branch currents) to out. */
int transient(Circuit *c, const Options *o, FILE *out, Stats *st, char *err, size_t errlen);

/* Write the CSV header / one row (shared with the OP printer). */
void write_header(const Circuit *c, FILE *out);
void write_row(const Circuit *c, double t, const double *x, FILE *out);

#endif
