/*
 * analysis.c - Modified Nodal Analysis (MNA), Newton-Raphson, and time stepping.
 *
 * MNA in one paragraph: write Kirchhoff's current law at every node except
 * ground. Unknowns are the node voltages, plus one extra unknown (a branch
 * current) for every voltage source and inductor, because their current can't
 * be written from the voltages alone. That gives a square system A x = b.
 * Each element "stamps" a few numbers into A and b; solving gives x.
 *
 * Nonlinear elements (diodes) are linearized around a guess and the system is
 * solved repeatedly until the answer stops changing (Newton-Raphson).
 *
 * Capacitors and inductors are replaced at each time step by a "companion
 * model": a resistor plus a current/voltage source that encodes the past.
 * With the trapezoidal rule this is exactly the trapezoid area rule from
 * calculus applied to i = C dv/dt and v = L di/dt.
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "analysis.h"
#include "devices.h"
#include "solver.h"

typedef enum { MODE_OP, MODE_TRAN } Mode;

typedef struct {
    int n;          /* number of unknowns */
    int nn;         /* number of non-ground nodes */
    double *A, *b, *x, *xprev;
    int *piv;
} Work;

Options default_options(void)
{
    Options o = {1e-3, 1e-6, 1e-9, 1e-12, 100, 1};
    return o;
}

/* ------------------------------------------------------------------ stamps */

static double volt(const Work *w, const double *x, int node)
{
    (void)w;
    return node > 0 ? x[node - 1] : 0.0;
}

static void add_a(Work *w, int row, int col, double v)    /* row/col are unknown indices */
{
    if (row >= 0 && col >= 0)
        w->A[row * w->n + col] += v;
}

static void add_b(Work *w, int row, double v)
{
    if (row >= 0)
        w->b[row] += v;
}

/* conductance g between nodes p and q */
static void stamp_g(Work *w, int p, int q, double g)
{
    add_a(w, p - 1, p - 1, g);
    add_a(w, q - 1, q - 1, g);
    add_a(w, p - 1, q - 1, -g);
    add_a(w, q - 1, p - 1, -g);
}

/* current i injected INTO node p and drawn OUT of node q */
static void stamp_i(Work *w, int p, int q, double i)
{
    add_b(w, p - 1, i);
    add_b(w, q - 1, -i);
}

/* branch current unknown k flows from node p, through the element, to node q */
static void stamp_branch_kcl(Work *w, int k, int p, int q)
{
    add_a(w, p - 1, k, 1.0);
    add_a(w, q - 1, k, -1.0);
}

typedef struct {
    Mode mode;
    double t;       /* time of the point being solved */
    double h;       /* step size (transient)          */
    int trap;       /* integration method for this step */
    double gmin;
} Point;

/* Build A and b around the guess x. */
static void load(Circuit *c, Work *w, const double *x, const Point *pt)
{
    memset(w->A, 0, sizeof(double) * (size_t)(w->n * w->n));
    memset(w->b, 0, sizeof(double) * (size_t)w->n);
    for (int k = 1; k <= w->nn; k++)            /* gmin keeps every node anchored */
        add_a(w, k - 1, k - 1, pt->gmin);

    for (int i = 0; i < c->nel; i++) {
        Element *e = &c->el[i];
        int a = e->n1, b = e->n2;
        int br = e->branch >= 0 ? w->nn + e->branch : -1;
        switch (e->type) {
        case EL_R:
            stamp_g(w, a, b, 1.0 / e->value);
            break;
        case EL_C:
            if (pt->mode == MODE_TRAN) {
                /* i_new = geq * v_new - ieq */
                double geq = (pt->trap ? 2.0 : 1.0) * e->value / pt->h;
                double ieq = geq * e->v_prev + (pt->trap ? e->i_prev : 0.0);
                stamp_g(w, a, b, geq);
                stamp_i(w, a, b, ieq);
            }
            break;                              /* DC: a capacitor is an open circuit */
        case EL_L:
            stamp_branch_kcl(w, br, a, b);
            add_a(w, br, a - 1, 1.0);
            add_a(w, br, b - 1, -1.0);
            if (pt->mode == MODE_TRAN) {
                /* v_new - req * i_new = -req * i_prev [- v_prev for trapezoidal] */
                double req = (pt->trap ? 2.0 : 1.0) * e->value / pt->h;
                add_a(w, br, br, -req);
                add_b(w, br, -req * e->i_prev - (pt->trap ? e->v_prev : 0.0));
            }
            break;                              /* DC: a short circuit (v = 0) */
        case EL_V:
            stamp_branch_kcl(w, br, a, b);
            add_a(w, br, a - 1, 1.0);
            add_a(w, br, b - 1, -1.0);
            add_b(w, br, source_value(&e->src, pt->t, c->tstep, c->tstop));
            break;
        case EL_I: {
            double val = source_value(&e->src, pt->t, c->tstep, c->tstop);
            stamp_i(w, b, a, val);              /* flows from n+ through the source to n- */
            break;
        }
        case EL_D: {
            const Model *m = &c->models[e->model];
            double vd = volt(w, x, a) - volt(w, x, b);
            vd = diode_limit(m, vd, e->vd_iter);
            e->vd_iter = vd;
            double id, gd;
            diode_eval(m, vd, &id, &gd);
            /* linear model around vd: i = gd * v + (id - gd * vd) */
            stamp_g(w, a, b, gd);
            stamp_i(w, b, a, id - gd * vd);
            break;
        }
        case EL_S: {
            const Model *m = &c->models[e->model];
            double vc = volt(w, x, e->nc1) - volt(w, x, e->nc2);
            stamp_g(w, a, b, vc > m->vt ? 1.0 / m->ron : 1.0 / m->roff);
            break;
        }
        case EL_E:
            stamp_branch_kcl(w, br, a, b);
            add_a(w, br, a - 1, 1.0);
            add_a(w, br, b - 1, -1.0);
            add_a(w, br, e->nc1 - 1, -e->value);
            add_a(w, br, e->nc2 - 1, e->value);
            break;
        }
    }
}

static int is_nonlinear(const Circuit *c)
{
    for (int i = 0; i < c->nel; i++)
        if (c->el[i].type == EL_D || c->el[i].type == EL_S)
            return 1;
    return 0;
}

/* Newton-Raphson at one time point. On entry w->x is the starting guess;
 * on success it holds the solution. Returns iterations used, or -1. */
static int newton(Circuit *c, Work *w, const Point *pt, const Options *o, Stats *st,
                  char *err, size_t errlen)
{
    int nonlinear = is_nonlinear(c);
    for (int it = 1; it <= o->maxiter; it++) {
        load(c, w, w->x, pt);
        int sing = lu_factor(w->A, w->n, w->piv);
        st->newton++;
        if (sing) {
            int k = sing - 1;
            if (k < w->nn)
                snprintf(err, errlen, "singular matrix at node '%s' (floating node or loop of voltage sources?)",
                         c->node_names[k + 1]);
            else
                snprintf(err, errlen, "singular matrix at a branch current (loop of voltage sources/inductors?)");
            return -2;
        }
        lu_solve(w->A, w->n, w->piv, w->b);     /* w->b now holds the new x */

        int converged = 1;
        for (int k = 0; k < w->n; k++) {
            double tol = o->reltol * fmax(fabs(w->b[k]), fabs(w->x[k]))
                         + (k < w->nn ? o->vntol : o->abstol);
            if (fabs(w->b[k] - w->x[k]) > tol)
                converged = 0;
        }
        memcpy(w->x, w->b, sizeof(double) * (size_t)w->n);
        if (!nonlinear || (converged && it > 1))
            return it;
    }
    return -1;
}

static Work work_new(const Circuit *c)
{
    Work w;
    w.nn = c->nnodes - 1;
    w.n = circuit_size(c);
    w.A = malloc(sizeof(double) * (size_t)(w.n * w.n));
    w.b = malloc(sizeof(double) * (size_t)w.n);
    w.x = calloc((size_t)w.n, sizeof(double));
    w.xprev = calloc((size_t)w.n, sizeof(double));
    w.piv = malloc(sizeof(int) * (size_t)w.n);
    return w;
}

static void work_free(Work *w)
{
    free(w->A);
    free(w->b);
    free(w->x);
    free(w->xprev);
    free(w->piv);
}

/* -------------------------------------------------------- operating point */

static int solve_op(Circuit *c, Work *w, const Options *o, Stats *st, char *err, size_t errlen)
{
    Point pt = {MODE_OP, 0.0, 0.0, 0, o->gmin};
    for (int i = 0; i < c->nel; i++)
        c->el[i].vd_iter = 0.0;
    int rc = newton(c, w, &pt, o, st, err, errlen);
    if (rc == -2)
        return -1;
    if (rc > 0)
        return 0;
    /* gmin stepping: start with big conductances to ground (an easy, almost
     * linear problem) and remove them gradually, each solution seeding the next */
    memset(w->x, 0, sizeof(double) * (size_t)w->n);
    for (double g = 1e-2; g >= o->gmin * 0.99; g /= 10.0) {
        pt.gmin = g;
        if (newton(c, w, &pt, o, st, err, errlen) < 0) {
            snprintf(err, errlen, "DC operating point did not converge (even with gmin stepping)");
            return -1;
        }
    }
    return 0;
}

int dc_operating_point(Circuit *c, const Options *o, double *x, Stats *st, char *err, size_t errlen)
{
    Work w = work_new(c);
    int rc = solve_op(c, &w, o, st, err, errlen);
    if (rc == 0)
        memcpy(x, w.x, sizeof(double) * (size_t)w.n);
    work_free(&w);
    return rc;
}

/* ----------------------------------------------------------------- output */

void write_header(const Circuit *c, FILE *out)
{
    fprintf(out, "time");
    for (int k = 1; k < c->nnodes; k++)
        fprintf(out, ",v(%s)", c->node_names[k]);
    for (int i = 0; i < c->nel; i++)
        if (c->el[i].branch >= 0 && c->el[i].type != EL_E)
            fprintf(out, ",i(%s)", c->el[i].name);
    fprintf(out, "\n");
}

void write_row(const Circuit *c, double t, const double *x, FILE *out)
{
    int nn = c->nnodes - 1;
    fprintf(out, "%.9e", t);
    for (int k = 0; k < nn; k++)
        fprintf(out, ",%.9e", x[k]);
    for (int i = 0; i < c->nel; i++)
        if (c->el[i].branch >= 0 && c->el[i].type != EL_E)
            fprintf(out, ",%.9e", x[nn + c->el[i].branch]);
    fprintf(out, "\n");
}

/* -------------------------------------------------------------- transient */

static int cmp_double(const void *p, const void *q)
{
    double a = *(const double *)p, b = *(const double *)q;
    return (a > b) - (a < b);
}

/* After a step is accepted, remember each energy-storing element's state. */
static void accept_step(Circuit *c, const Work *w, const Point *pt)
{
    for (int i = 0; i < c->nel; i++) {
        Element *e = &c->el[i];
        double v = volt(w, w->x, e->n1) - volt(w, w->x, e->n2);
        if (e->type == EL_C) {
            double geq = (pt->trap ? 2.0 : 1.0) * e->value / pt->h;
            double ieq = geq * e->v_prev + (pt->trap ? e->i_prev : 0.0);
            e->i_prev = geq * v - ieq;
            e->v_prev = v;
        } else if (e->type == EL_L) {
            e->i_prev = w->x[w->nn + e->branch];
            e->v_prev = v;
        }
    }
}

int transient(Circuit *c, const Options *o, FILE *out, Stats *st, char *err, size_t errlen)
{
    if (!c->has_tran || c->tstep <= 0 || c->tstop <= 0) {
        snprintf(err, errlen, "no valid .tran line");
        return -1;
    }
    Work w = work_new(c);
    if (solve_op(c, &w, o, st, err, errlen)) {
        work_free(&w);
        return -1;
    }
    /* initial conditions from the operating point */
    for (int i = 0; i < c->nel; i++) {
        Element *e = &c->el[i];
        double v = volt(&w, w.x, e->n1) - volt(&w, w.x, e->n2);
        e->v_prev = v;
        e->i_prev = e->type == EL_L ? w.x[w.nn + e->branch] : 0.0;
        e->vd_iter = v;
    }

    /* breakpoints: corners of the source waveforms, sorted */
    int nbp = 0, cap = 1 << 16;
    double *bp = malloc(sizeof(double) * (size_t)cap);
    for (int i = 0; i < c->nel; i++)
        if (c->el[i].type == EL_V || c->el[i].type == EL_I)
            nbp += source_breakpoints(&c->el[i].src, c->tstep, c->tstop, bp + nbp, cap - nbp);
    bp[nbp++] = c->tstop;
    qsort(bp, (size_t)nbp, sizeof(double), cmp_double);

    write_header(c, out);
    if (c->tstart <= 0.0)
        write_row(c, 0.0, w.x, out);

    double t = 0.0, hmax = c->tstep, h = hmax;
    double hmin = c->tstep * 1e-9;
    int ib = 0, after_corner = 1;           /* first step uses backward Euler */
    int rc = 0;
    while (t < c->tstop * (1.0 - 1e-12)) {
        while (ib < nbp && bp[ib] <= t * (1.0 + 1e-12) + 1e-21)
            ib++;
        double step = fmin(h, c->tstop - t);
        int lands_on_corner = 0;
        if (ib < nbp && t + step >= bp[ib] * (1.0 - 1e-12)) {
            step = bp[ib] - t;
            lands_on_corner = 1;
        }
        Point pt = {MODE_TRAN, t + step, step, o->trap && !after_corner, o->gmin};
        memcpy(w.xprev, w.x, sizeof(double) * (size_t)w.n);
        int it = newton(c, &w, &pt, o, st, err, errlen);
        if (it == -2) {
            rc = -1;
            break;
        }
        if (it < 0) {                       /* no convergence: retry with half the step */
            memcpy(w.x, w.xprev, sizeof(double) * (size_t)w.n);
            st->rejected++;
            h = step / 2.0;
            if (h < hmin) {
                snprintf(err, errlen, "time step too small at t = %g s", t);
                rc = -1;
                break;
            }
            continue;
        }
        accept_step(c, &w, &pt);
        t += step;
        st->steps++;
        if (t >= c->tstart)
            write_row(c, t, w.x, out);
        /* integrate across a waveform corner with backward Euler for one
         * step: the trapezoidal rule rings ("trap ringing") at discontinuities */
        after_corner = lands_on_corner;
        if (h < hmax)
            h = fmin(hmax, 2.0 * h);
    }
    free(bp);
    work_free(&w);
    return rc;
}
