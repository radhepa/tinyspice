/*
 * devices.c - source waveforms and the diode model.
 */
#include <math.h>
#include "devices.h"

#define PI 3.14159265358979323846

double source_value(const Source *s, double t, double default_edge, double tstop)
{
    switch (s->kind) {
    case SRC_DC:
        return s->dc;
    case SRC_PULSE: {
        double v1 = s->p[0], v2 = s->p[1], td = s->p[2];
        double tr = s->p[3] > 0 ? s->p[3] : default_edge;
        double tf = s->p[4] > 0 ? s->p[4] : default_edge;
        double pw = s->p[5] > 0 ? s->p[5] : tstop;
        double per = s->p[6] > 0 ? s->p[6] : tstop + tr + pw + tf;
        if (t < td)
            return v1;
        double tt = fmod(t - td, per);
        if (tt < tr)
            return v1 + (v2 - v1) * tt / tr;
        if (tt < tr + pw)
            return v2;
        if (tt < tr + pw + tf)
            return v2 + (v1 - v2) * (tt - tr - pw) / tf;
        return v1;
    }
    case SRC_SIN: {
        double vo = s->p[0], va = s->p[1], f = s->p[2], td = s->p[3];
        if (t < td)
            return vo;
        return vo + va * sin(2.0 * PI * f * (t - td));
    }
    case SRC_PWL: {
        if (t <= s->pwl_t[0])
            return s->pwl_v[0];
        for (int i = 1; i < s->npwl; i++) {
            if (t <= s->pwl_t[i]) {
                double f = (t - s->pwl_t[i - 1]) / (s->pwl_t[i] - s->pwl_t[i - 1]);
                return s->pwl_v[i - 1] + f * (s->pwl_v[i] - s->pwl_v[i - 1]);
            }
        }
        return s->pwl_v[s->npwl - 1];
    }
    }
    return 0.0;
}

int source_breakpoints(const Source *s, double default_edge, double tstop, double *out, int max)
{
    int n = 0;
    if (s->kind == SRC_PULSE) {
        double td = s->p[2];
        double tr = s->p[3] > 0 ? s->p[3] : default_edge;
        double tf = s->p[4] > 0 ? s->p[4] : default_edge;
        double pw = s->p[5] > 0 ? s->p[5] : tstop;
        double per = s->p[6] > 0 ? s->p[6] : tstop + tr + pw + tf;
        for (double t0 = td; t0 < tstop && n + 4 <= max; t0 += per) {
            out[n++] = t0;
            out[n++] = t0 + tr;
            out[n++] = t0 + tr + pw;
            out[n++] = t0 + tr + pw + tf;
        }
    } else if (s->kind == SRC_PWL) {
        for (int i = 0; i < s->npwl && n < max; i++)
            out[n++] = s->pwl_t[i];
    }
    return n;
}

/* exp() that grows linearly instead of overflowing for very large arguments */
static double safe_exp(double x, double *dexp)
{
    const double xmax = 80.0;
    if (x > xmax) {
        double e = exp(xmax);
        *dexp = e;
        return e * (1.0 + (x - xmax));
    }
    double e = exp(x);
    *dexp = e;
    return e;
}

void diode_eval(const Model *m, double vd, double *id, double *gd)
{
    double nvt = m->n * VT_300K;
    double de;
    double e = safe_exp(vd / nvt, &de);
    *id = m->is * (e - 1.0);
    *gd = m->is * de / nvt;
}

double diode_limit(const Model *m, double vnew, double vold)
{
    double nvt = m->n * VT_300K;
    double vcrit = nvt * log(nvt / (sqrt(2.0) * m->is));
    if (vnew > vcrit && fabs(vnew - vold) > 2.0 * nvt) {
        if (vold > 0.0) {
            double arg = 1.0 + (vnew - vold) / nvt;
            vnew = arg > 0.0 ? vold + nvt * log(arg) : vcrit;
        } else {
            vnew = nvt * log(vnew / nvt);
        }
    }
    return vnew;
}
