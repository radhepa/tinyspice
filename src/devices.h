/*
 * devices.h - device equations: source waveforms and the diode.
 */
#ifndef DEVICES_H
#define DEVICES_H

#include "circuit.h"

#define VT_300K 0.025852   /* thermal voltage kT/q at 300 K, volts */

/* Value of an independent source at time t. default_edge is used for PULSE
 * rise/fall times given as 0 (SPICE uses the print step). */
double source_value(const Source *s, double t, double default_edge, double tstop);

/* Collect waveform corners (PULSE edges, PWL points) up to tstop so the
 * time stepper can land exactly on them. Returns how many were written. */
int source_breakpoints(const Source *s, double default_edge, double tstop, double *out, int max);

/* Diode current and its derivative (conductance) at voltage vd. */
void diode_eval(const Model *m, double vd, double *id, double *gd);

/* SPICE-style junction voltage limiting: stops Newton from jumping to huge
 * forward voltages where exp() would overflow. */
double diode_limit(const Model *m, double vnew, double vold);

#endif
