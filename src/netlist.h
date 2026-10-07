/*
 * netlist.h - read a SPICE-style text netlist into a Circuit.
 */
#ifndef NETLIST_H
#define NETLIST_H

#include <stddef.h>
#include "circuit.h"

/* Returns 0 on success; on failure fills err with a message and returns -1. */
int parse_netlist(const char *path, Circuit *c, char *err, size_t errlen);

/* Parse a number with an optional SPICE suffix (1k, 4.7u, 2meg, 10pF ...).
 * Returns 0 on success. */
int parse_value(const char *s, double *out);

#endif
