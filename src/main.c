/*
 * tinyspice - a small SPICE-style circuit simulator.
 *
 *   tinyspice circuit.cir            print the DC operating point, run .tran to stdout
 *   tinyspice circuit.cir -o out.csv write the transient results to out.csv
 *   tinyspice circuit.cir --be       use backward Euler instead of the trapezoidal rule
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "analysis.h"
#include "netlist.h"

static void usage(void)
{
    fprintf(stderr, "usage: tinyspice <netlist.cir> [-o out.csv] [--be] [--quiet]\n");
}

int main(int argc, char **argv)
{
    const char *in = NULL, *outpath = NULL;
    int quiet = 0;
    Options o = default_options();
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-o") == 0 && i + 1 < argc)
            outpath = argv[++i];
        else if (strcmp(argv[i], "--be") == 0)
            o.trap = 0;
        else if (strcmp(argv[i], "--quiet") == 0)
            quiet = 1;
        else if (argv[i][0] == '-') {
            usage();
            return 2;
        } else
            in = argv[i];
    }
    if (!in) {
        usage();
        return 2;
    }

    Circuit c;
    char err[512] = "";
    if (parse_netlist(in, &c, err, sizeof err)) {
        fprintf(stderr, "tinyspice: %s\n", err);
        return 1;
    }
    if (!quiet)
        fprintf(stderr, "%s\n  %d nodes, %d elements, %d unknowns\n", c.title, c.nnodes - 1, c.nel,
                circuit_size(&c));

    Stats st = {0, 0, 0};
    clock_t t0 = clock();
    if (!c.has_tran) {
        double *x = calloc((size_t)circuit_size(&c), sizeof(double));
        if (dc_operating_point(&c, &o, x, &st, err, sizeof err)) {
            fprintf(stderr, "tinyspice: %s\n", err);
            return 1;
        }
        printf("DC operating point\n");
        for (int k = 1; k < c.nnodes; k++)
            printf("  v(%s) = %.6g V\n", c.node_names[k], x[k - 1]);
        for (int i = 0; i < c.nel; i++)
            if (c.el[i].branch >= 0 && c.el[i].type != EL_E)
                printf("  i(%s) = %.6g A\n", c.el[i].name, x[c.nnodes - 1 + c.el[i].branch]);
        free(x);
    } else {
        FILE *out = outpath ? fopen(outpath, "w") : stdout;
        if (!out) {
            fprintf(stderr, "tinyspice: cannot write %s\n", outpath);
            return 1;
        }
        int rc = transient(&c, &o, out, &st, err, sizeof err);
        if (outpath)
            fclose(out);
        if (rc) {
            fprintf(stderr, "tinyspice: %s\n", err);
            return 1;
        }
    }
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    if (!quiet)
        fprintf(stderr, "  %ld time steps (%ld retried), %ld Newton iterations, %.3f s\n",
                st.steps, st.rejected, st.newton, secs);
    circuit_free(&c);
    return 0;
}
