/*
 * solver.c - Gaussian elimination as an LU factorization.
 *
 * For every column k:
 *   1. pick the row with the largest |value| in column k (partial pivoting,
 *      which keeps round-off errors small) and swap it to row k;
 *   2. subtract multiples of row k from the rows below so column k becomes
 *      zero under the diagonal. The multipliers are stored where the zeros
 *      would be: that is the L matrix; what is left on and above the diagonal
 *      is U.
 * Then A x = b is solved in two cheap passes: L y = P b (forward), U x = y (back).
 */
#include <math.h>
#include <stdlib.h>
#include "solver.h"

int lu_factor(double *a, int n, int *piv)
{
    for (int i = 0; i < n; i++)
        piv[i] = i;

    for (int k = 0; k < n; k++) {
        /* 1. pivot search */
        int p = k;
        double best = fabs(a[k * n + k]);
        for (int i = k + 1; i < n; i++) {
            double v = fabs(a[i * n + k]);
            if (v > best) {
                best = v;
                p = i;
            }
        }
        if (best < 1e-300)
            return k + 1;                     /* singular: e.g. a floating node */
        if (p != k) {
            for (int j = 0; j < n; j++) {
                double t = a[k * n + j];
                a[k * n + j] = a[p * n + j];
                a[p * n + j] = t;
            }
            int t = piv[k];
            piv[k] = piv[p];
            piv[p] = t;
        }
        /* 2. eliminate below the pivot */
        double d = a[k * n + k];
        for (int i = k + 1; i < n; i++) {
            double m = a[i * n + k] / d;
            a[i * n + k] = m;                 /* store the multiplier (L) */
            if (m == 0.0)
                continue;
            for (int j = k + 1; j < n; j++)
                a[i * n + j] -= m * a[k * n + j];
        }
    }
    return 0;
}

void lu_solve(const double *a, int n, const int *piv, double *b)
{
    /* apply the row swaps to b */
    double *tmp = malloc(sizeof(double) * (size_t)(n > 0 ? n : 1));
    for (int i = 0; i < n; i++)
        tmp[i] = b[piv[i]];
    /* forward substitution: L y = P b  (L has ones on the diagonal) */
    for (int i = 0; i < n; i++) {
        double s = tmp[i];
        for (int j = 0; j < i; j++)
            s -= a[i * n + j] * tmp[j];
        tmp[i] = s;
    }
    /* back substitution: U x = y */
    for (int i = n - 1; i >= 0; i--) {
        double s = tmp[i];
        for (int j = i + 1; j < n; j++)
            s -= a[i * n + j] * tmp[j];
        tmp[i] = s / a[i * n + i];
    }
    for (int i = 0; i < n; i++)
        b[i] = tmp[i];
    free(tmp);
}
