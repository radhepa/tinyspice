/*
 * test_core.c - unit tests for the building blocks.
 * Build and run with:  make test
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "devices.h"
#include "netlist.h"
#include "solver.h"

static int failures = 0, checks = 0;

#define CHECK(cond, ...)                                  \
    do {                                                  \
        checks++;                                         \
        if (!(cond)) {                                    \
            failures++;                                   \
            printf("FAIL %s:%d: ", __FILE__, __LINE__);   \
            printf(__VA_ARGS__);                          \
            printf("\n");                                 \
        }                                                 \
    } while (0)

static int close_to(double a, double b, double rel)
{
    return fabs(a - b) <= rel * fmax(fabs(a), fabs(b)) + 1e-300;
}

static void test_lu_small(void)
{
    /* 2x + y - z = 8, -3x - y + 2z = -11, -2x + y + 2z = -3  ->  x=2, y=3, z=-1 */
    double a[9] = {2, 1, -1, -3, -1, 2, -2, 1, 2};
    double b[3] = {8, -11, -3};
    int piv[3];
    CHECK(lu_factor(a, 3, piv) == 0, "3x3 factor failed");
    lu_solve(a, 3, piv, b);
    CHECK(close_to(b[0], 2, 1e-12) && close_to(b[1], 3, 1e-12) && close_to(b[2], -1, 1e-12),
          "3x3 solution wrong: %g %g %g", b[0], b[1], b[2]);
}

static void test_lu_random(void)
{
    const int n = 60;
    double *a = malloc(sizeof(double) * n * n), *a0 = malloc(sizeof(double) * n * n);
    double *x = malloc(sizeof(double) * n), *b = malloc(sizeof(double) * n);
    int *piv = malloc(sizeof(int) * n);
    srand(42);
    for (int i = 0; i < n * n; i++)
        a0[i] = a[i] = (double)rand() / RAND_MAX - 0.5;
    for (int i = 0; i < n; i++) {
        x[i] = (double)rand() / RAND_MAX;
        b[i] = 0;
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            b[i] += a0[i * n + j] * x[j];
    CHECK(lu_factor(a, n, piv) == 0, "random 60x60 reported singular");
    lu_solve(a, n, piv, b);
    double worst = 0;
    for (int i = 0; i < n; i++)
        worst = fmax(worst, fabs(b[i] - x[i]));
    CHECK(worst < 1e-9, "random 60x60 error %g", worst);
    free(a); free(a0); free(x); free(b); free(piv);
}

static void test_lu_singular(void)
{
    double a[4] = {1, 2, 2, 4};          /* row 2 = 2 x row 1 */
    int piv[2];
    CHECK(lu_factor(a, 2, piv) != 0, "singular matrix not detected");
}

static void test_values(void)
{
    double v;
    CHECK(!parse_value("1k", &v) && close_to(v, 1e3, 1e-12), "1k -> %g", v);
    CHECK(!parse_value("4.7u", &v) && close_to(v, 4.7e-6, 1e-12), "4.7u -> %g", v);
    CHECK(!parse_value("2meg", &v) && close_to(v, 2e6, 1e-12), "2meg -> %g", v);
    CHECK(!parse_value("10pF", &v) && close_to(v, 1e-11, 1e-12), "10pF -> %g", v);
    /* the classic SPICE trap: M means milli, not mega */
    CHECK(!parse_value("1M", &v) && close_to(v, 1e-3, 1e-12), "1M -> %g", v);
    CHECK(!parse_value("1e-3", &v) && close_to(v, 1e-3, 1e-12), "1e-3 -> %g", v);
    CHECK(parse_value("abc", &v) != 0, "abc should fail");
}

static void test_diode_derivative(void)
{
    Model m = {"d", MOD_D, 1e-14, 1.0, 0, 0, 0};
    for (double vd = -1.0; vd <= 0.8; vd += 0.1) {
        double i1, g1, i2, g2, id, gd, h = 1e-6;
        diode_eval(&m, vd, &id, &gd);
        diode_eval(&m, vd + h, &i1, &g1);
        diode_eval(&m, vd - h, &i2, &g2);
        double numeric = (i1 - i2) / (2 * h);
        CHECK(close_to(gd, numeric, 1e-5) || fabs(gd - numeric) < 1e-15,
              "diode gd %g vs numeric %g at %g V", gd, numeric, vd);
    }
}

static void test_pulse(void)
{
    Source s = {SRC_PULSE, 0, {0, 5, 1e-6, 1e-9, 1e-9, 4e-6, 10e-6}, 0, NULL, NULL};
    CHECK(source_value(&s, 0.5e-6, 1e-9, 1e-3) == 0, "before delay");
    CHECK(close_to(source_value(&s, 1e-6 + 0.5e-9, 1e-9, 1e-3), 2.5, 1e-9), "mid-rise");
    CHECK(source_value(&s, 3e-6, 1e-9, 1e-3) == 5, "high");
    CHECK(source_value(&s, 8e-6, 1e-9, 1e-3) == 0, "low");
    CHECK(source_value(&s, 13e-6, 1e-9, 1e-3) == 5, "second period high");
}

int main(void)
{
    test_lu_small();
    test_lu_random();
    test_lu_singular();
    test_values();
    test_diode_derivative();
    test_pulse();
    printf("%d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}
