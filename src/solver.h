/*
 * solver.h - dense LU decomposition with partial pivoting.
 *
 * Solves A x = b. A is stored row-major: A[i*n + j] is row i, column j.
 * Circuits here have at most a few hundred unknowns, so a dense O(n^3)
 * solver is simple and fast enough.
 */
#ifndef SOLVER_H
#define SOLVER_H

/* Factor A in place into L and U (PA = LU). piv receives the row order.
 * Returns 0 on success, or 1 + the column index where A is singular. */
int lu_factor(double *a, int n, int *piv);

/* Solve using a factored matrix; b is overwritten with the solution x. */
void lu_solve(const double *a, int n, const int *piv, double *b);

#endif
