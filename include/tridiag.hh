#ifndef _TRIDIAG_HH_
#define _TRIDIAG_HH_

#include <cstdlib>

int eval_M2_tridiag(double *T, size_t N);
int eval_D2_tridiag(double *T, double h, size_t N);
int tridiag_mul_diag(double *T, double *D, double *TD, size_t N);

#endif // _TRIDIAG_HH_
