#ifndef _TRIDIAG_H_
#define _TRIDIAG_H_

#include <cstdlib>
#include <complex>

#include "matrix.hh"

int eval_M2_tridiag(double *T, size_t N);
int eval_D2_tridiag(double *T, double h, size_t N, double coef=1.0);

int eval_M1_tridiag(double *T, size_t N);
int eval_D1_tridiag(double *T, double h, size_t N);

int tridiag_mul_diag(double *T, double *D, double *TD, size_t N);

template <typename TD, typename type_vec>
int tridiag_forward(TD *td, type_vec *v, type_vec *b, size_t N)
{ return tridiag_mul_forward(td, td+N, td+2*N, v, b, N); }

template <typename TD, typename type_vec>
int tridiag_backward(TD *td, type_vec *v, type_vec *b, size_t N) 
{ return tridiag_mul_backward(td, td+N, td+2*N, v, b, N); }

void print_tridiag(std::complex<double> *td, size_t N);

#endif // _TRIDIAG_H_
