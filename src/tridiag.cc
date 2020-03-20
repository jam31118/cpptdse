#include "../include/tridiag.h"
#include <complex>


int eval_M2_tridiag(double *T, size_t N) {

	const double val_diag = 10./12., val_offdiag = 1./12.;

	double *pT = T+1, *pT_max = T;
	for (pT_max += N; pT < pT_max; ++pT) { *pT = val_offdiag; }
	for (pT_max += N; pT < pT_max; ++pT) { *pT = val_diag; }
	for (pT_max += N-1; pT < pT_max; ++pT) { *pT = val_offdiag; }

	return EXIT_SUCCESS;
}


int eval_M1_tridiag(double *T, size_t N) {

	const double val_diag = 2./3., val_offdiag = 1./6.;

	double *pT = T+1, *pT_max = T;
	for (pT_max += N; pT < pT_max; ++pT) { *pT = val_offdiag; }
	for (pT_max += N; pT < pT_max; ++pT) { *pT = val_diag; }
	for (pT_max += N-1; pT < pT_max; ++pT) { *pT = val_offdiag; }

	return EXIT_SUCCESS;
}


int eval_D2_tridiag(double *T, double h, size_t N, double coef) {

	if (!(h > 0)) { return EXIT_FAILURE; }
	
	const double val_offdiag = coef / (h*h);
	const double val_diag = -2. * val_offdiag;

	double *pT = T+1, *pT_max = T;
	for (pT_max += N; pT < pT_max; ++pT) { *pT = val_offdiag; }
	for (pT_max += N; pT < pT_max; ++pT) { *pT = val_diag; }
	for (pT_max += N-1; pT < pT_max; ++pT) { *pT = val_offdiag; }

	return EXIT_SUCCESS;
}


int eval_D1_tridiag(double *T, double h, size_t N) {

	if (!(h > 0)) { return EXIT_FAILURE; }
	
	const double val_lower_offdiag = - 1. / (2.*h);
	const double val_upper_offdiag = 1. / (2.*h);
	const double val_diag = 0.;

	double *pT = T+1, *pT_max = T;
	for (pT_max += N; pT < pT_max; ++pT) { *pT = val_lower_offdiag; }
	for (pT_max += N; pT < pT_max; ++pT) { *pT = val_diag; }
	for (pT_max += N-1; pT < pT_max; ++pT) { *pT = val_upper_offdiag; }

	return EXIT_SUCCESS;
}


int tridiag_mul_diag(double *T, double *D, double *TD, size_t N) {

	double *pT = T+1, *pTD = TD+1, *pD, *pD_max;
	for (pD=D,pD_max=D+N-1; pD<pD_max; ++pD,++pT,++pTD) { *pTD = (*pT) * (*pD); }
	for (pD=D,pD_max=D+N; pD<pD_max; ++pD,++pT,++pTD) { *pTD = (*pT) * (*pD); }
	for (pD=D+1,pD_max=D+N; pD<pD_max; ++pD,++pT,++pTD) { *pTD = (*pT) * (*pD); }

	return EXIT_SUCCESS;
}


void print_tridiag(std::complex<double> *td, size_t N) {
	std::complex<double> *ptd=td, *ptd_max=td;
	for (int i=0; i<3; ++i) {
		for (ptd_max+=N; ptd<ptd_max; ++ptd) {
			std::cout << *ptd;
		} std::cout << std::endl;
	}
}

