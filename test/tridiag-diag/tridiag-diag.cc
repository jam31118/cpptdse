#include <iostream>
#include <cstdlib>

#include "../../include/tridiag.hh"

int main() {

	const size_t N = 7;

	double *M2 = new double[3*N];
	eval_M2_tridiag(M2, N);

	double *V = new double[N];
	for (double *pV=V, *pVmax=V+N; pV<pVmax; ++pV) { *pV = 12.; }
	V[0] = 24.;

	double *M2V = new double[3*N];
	tridiag_mul_diag(M2, V, M2V, N);

	std::cout << "M2:\n";
	for (int irow=0; irow<3; ++irow) {
		for (size_t i=0; i<N; ++i) {
			std::cout << M2[irow*N+i] << " ";
		} std::cout << std::endl;
	}

	std::cout << "V:\n";
	for (size_t i=0; i<N; ++i) { std::cout << V[i] << " "; }
	std::cout << std::endl;

	std::cout << "M2V:\n";
	for (int irow=0; irow<3; ++irow) {
		for (size_t i=0; i<N; ++i) {
			std::cout << M2V[irow*N+i] << " ";
		} std::cout << std::endl;
	}
	

	delete [] M2;
	delete [] V;
	delete [] M2V;

	return EXIT_SUCCESS;
}

