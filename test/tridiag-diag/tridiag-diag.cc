#include <iostream>
#include <cstdlib>

#include "../../include/tridiag.hh"

int main() {

	const size_t N = 7;
	double *M2 = new double[3*N];
	eval_M2_tridiag(M2, N);
	
	for (size_t i=0; i<3*N; i++) {
		std::cout << M2[i] << " ";
	} std::cout << std::endl;

	delete [] M2;

	return EXIT_SUCCESS;
}

