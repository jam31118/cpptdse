#include <iostream>

#include "../../../include/array-threaded.h"

int main() {

	const size_t Ns = 100;
	double *seq = new double[Ns];
	for (size_t i=0; i<Ns; ++i) { seq[i] = i+1; }

	double total_sum = 0;
	const size_t Nths = 3;
	sum_th(seq, Ns, &total_sum, Nths);
	delete [] seq;

	fprintf(stdout, "Total sum = %f\n", total_sum);
	
	return EXIT_SUCCESS;
}

