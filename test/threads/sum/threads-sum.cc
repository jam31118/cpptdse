#include <iostream>

#include <pthread.h>


struct psum_args {
	int *arr;
	size_t start_index;
	size_t N;
	int *psum;
};


void *partial_sum(void *args) {
	struct psum_args *pargs = (struct psum_args *) args;
	int *const pastart = pargs->arr + pargs->start_index;
	int *const pamax = pastart + pargs->N;
	int sum = 0;
	for (int *pa = pastart; pa < pamax; pa++) { sum += *pa; }
	*(pargs->psum) = sum;
	return NULL;
}


int main() {

	const size_t Ns = 10;
	int *seq = new int[Ns];	
	for (size_t i=0; i<Ns; ++i) { seq[i] = i+1; }

	const size_t Nths = 3;
	pthread_t *threads = new pthread_t[Nths];
	
	int *psum_arr = new int[Nths];

	const size_t min_N = Ns / Nths, remainder_N = Ns % Nths;
	struct psum_args *pargs_arr = new struct psum_args[Nths];
	for (size_t iths=0, istart=0; iths<Nths; iths++) {
		size_t my_N = min_N + (iths < remainder_N);
		pargs_arr[iths] = {seq, istart, my_N, psum_arr + iths};
		pthread_create(threads+iths, NULL, &partial_sum, (void *) (pargs_arr+iths));
		istart += my_N;
	}
	delete[] pargs_arr;

	for (size_t iths=0; iths<Nths; iths++) {
		pthread_join(threads[iths], NULL);
	}

	int total_sum = 0;
	for (size_t iths=0; iths<Nths; iths++) {
		int psum = psum_arr[iths];
		fprintf(stdout, "[thread=%02ld] partial sum = %d\n", iths, psum);
		total_sum += psum;
	}	fprintf(stdout, "Total sum = %d\n", total_sum); 

	delete [] psum_arr;
	delete [] threads;

	delete [] seq;
	
	return EXIT_SUCCESS;
}
