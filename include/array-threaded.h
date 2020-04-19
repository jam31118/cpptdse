#ifndef _ARRAY_THREADED_H_
#define _ARRAY_THREADED_H_

#include <iostream>

#include <pthread.h>


template<typename T>
struct psum_args {
	T *arr;
	size_t start_index;
	size_t N_partial;
	T *psum;
};


template<typename T>
void *partial_sum(void *args) {
	struct psum_args<T> *pargs = (struct psum_args<T> *) args;
	T *const pastart = pargs->arr + pargs->start_index;
	T *const pamax = pastart + pargs->N_partial;
	T partial_sum = 0;
	for (T *pa = pastart; pa < pamax; pa++) { partial_sum += *pa; }
	*(pargs->psum) = partial_sum;
	return NULL;
}


template<typename T>
int sum_th(T *const arr, const size_t N, T *sum, const size_t Nths) {

	pthread_t *threads = new pthread_t[Nths];
	
	T *psums = new T[Nths];

	const size_t min_N = N / Nths, remainder_N = N % Nths;
	struct psum_args<T> *pargs_arr = new struct psum_args<T>[Nths];
	for (size_t iths=0, istart=0; iths<Nths; iths++) {
		size_t my_N = min_N + (iths < remainder_N);
		pargs_arr[iths] = {arr, istart, my_N, psums + iths};
		pthread_create(threads+iths, NULL, &partial_sum<T>, (void *)(pargs_arr+iths));
		istart += my_N;
	}
	for (size_t iths=0; iths<Nths; iths++) { pthread_join(threads[iths], NULL); }
	delete [] pargs_arr;
	delete [] threads;

	T total_sum = 0;
	for (T *pa=psums,*pamax=psums+Nths; pa < pamax; ++pa) { total_sum += *pa; }
	delete [] psums;

	*sum = total_sum;

	return EXIT_SUCCESS;
}



#endif // _ARRAY_THREADED_H_
