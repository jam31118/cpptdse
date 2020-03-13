#ifndef _ARRAY_HH_
#define _ARRAY_HH_

#include <cstdlib>
#include <iostream>

int add(double *a, double *b, double *c, size_t N);

template <typename T>
void print_array(T *a, size_t N) {
	for (T *pT=a, *pTmax=a+N; pT < pTmax; ++pT) {
		std::cout << *pT;
	} std::cout << std::endl;
}

template <typename T1, typename T2>
int array_mul_scalar(T1 *a, size_t N, T2 c) {
	for (T1 *pa=a, *pamax=a+N; pa<pamax; ++pa) { *pa *= c; }
	return EXIT_SUCCESS;
}

#endif // _ARRAY_HH_
