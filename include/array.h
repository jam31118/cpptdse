#ifndef _ARRAY_H_
#define _ARRAY_H_

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

template <typename T1, typename T2, typename T>
int substract(T1 *a1, T2 *a2, T *a, size_t N) {
	T1 *pa1=a1; T2 *pa2=a2; T *pa=a, *pamax=a+N;
	for (; pa<pamax; ++pa, ++pa1, ++pa2) 
	{ *pa = (*pa1) - (*pa2); }
	return EXIT_SUCCESS;
}

#endif // _ARRAY_H_
