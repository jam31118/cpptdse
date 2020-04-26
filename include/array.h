#ifndef _ARRAY_H_
#define _ARRAY_H_

#include <cstdlib>
#include <iostream>
#include <random>

int add(double *a, double *b, double *c, size_t N);

template <typename type_v, typename type_v1, typename type_c, typename type_v2>
int v1_add_c_mul_v2(type_v *v, 
		const type_v1 *v1, const type_c c, const type_v2 *v2, const size_t N) 
{
	const type_v1 *pv1=v1; 
	const type_v2 *pv2=v2; 
	type_v *pv=v, *const pvmax=v+N;
	for (; pv<pvmax; ++pv, ++pv1, ++pv2) { *pv = (*pv1) + c * (*pv2); }
	return EXIT_SUCCESS;
}

template <typename T>
void print_array(T *a, size_t N) {
	for (T *pT=a, *pTmax=a+N; pT < pTmax; ++pT) {
		std::cout << *pT << " ";
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


template <typename T>
int set_to_zeros(T *a, size_t N) {
	for (T *pa=a, *pa_max=a+N; pa<pa_max; ++pa) { *pa = 0.; }
	return EXIT_SUCCESS;
}


// [NOTE] For an array of complex numbers, 
// .. imaginary part is set to zero by default
template <typename T>
int set_to_randoms(
		T *const a, const size_t N, const double l=-1., const double u=1.) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> unidis(l, u);
	for (T *pa=a, *pamax=a+N; pa<pamax; ++pa)
	{ *pa = unidis(gen); } 
	return EXIT_SUCCESS;
}


template <typename T>
int set_2d_view_of_1d(T **arr2d, T *arr1d, size_t N0, size_t N1) {
	for (
			T **parr2d = arr2d, **parr2d_max = arr2d + N0, *parr1d = arr1d; 
			parr2d < parr2d_max; 
			++parr2d, parr1d += N1) 
	{ *parr2d = parr1d; }
	return EXIT_SUCCESS;
}


#endif // _ARRAY_H_
