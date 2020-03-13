#include "../include/array.h"

int add(double *a, double *b, double *c, size_t N) {
	for (double *pa=a, *pb=b, *pc=c, *pcmax=c+N; pc < pcmax; ++pa, ++pb, ++pc) {
		*pc = *pa + *pb;
	}
	return EXIT_SUCCESS;
}
