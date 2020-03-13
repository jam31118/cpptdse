#ifndef _WAVEFUNCTION_ON_BOX_1D_H_
#define _WAVEFUNCTION_ON_BOX_1D_H_

#include <complex>

class Wavefunction_on_Box_1D {
	size_t Nx;
	double dx;
public:
	Wavefunction_on_Box_1D(size_t Nx, double dx);
	static double norm_sq(std::complex<double> *wf, size_t Nx, double dx);
	double norm_sq(std::complex<double> *wf);
	static int normalize(std::complex<double> *wf, size_t Nx, double dx);
	int normalize(std::complex<double> *wf);
};

#endif // _WAVEFUNCTION_ON_BOX_1D_H_
