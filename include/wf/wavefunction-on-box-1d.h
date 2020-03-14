#ifndef _WAVEFUNCTION_ON_BOX_1D_H_
#define _WAVEFUNCTION_ON_BOX_1D_H_

#include <complex>


/*
 * Members of `Wavefunction_on_Box_1D`
 * */

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



/*
 * A set of analytical expressions for possible states in Box_1D
 * */

int eval_ground_state_in_box_1d(
		std::complex<double> *wf, size_t Nx, double dx);


#endif // _WAVEFUNCTION_ON_BOX_1D_H_
