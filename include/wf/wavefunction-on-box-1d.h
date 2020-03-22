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
	static const size_t Ndim = 1;
	Wavefunction_on_Box_1D(size_t Nx, double dx);
	static double norm_sq(std::complex<double> *wf, size_t Nx, double dx);
	double norm_sq(std::complex<double> *wf);
	static int normalize(std::complex<double> *wf, size_t Nx, double dx);
	int normalize(std::complex<double> *wf);
	static int eval_x_tot_arr(double *const x_tot_arr, 
			const size_t Nx_tot, const double dx, const double xmin);
	size_t get_Nx_tot();
	double get_xmax(double xmin);
	double get_dx();
};



/*
 * A set of analytical expressions for possible states in Box_1D
 * */

int eval_ground_state_in_box_1d(
		std::complex<double> *wf, size_t Nx, double dx);


#endif // _WAVEFUNCTION_ON_BOX_1D_H_
