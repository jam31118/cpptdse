#ifndef _WAVEFUNCTION_FOR_SPHARM_ZERO_M_H_
#define _WAVEFUNCTION_FOR_SPHARM_ZERO_M_H_

#include <complex>

class Wavefunction_for_spharm_zero_m {

	size_t Nr;
	double dr;
	size_t Nl;

public:
	Wavefunction_for_spharm_zero_m(size_t Nr, double dr, size_t Nl);
	static double norm_sq(
			std::complex<double> *wf, size_t Nr, double dr, size_t Nl);
	double norm_sq(std::complex<double> *wf);
	static int normalize(
			std::complex<double> *wf, size_t Nr, double dr, size_t Nl);
	int normalize(std::complex<double> *wf);
	static size_t length(size_t Nr, size_t Nl);
	size_t length();
	static int eval_radial_arr(
			double *const rarr, const size_t Nr, const double dr);
	int eval_radial_arr(double *const rarr);
};

#endif // _WAVEFUNCTION_FOR_SPHARM_ZERO_M_H_
