#include "../../include/wf/wavefunction-for-spharm-zero-m.h"

#include "../../include/array.h"

Wavefunction_for_spharm_zero_m
::Wavefunction_for_spharm_zero_m(size_t Nr, double dr, size_t Nl): 
	Nr(Nr), dr(dr), Nl(Nl) {}

double Wavefunction_for_spharm_zero_m::norm_sq(
		std::complex<double> *wf, size_t Nr, double dr, size_t Nl) {
	
	const size_t N_total = Nr * Nl;
	std::complex<double> *const pwfmax = wf + N_total;
	double _norm_sq = 0.;
	std::complex<double> _wf;
	for (std::complex<double> *pwf=wf; pwf < pwfmax; ++pwf)
	{ 
		_wf = *pwf;
		_norm_sq += (std::conj(_wf) * _wf ).real(); 
	} _norm_sq *= dr;
	return _norm_sq;

}


double Wavefunction_for_spharm_zero_m::norm_sq(std::complex<double> *wf)
{ return this->norm_sq(wf, Nr, dr, Nl); }


int Wavefunction_for_spharm_zero_m::normalize(
		std::complex<double> *wf, size_t Nr, double dr, size_t Nl) {
	double _norm_sq = Wavefunction_for_spharm_zero_m::norm_sq(wf, Nr, dr, Nl);
	const size_t N_total = Nr * Nl;
	return array_mul_scalar(wf, N_total, 1./std::sqrt(_norm_sq));	
}


int Wavefunction_for_spharm_zero_m::normalize(std::complex<double> *wf)
{ return this->normalize(wf, Nr, dr, Nl); }


size_t Wavefunction_for_spharm_zero_m
::length(size_t Nr, size_t Nl) { return Nr * Nl; }

size_t Wavefunction_for_spharm_zero_m
::length() { return length(Nr, Nl); }


int Wavefunction_for_spharm_zero_m::eval_radial_arr(
		double *const rarr, const size_t Nr, const double dr) {

	for (double *pa=rarr, *const pamax=rarr+Nr, r=dr; pa<pamax; ++pa, r+=dr) 
	{ *pa = r; }

	return EXIT_SUCCESS;
}


int Wavefunction_for_spharm_zero_m::eval_radial_arr(double *const rarr) {
	return eval_radial_arr(rarr, Nr, dr);
}

