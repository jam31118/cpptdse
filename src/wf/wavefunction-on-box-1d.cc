#include "../../include/wf/wavefunction-on-box-1d.h"
#include "../../include/array.h"

#include <cmath>


Wavefunction_on_Box_1D::Wavefunction_on_Box_1D(size_t Nx, double dx): 
	Nx(Nx), dx(dx) {
	if (!(Nx > 0)) { throw "`Nx` should be positive"; }
	if (!(dx > 0)) { throw "`dx` should be real"; }
};


double Wavefunction_on_Box_1D::norm_sq(std::complex<double> *wf) {
	return this->norm_sq(wf, this->Nx, this->dx);
}


double Wavefunction_on_Box_1D::norm_sq(
		std::complex<double> *wf, size_t Nx, double dx) {
	double _norm_sq = 0.0;
	std::complex<double> _wf;
	for (std::complex<double> *pwf=wf, *pwfmax=wf+Nx; pwf<pwfmax; ++pwf)
	{ 
		_wf = *pwf;
		_norm_sq += std::real( std::conj(_wf) * _wf ); 
	} _norm_sq *= dx;
	return _norm_sq;
}


int Wavefunction_on_Box_1D::normalize(
		std::complex<double> *wf, size_t Nx, double dx) {
	double _norm_sq = Wavefunction_on_Box_1D::norm_sq(wf, Nx, dx);
	return array_mul_scalar(wf, Nx, 1./std::sqrt(_norm_sq));	
}

int Wavefunction_on_Box_1D::normalize(std::complex<double> *wf) {
	return this->normalize(wf, this->Nx, this->dx);
}



int eval_ground_state_in_box_1d(
		std::complex<double> *wf, size_t Nx, double dx) {
	const size_t Nx_tot = 1+Nx+1;
	const double L = dx * (Nx_tot-1);
	const double c = std::sqrt(2./L);  // a normalizing constant
	for (size_t i=0; i<Nx; ++i) { wf[i] = c * std::sin(M_PI / L * (i+1)*dx); }
	return EXIT_SUCCESS;	
}


// Construct spatial array
//
int Wavefunction_on_Box_1D::eval_x_tot_arr(double *const x_tot_arr, 
		const size_t Nx_tot, const double dx, const double xmin) {
	
  for (double *px=x_tot_arr, *pxmax=x_tot_arr+Nx_tot, val=xmin; 
			px < pxmax; ++px, val+=dx)
  { *px = val; }

	return EXIT_SUCCESS;
}

