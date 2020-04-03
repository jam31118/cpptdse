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

/**
 * The number of grid points `Nx` is the total number except two at both ends.
 */
double Wavefunction_on_Box_1D::mean_x(
		std::complex<double> *wf, size_t Nx, double dx, double xmin) {
	
	double _mean_x = 0., _x = xmin;
	std::complex<double> _wf;
	for (std::complex<double> *pwf=wf, *pwfmax=wf+Nx; pwf<pwfmax; ++pwf)
	{ 
		_wf = *pwf; 
		_x += dx;
		_mean_x += _x * (std::conj(_wf) * _wf).real();
	} _mean_x *= dx;
	return _mean_x;
}


double Wavefunction_on_Box_1D::mean_x(std::complex<double> *wf, double xmin) {
	return Wavefunction_on_Box_1D::mean_x(wf, Nx, dx, xmin);
}


int Wavefunction_on_Box_1D::mean_and_stdev_x(
		std::complex<double> *wf, size_t Nx, double dx, double xmin, 
		double *p_mean_x, double *p_stdev_x) {
	
	double _mean_x = 0., _mean_x_sq=0., _x = xmin;
	double _piece;
	std::complex<double> _wf;
	for (std::complex<double> *pwf=wf, *pwfmax=wf+Nx; pwf<pwfmax; ++pwf)
	{ 
		_wf = *pwf; 
		_x += dx;
		_piece = _x * (std::conj(_wf) * _wf).real();
		_mean_x += _piece;
		_piece *= _x;
		_mean_x_sq += _piece;
	} _mean_x *= dx; _mean_x_sq *= dx;
	double _var_x = _mean_x_sq - _mean_x*_mean_x;
	double _stdev_x = std::sqrt(_var_x);
	*p_stdev_x = _stdev_x; *p_mean_x = _mean_x;

	return EXIT_SUCCESS;
}



int Wavefunction_on_Box_1D::mean_and_stdev_x(
		std::complex<double> *wf, double xmin, double *p_mean_x, double *p_stdev_x)
{
	return Wavefunction_on_Box_1D::mean_and_stdev_x(
			wf, Nx, dx, xmin, p_mean_x, p_stdev_x);
}



int Wavefunction_on_Box_1D::stdev_x(
		std::complex<double> *wf, size_t Nx, double dx, double xmin, 
		double *p_stdev_x) 
{
	double _mean_x, _stdev_x;
	int status = Wavefunction_on_Box_1D::mean_and_stdev_x(
			wf, Nx, dx, xmin, &_mean_x, &_stdev_x);
	if (status != EXIT_SUCCESS) { return EXIT_FAILURE; }
	*p_stdev_x = _stdev_x;
	return EXIT_SUCCESS;
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

size_t Wavefunction_on_Box_1D::get_Nx_tot() 
{ return 1 + Nx + 1; }

double Wavefunction_on_Box_1D::get_xmax(double xmin)
{ return xmin + (Nx+1)*dx; }

double Wavefunction_on_Box_1D::get_dx() { return dx; }

