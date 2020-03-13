#include "../../include/wf/wavefunction-on-box-1d.h"
#include "../../include/array.hh"


Wavefunction_on_Box_1D::Wavefunction_on_Box_1D(size_t Nx, double dx): 
	Nx(Nx), dx(dx) {
//	std::cerr << "in:Wavefunction_on_Box_1D: first\n";
	if (!(Nx > 0)) { throw "`Nx` should be positive"; }
	if (!(dx > 0)) { throw "`dx` should be real"; }
//	std::cerr << "in:Wavefunction_on_Box_1D: end\n";
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

