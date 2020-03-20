#include <iostream>

#include "../../include/propagator/propagator-on-box-1d.h"
#include "../../include/tridiag.h"
#include "../../include/array.h"
#include "../../include/wf/wavefunction-on-box-1d.h"


Propagator_on_Box_1D::Propagator_on_Box_1D() {}


Propagator_on_Box_1D::Propagator_on_Box_1D(
		size_t Nx, double dx, double *Vx, double hbar, double mass): 
	Nx(Nx), dx(dx), Vx(Vx), hbar(hbar), mass(mass) {

	if (dx <= 0.) { throw "[ERROR] Negative `dx` found"; }

	N_tridiag = 3*Nx;

	U_forward = new std::complex<double>[N_tridiag];
	U_backward = new std::complex<double>[N_tridiag];

	M2 = new double[N_tridiag];
	eval_M2_tridiag(M2, Nx);

	// Construct `M2Hreal`
	M2ReH = new double[N_tridiag];

	double *M2ReVx = new double[N_tridiag];
	tridiag_mul_diag(M2, Vx, M2ReVx, Nx);

	double *Kx = new double[N_tridiag];  // == -hbar*hbar/(2*mass) * D2
	eval_D2_tridiag(Kx, dx, Nx, -hbar*hbar/(2.*mass));

	add(Kx, M2ReVx, M2ReH, N_tridiag);

	delete [] Kx;
	delete [] M2ReVx;
}


Propagator_on_Box_1D::~Propagator_on_Box_1D() {
	delete [] U_forward;
	delete [] U_backward;
	delete [] M2;
	delete [] M2ReH;
}


int Propagator_on_Box_1D::eval_time_evol_unitary_for_real_timestep(double dt) {
	
	// Construct unitary propagators
	const double _c = -0.5*dt/hbar;
	std::complex<double> *pUf, *pUf_max, *pUb;
	double *pM2, *pM2ReH;
	double _M2, _c_M2ReH;
	for (pUf=U_forward, pUf_max=U_forward+3*Nx, pUb=U_backward, 
			pM2=M2, pM2ReH=M2ReH;
			pUf < pUf_max;
			++pUf, ++pUb, ++pM2, ++pM2ReH) {
		_M2 = *pM2, _c_M2ReH = _c * (*pM2ReH);
		*pUf = {_M2, _c_M2ReH};
		*pUb = {_M2, -_c_M2ReH};
	}

	// Set both ends of each time evol unitary tridiagonals
	U_forward[3*Nx-1] = 0.0, U_backward[3*Nx-1] = 0.0;
	U_forward[0] = 0.0, U_backward[0] = 0.0;

	return EXIT_SUCCESS;
}


int Propagator_on_Box_1D::eval_time_evol_unitary_for_imag_timestep(
		double dt_imag) {

	// Construct unitary propagators
	const double _c = -0.5*dt_imag/hbar;
	std::complex<double> *pUf, *pUb;
	double *pM2, *pM2max, *pM2ReH;
	double _M2, _c_M2ReH;
	for (pUf=U_forward, pUb=U_backward, pM2max=M2+3*Nx,
			pM2=M2, pM2ReH=M2ReH;
			pM2 < pM2max;
			++pUf, ++pUb, ++pM2, ++pM2ReH) {
		_M2 = *pM2, _c_M2ReH = _c * (*pM2ReH);
		// since these are pure real numbers, may optimized further
		*pUf = {_M2 + _c_M2ReH, 0.};  
		*pUb = {_M2 - _c_M2ReH, 0.};
	}
	
	// Set both ends of each time evol unitary tridiagonals
	U_forward[3*Nx-1] = 0.0, U_backward[3*Nx-1] = 0.0;
	U_forward[0] = 0.0, U_backward[0] = 0.0;

	return EXIT_SUCCESS;
}


int Propagator_on_Box_1D::propagate(
		std::complex<double> *wf, double dt, size_t Nt) 
{
	// Construct unitary propagators
	eval_time_evol_unitary_for_real_timestep(dt);

	// Propagate
	std::complex<double> *wf_mid = new std::complex<double>[Nx];
	for (size_t it=0; it<Nt; ++it) {
		tridiag_forward(U_forward, wf, wf_mid, Nx);
		tridiag_backward(U_backward, wf, wf_mid, Nx);
	}
	delete [] wf_mid;

	return EXIT_SUCCESS;
}


int Propagator_on_Box_1D::propagate_to_ground_state(
		std::complex<double> *wf, double dt, size_t Nt_max, double stop_thres) {
	
	std::complex<double> *wf_max = wf + Nx;
	std::complex<double> *wf_prev = new std::complex<double>[Nx];
	Wavefunction_on_Box_1D::normalize(wf, Nx, dx);	
	std::copy(wf, wf_max, wf_prev);

	double norm_diff;

	std::complex<double> *wf_diff = new std::complex<double>[Nx];
	
	// Propagate
	eval_time_evol_unitary_for_imag_timestep(dt);
	std::complex<double> *wf_mid = new std::complex<double>[Nx];
	size_t it;
	for (it=0; it<Nt_max; ++it) {

		tridiag_forward(U_forward, wf, wf_mid, Nx);
		tridiag_backward(U_backward, wf, wf_mid, Nx);

		Wavefunction_on_Box_1D::normalize(wf, Nx, dx);	

		substract(wf, wf_prev, wf_diff, Nx);
		norm_diff = Wavefunction_on_Box_1D::norm_sq(wf_diff, Nx, dx);
		if (norm_diff < stop_thres) { break; }

		std::copy(wf, wf_max, wf_prev);
	}
	int return_code = EXIT_SUCCESS;
	if (it >= Nt_max) { 
		std::cerr << "[ERROR] Maximum iteration (=" << Nt_max << ") exceeded\n"; 
		// should not return from here 
		// since we need to proceed to complete resource release process
		return_code = EXIT_FAILURE;  
	}

	delete [] wf_mid;
	delete [] wf_prev;
	delete [] wf_diff;

	return return_code;
}



