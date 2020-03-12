#include <iostream>

#include "../../include/propagator/1d.hh"
#include "../../include/tridiag.hh"
#include "../../include/array.hh"
#include "../../include/tridiag.hh"

#include "matrix.hh"

Propagator_on_Box_1D::Propagator_on_Box_1D(
		size_t Nx, double dx, double *Vx, double hbar, double mass): 
	Nx(Nx), dx(dx), Vx(Vx), hbar(hbar), mass(mass) {

	if (dx <= 0.) { throw "[ERROR] Negative `dx` found"; }

	const size_t N_tridiag = 3*Nx;

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
//	U_forward[3*Nx-1] = 0.0, U_backward[3*Nx-1] = 0.0;
//	U_forward[0] = 0.0, U_backward[0] = 0.0;
	return EXIT_SUCCESS;
}

int Propagator_on_Box_1D::eval_time_evol_unitary_for_imag_timestep(double dt_imag) {
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
	return EXIT_SUCCESS;
}

int Propagator_on_Box_1D::propagate(
		std::complex<double> *wf, double dt, size_t Nt) 
{

	eval_time_evol_unitary_for_real_timestep(dt);

	std::cout << "U_forward: \n";
	print_tridiag(U_forward, Nx);
	std::cout << "U_backward: \n";
	print_tridiag(U_backward, Nx);
	
	// Construct unitary propagators
//	const double _c = -0.5*dt/hbar;
//	std::complex<double> *pUf, *pUf_max, *pUb;
//	double *pM2, *pM2ReH;
//	double _M2, _c_M2ReH;
//	for (pUf=U_forward, pUf_max=U_forward+3*Nx, pUb=U_backward, 
//			pM2=M2, pM2ReH=M2ReH;
//			pUf < pUf_max;
//			++pUf, ++pUb, ++pM2, ++pM2ReH) {
//		_M2 = *pM2, _c_M2ReH = _c * (*pM2ReH);
//		*pUf = {_M2, _c_M2ReH};
//		*pUb = {_M2, -_c_M2ReH};
//	}

	// Propagate
	std::complex<double> *wf_mid = new std::complex<double>[Nx];
	for (size_t it=0; it<Nt; ++it) {
		tridiag_forward(U_forward, wf, wf_mid, Nx);
		tridiag_backward(U_backward, wf, wf_mid, Nx);
	}
	delete [] wf_mid;

	return EXIT_SUCCESS;
}


