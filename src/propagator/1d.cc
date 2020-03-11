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

	U_forward = new std::complex<double>[3*Nx];
	U_backward = new std::complex<double>[3*Nx];

	M2 = new double[3*Nx];
	eval_M2_tridiag(M2, Nx);

	// Construct `M2Hreal`
	M2ReH = new double[3*Nx];

	double *M2ReVx = new double[3*Nx];
	tridiag_mul_diag(M2, Vx, M2ReVx, Nx);

	double *Kx = new double[3*Nx];  // == -hbar*hbar/(2*mass) * D2
	eval_D2_tridiag(Kx, dx, Nx, -hbar*hbar/(2.*mass));

	add(Kx, M2ReVx, M2ReH, Nx);

	delete [] Kx;
	delete [] M2ReVx;
}

Propagator_on_Box_1D::~Propagator_on_Box_1D() {
	delete [] U_forward;
	delete [] U_backward;
	delete [] M2;
	delete [] M2ReH;
}

int Propagator_on_Box_1D::propagate(
		std::complex<double> *wf, double dt, size_t Nt) 
{
	
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

	// Propagate
	std::complex<double> *wf_mid = new std::complex<double>[Nx];
	for (size_t it=0; it<Nt; ++it) {
		tridiag_forward(U_forward, wf, wf_mid, Nx);
		tridiag_backward(U_backward, wf, wf_mid, Nx);
	}
	delete [] wf_mid;

	return EXIT_SUCCESS;
}


