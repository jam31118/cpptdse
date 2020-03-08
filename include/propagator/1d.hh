#ifndef _PROPAGATOR_HH_
#define _PROPAGATOR_HH_

#include <cstdlib>
#include <complex>

//#include "../tridiag.hh"
//#include "../array.hh"

class Propagator_on_Box_1D {

	size_t Nx;
	double dx, *Vx, hbar, mass;
	std::complex<double> *U_forward, *U_backward;
	double *M2;
	double *M2ReH;

public:
	Propagator_on_Box_1D(
			size_t Nx, double dx, double *Vx, double hbar=1, double mass=1);
	~Propagator_on_Box_1D();
	int propagate(std::complex<double> *wf, double dt, size_t Nt);
};

#endif // _PROPAGATOR_HH_
