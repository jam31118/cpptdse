#ifndef _PROPAGATOR_HH_
#define _PROPAGATOR_HH_

#include <cstdlib>
#include <complex>

#include "../wf/wavefunction-on-box-1d.h"

class Propagator_on_Box_1D {

	size_t Nx;
	double dx, *Vx, hbar, mass;
	std::complex<double> *U_forward, *U_backward;
	double *M2;
	double *M2ReH;

public:

	Wavefunction_on_Box_1D *wf = NULL;

	Propagator_on_Box_1D(
			size_t Nx, double dx, double *Vx, double hbar=1, double mass=1);
	~Propagator_on_Box_1D();
	int eval_time_evol_unitary_for_real_timestep(double dt);
	int eval_time_evol_unitary_for_imag_timestep(double dt_imag);
	int propagate(std::complex<double> *wf, double dt, size_t Nt);

};


// [NOTE] In order to use some of the members of parent class,
// those members should be declared as `private` in the parent class.
class Propagator_on_Box_1D_with_imag_pot {
	Propagator_on_Box_1D_with_imag_pot(
			size_t Nx, double dx, double *Vx, double *ImVx, double hbar=1, double mass=1);
	~Propagator_on_Box_1D_with_imag_pot();
	int eval_time_evol_unitary_for_real_timestep(double dt);
	int eval_time_evol_unitary_for_imag_timestep(double dt_imag);
};

#endif // _PROPAGATOR_HH_
