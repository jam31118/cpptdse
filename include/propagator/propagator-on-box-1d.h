#ifndef _PROPAGATOR_ON_BOX_1D_H_
#define _PROPAGATOR_ON_BOX_1D_H_

#include <cstdlib>
#include <complex>

#include "../wf/wavefunction-on-box-1d.h"


class Propagator_on_Box_1D {

protected:
	size_t Nx, N_tridiag;
	double dx;
	const size_t Ndim = 1;
	std::complex<double> *U_forward, *U_backward;
	double *Vx;
	double *M2;
	double *M2ReH;
	double hbar, mass;

public:
	Propagator_on_Box_1D();
	Propagator_on_Box_1D(size_t Nx, double dx, double *Vx, 
			double hbar=1, double mass=1);
	virtual ~Propagator_on_Box_1D();
	virtual int eval_time_evol_unitary_for_real_timestep(double dt);
	int eval_time_evol_unitary_for_imag_timestep(double dt_imag);
	int propagate(std::complex<double> *wf, double dt, size_t Nt);
	int propagate_to_ground_state(
			std::complex<double> *wf, double dt, size_t Nt_max, double stop_thres);
};


class Propagator_on_Box_1D_with_imag_pot : public Propagator_on_Box_1D {

	double *Vx_imag;
	double *M2ImV;

public:
	Propagator_on_Box_1D_with_imag_pot(
			size_t Nx, double dx, double *Vx, double *Vx_imag, 
			double hbar=1, double mass=1);
	~Propagator_on_Box_1D_with_imag_pot();
	int eval_time_evol_unitary_for_real_timestep(double dt);
//	int eval_time_evol_unitary_for_imag_timestep(double dt_imag);
};



#endif // _PROPAGATOR_ON_BOX_1D_H_
